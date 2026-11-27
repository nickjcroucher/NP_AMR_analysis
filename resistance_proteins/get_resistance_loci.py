#!/usr/bin/env python3

import argparse
import tarfile
import re
import pandas as pd
from Bio import SeqIO


def normalize(seq):
    """Normalize protein sequences for exact matching."""
    return str(seq).strip().upper().rstrip("*")


def read_csv_from_tar(tar_path, filename):
    """
    Read a CSV file directly from a .tar.gz archive
    without extracting it to disk.
    """

    with tarfile.open(tar_path, mode="r:gz") as tar:

        members = [
            m for m in tar.getmembers()
            if m.isfile() and m.name.endswith(filename)
        ]

        if not members:
            raise FileNotFoundError(
                f"{filename} not found in archive: {tar_path}"
            )

        if len(members) > 1:
            raise ValueError(
                f"Multiple instances of {filename} found "
                f"in archive: {tar_path}"
            )

        with tar.extractfile(members[0]) as f:
            df = pd.read_csv(f)

    return df


def main():

    # --------------------------------------------------
    # Command-line arguments
    # --------------------------------------------------

    parser = argparse.ArgumentParser(
        description=(
            "Map protein sequences to Panaroo gene clusters "
            "using exact sequence matching."
        )
    )

    parser.add_argument(
        "--gene-data",
        required=True,
        help="Path to gene_data.csv.tar.gz"
    )

    parser.add_argument(
        "--presence",
        required=True,
        help="Path to gene_presence_absence.csv.tar.gz"
    )

    parser.add_argument(
        "--query",
        required=True,
        help="Input protein FASTA file"
    )

    parser.add_argument(
        "--output",
        default="protein_cluster_mapping.csv",
        help="Output CSV filename "
             "(default: protein_cluster_mapping.csv)"
    )

    args = parser.parse_args()

    # --------------------------------------------------
    # Load gene_data.csv from compressed archive
    # --------------------------------------------------

    print("Loading Panaroo gene data...")

    gd = read_csv_from_tar(
        args.gene_data,
        "gene_data.csv"
    )

    print(f"Loaded {len(gd)} Panaroo gene records")

    # --------------------------------------------------
    # Build protein sequence -> annotation ID lookup
    # --------------------------------------------------

    seq_to_annotation = {}

    for seq, annotation_id in zip(
        gd["prot_sequence"],
        gd["annotation_id"]
    ):

        seq = normalize(seq)

        if seq:
            seq_to_annotation.setdefault(seq, []).append(
                annotation_id
            )

    print(
        f"Indexed {len(seq_to_annotation)} unique protein sequences"
    )

    # --------------------------------------------------
    # Load gene_presence_absence.csv from archive
    # --------------------------------------------------

    print("Loading Panaroo presence/absence data...")

    pa = read_csv_from_tar(
        args.presence,
        "gene_presence_absence.csv"
    )

    print(f"Loaded {len(pa)} gene clusters")

    # --------------------------------------------------
    # Identify genome columns
    # --------------------------------------------------

    metadata_cols = {
        "Gene",
        "Non-unique Gene name",
        "Annotation",
        "No. isolates",
        "No. sequences",
        "Avg sequences per isolate"
    }

    genome_cols = [
        col for col in pa.columns
        if col not in metadata_cols
    ]

    # --------------------------------------------------
    # Build annotation ID -> cluster ID lookup
    # --------------------------------------------------

    annotation_to_cluster = {}

    for _, row in pa.iterrows():

        cluster = row["Gene"]

        for col in genome_cols:

            cell = row[col]

            if pd.isna(cell):
                continue

            # Split multiple gene IDs if present
            annotation_ids = re.split(r"[;,]", str(cell))

            for annotation_id in annotation_ids:

                annotation_id = annotation_id.strip()

                if annotation_id:
                    annotation_to_cluster[annotation_id] = cluster

    print(
        f"Mapped {len(annotation_to_cluster)} annotation IDs "
        f"to clusters"
    )

    # --------------------------------------------------
    # Match query protein sequences
    # --------------------------------------------------

    results = []

    unmatched_records = []

    n_queries = 0
    n_matched = 0

    print("Matching query proteins...")

    for record in SeqIO.parse(args.query, "fasta"):

        n_queries += 1

        seq = normalize(record.seq)

        matched_annotation_ids = seq_to_annotation.get(seq, [])

        if matched_annotation_ids:

            n_matched += 1

            for annotation_id in matched_annotation_ids:

                cluster = annotation_to_cluster.get(annotation_id)

                results.append({
                    "query_id": record.id,
                    "annotation_id": annotation_id,
                    "cluster_id": cluster,
                    "status": (
                        "matched"
                        if cluster is not None
                        else "annotation_matched_cluster_not_found"
                    )
                })

        else:

            results.append({
                "query_id": record.id,
                "annotation_id": None,
                "cluster_id": None,
                "status": "unmatched"
            })

            unmatched_records.append(record)

    # --------------------------------------------------
    # Save mapping output
    # --------------------------------------------------

    out = pd.DataFrame(
        results,
        columns=[
            "query_id",
            "annotation_id",
            "cluster_id",
            "status"
        ]
    )

    out.to_csv(args.output, index=False)

    # --------------------------------------------------
    # Save unmatched proteins to FASTA
    # --------------------------------------------------

    if unmatched_records:

        unmatched_fasta = args.output.replace(
            ".csv",
            "_unmatched.fasta"
        )

        SeqIO.write(
            unmatched_records,
            unmatched_fasta,
            "fasta"
        )

        print(
            f"\nUnmatched proteins saved to: "
            f"{unmatched_fasta}"
        )

    else:

        print("\nAll query proteins matched successfully!")

    # --------------------------------------------------
    # Summary
    # --------------------------------------------------

    print("\nMapping complete")
    print("------------------------------")
    print(f"Queries processed: {n_queries}")
    print(f"Matched queries:   {n_matched}")
    print(f"Unmatched queries: {n_queries - n_matched}")
    print(f"Output file:       {args.output}")

    # --------------------------------------------------
    # Report unmatched protein IDs
    # --------------------------------------------------

    if unmatched_records:

        print("\nUnmatched protein IDs:")

        for record in unmatched_records:
            print(f"  {record.id}")

    else:

        print("\nNo unmatched proteins.")


if __name__ == "__main__":
    main()