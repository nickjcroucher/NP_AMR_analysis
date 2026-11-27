import pandas as pd
import tarfile
import re

# Input files
panaroo_file = "hmrM_proteins.csv"
tar_file = "Hflu_gene_data.csv"

# 1. Read the Panaroo file
panaroo = pd.read_csv(panaroo_file, header=None, dtype=str)

# Extract gene IDs, e.g. MHINDH_00539
pattern = re.compile(r"^[A-Z]{6}_\d{5}$")

gene_ids = []
for value in panaroo.to_numpy().flatten():
    if pd.notna(value):
        value = value.strip()
        if pattern.match(value):
            gene_ids.append(value)

# Remove duplicates while preserving order
gene_ids = list(dict.fromkeys(gene_ids))

print(f"Gene IDs requested: {len(gene_ids)}")

# 2. Read Hflu_gene_data.csv directly from the tar.gz archive

data = pd.read_csv(tar_file, low_memory=False)

# 3. Match IDs using annotation_id
matched = data[data["annotation_id"].isin(gene_ids)].copy()

# Preserve the order of IDs in the Panaroo file
matched["order"] = matched["annotation_id"].map(
    {gene: i for i, gene in enumerate(gene_ids)}
)
matched = matched.sort_values("order")

# 4. Write protein sequences to FASTA
with open("hmrM_proteins.fasta", "w") as f:
    for _, row in matched.iterrows():
        f.write(f">{row['annotation_id']} {row['gene_name']} | {row['description']}\n")
        f.write(f"{row['prot_sequence']}\n")

# 5. (Optional) Write nucleotide sequences to FASTA
with open("hmrM_proteins.fasta", "w") as f:
    for _, row in matched.iterrows():
        f.write(f">{row['annotation_id']} {row['gene_name']}\n")
        f.write(f"{row['dna_sequence']}\n")

# 6. Report missing IDs
found_ids = set(matched["annotation_id"])
missing_ids = set(gene_ids) - found_ids

print(f"Sequences extracted: {len(matched)}")
print(f"IDs not found: {len(missing_ids)}")

if missing_ids:
    with open("missing_hmrM_ids.txt", "w") as f:
        f.write("\n".join(sorted(missing_ids)))
    print("Missing IDs saved to missing_hmrM_ids.txt")

print("Finished!")

# Diagnose why more rows were extracted than requested

print("\n--- DIAGNOSTICS ---")

print(f"Unique gene IDs requested: {len(gene_ids)}")
print(f"Rows returned: {len(matched)}")
print(f"Unique annotation IDs returned: {matched['annotation_id'].nunique()}")

# Check for duplicated annotation IDs in the extracted rows
duplicates = matched[
    matched["annotation_id"].duplicated(keep=False)
].sort_values("annotation_id")

print(f"\nDuplicated annotation IDs: {duplicates['annotation_id'].nunique()}")

if not duplicates.empty:
    print("\nDuplicated entries:")
    print(duplicates[[
        "annotation_id",
        "gff_file",
        "scaffold_name",
        "clustering_id",
        "gene_name",
        "description"
    ]].to_string(index=False))

# Check how many times each requested ID occurs in the source data
counts = data[data["annotation_id"].isin(gene_ids)]["annotation_id"].value_counts()

multiple_entries = counts[counts > 1]

print(f"\nIDs occurring multiple times in source data: {len(multiple_entries)}")

if not multiple_entries.empty:
    print(multiple_entries.to_string())

# Check whether there are any unexpected IDs
unexpected = set(matched["annotation_id"]) - set(gene_ids)

print(f"\nUnexpected IDs returned: {len(unexpected)}")
if unexpected:
    print(unexpected)

# Check whether sequences are identical for duplicated IDs
if not duplicates.empty:
    print("\nSequence differences among duplicate IDs:")

    for gene_id in duplicates["annotation_id"].unique():
        subset = duplicates[duplicates["annotation_id"] == gene_id]

        print(
            gene_id,
            "rows =", len(subset),
            "unique protein sequences =", subset["prot_sequence"].nunique(),
            "unique DNA sequences =", subset["dna_sequence"].nunique()
        )