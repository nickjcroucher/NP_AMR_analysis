#' Plot a gggenomes locus comparison with flip support and separate red/blue
#' identity-gradient legends (red = same, blue = opposite direction).
#'
#' @param seqs         sequence tibble (from read_seqs())
#' @param feats        feature tibble (from read_feats())
#' @param links        link tibble (from read_paf()), must contain
#'                      seq_id, seq_id2, strand, de
#' @param flip_bins    character vector of bin_id(s) to flip (default: none).
#'                      With sync = TRUE, flipping the first plotted bin before
#'                      sync makes it the orientation reference for other bins.
#' @param sync         use gggenomes::sync() to orient bins automatically (default TRUE).
#' @param sync_link_track link track for sync; 1 is the PAF alignment track,
#'                      2 is the cluster track when clusters are supplied.
#' @param clusters     optional gggenomes cluster-membership table with feat_id
#'                      and cluster_id (one row per gene); not raw BLAST hits.
#' @param min_identity lower bound for both fill colourbars. Colourbars start
#'                      at the lowest observed match at or above this bound;
#'                      lower identities are clamped (default 0.90).
#' @param legend_y     vertical position of the legends within the panel (0 =
#'                      bottom, 1 = top; default 0.03).
#' @param scalebar_y   vertical position of the right-hand scale bar. NULL
#'                      (default) tracks legend_y, offset by 0.05.
#' @param legend_space space below the lowest sequence in bin-row units. NULL
#'                      (default) reserves space based on the number of bins;
#'                      increase manually for very short export heights.
#' @param sequence_names named list/vector: bin_id -> displayed bin label.
#'                        Plain text is displayed literally; valid plotmath
#'                        strings such as 'italic("Genus species")~"isolate"'
#'                        are formatted. Does not modify genomic IDs.
#' @param gene_annotations named list/vector: feat_id -> label shown above gene.
#'                         Plain text is displayed literally; plotmath labels
#'                         such as "italic(gene)" are formatted.
#' @param gene_colours named list/vector: feat_id -> fill colour. When clusters
#'                     are supplied, only these colours propagate within their
#'                     cluster; otherwise uncoloured genes stay white.
#' @param gene_label_nudge_y vertical offset of gene labels in bin-row units
#'                           (default 0.12); increase to raise them further.
plot_locus_comparison <- function(seqs, feats, links,
                                  flip_bins = character(0),
                                  min_identity = 0.90,
                                  legend_y = 0.03,
                                  scalebar_y = NULL,
                                  legend_space = NULL,
                                  sequence_names = NULL,
                                  gene_annotations = NULL,
                                  gene_colours = NULL,
                                  gene_label_nudge_y = 0.12,
                                  sync = TRUE,
                                  clusters = NULL,
                                  sync_link_track = 1) {
  
  if (is.null(scalebar_y)) scalebar_y <- legend_y + 0.05
  stopifnot(
    length(legend_y) == 1, is.finite(legend_y), legend_y >= 0, legend_y <= 1,
    length(scalebar_y) == 1, is.finite(scalebar_y), scalebar_y >= 0, scalebar_y <= 1,
    length(gene_label_nudge_y) == 1, is.finite(gene_label_nudge_y),
    gene_label_nudge_y >= 0,
    is.logical(sync), length(sync) == 1L, !is.na(sync),
    length(sync_link_track) == 1L, sync_link_track %in% c(1, 2)
  )
  if (!is.null(clusters)) {
    stopifnot(is.data.frame(clusters),
              all(c("feat_id", "cluster_id") %in% names(clusters)))
    clusters <- dplyr::distinct(clusters, feat_id, cluster_id)
    if (anyDuplicated(clusters$feat_id)) {
      stop("Each feat_id must belong to at most one cluster.")
    }
  }
  if (sync && sync_link_track == 2 && is.null(clusters)) {
    stop("sync_link_track = 2 requires clusters.")
  }
  if (!is.null(legend_space)) {
    stopifnot(length(legend_space) == 1, is.finite(legend_space), legend_space > 0)
  }
  if (is.null(legend_space)) {
    legend_space <- max(1.2, 0.3 * dplyr::n_distinct(seqs$bin_id))
  }
  
  # Mappings change display only; genomic IDs remain intact for links/flip.
  as_named_map <- function(x) {
    if (is.null(x)) return(NULL)
    x <- unlist(x, use.names = TRUE)
    stopifnot(is.character(x), !is.null(names(x)),
              all(nzchar(names(x))), !anyDuplicated(names(x)))
    x
  }
  bin_labels <- as_named_map(sequence_names)
  if (!is.null(bin_labels)) {
    bin_labels <- vapply(bin_labels, function(label) {
      # Accept italic(genus species) as shorthand for italic("genus species").
      label <- gsub(
        'italic\\(([[:alnum:]_]+(?:[[:space:]]+[[:alnum:]_]+)+)\\)',
        'italic("\\1")', label, perl = TRUE
      )
      expr <- tryCatch(parse(text = label), error = function(e) NULL)
      if (length(expr) == 1L && is.call(expr[[1L]])) label
      else encodeString(label, quote = '"')
    }, character(1))
  }
  annotation_map <- as_named_map(gene_annotations)
  colour_map <- as_named_map(gene_colours)
  
  if (!is.null(annotation_map)) {
    feats$gene_annotation <- unname(annotation_map[as.character(feats$feat_id)])
    feats$gene_annotation <- vapply(feats$gene_annotation, function(label) {
      if (is.na(label)) return(NA_character_)
      expr <- tryCatch(parse(text = label), error = function(e) NULL)
      if (length(expr) == 1L && is.call(expr[[1L]])) label
      else encodeString(label, quote = '"')
    }, character(1))
  }
  if (!is.null(colour_map)) {
    feats$gene_fill <- unname(colour_map[as.character(feats$feat_id)])
    if (!is.null(clusters)) {
      seeds <- clusters %>%
        mutate(gene_fill = unname(colour_map[as.character(feat_id)])) %>%
        filter(!is.na(gene_fill)) %>%
        distinct(cluster_id, gene_fill)
      if (anyDuplicated(seeds$cluster_id)) {
        stop("Conflicting gene_colours assigned to the same cluster.")
      }
      feats$gene_fill <- dplyr::coalesce(
        feats$gene_fill,
        unname(setNames(seeds$gene_fill, as.character(seeds$cluster_id))[
          as.character(clusters$cluster_id[match(feats$feat_id, clusters$feat_id)])
        ])
      )
    }
    feats$gene_fill[is.na(feats$gene_fill)] <- "white"
  }
  
  # Keep PAF strand unchanged for ribbon geometry. Colour is computed after
  # automatic and manual flips, from the actual final sequence orientations.
  links <- links %>% mutate(identity = 1 - de)
  
  observed_identity <- links$identity[is.finite(links$identity)]
  colour_min <- if (length(observed_identity)) {
    min(min_identity, min(observed_identity))
  } else {
    min_identity
  }
  if (colour_min >= 1) colour_min <- min_identity
  
  make_plot <- function(link_data) {
    q <- gggenomes::gggenomes(
      seqs  = seqs,
      genes = feats,
      links = link_data,
      adjacent_only = TRUE
    )
    if (!is.null(clusters)) q <- gggenomes::add_clusters(q, clusters)
    
    # Anchor sync to the requested orientation of the first plotted bin.
    # gggenomes::sync() orients later bins relative to the first bin.
    plotted_seqs <- gggenomes::pull_seqs(q)
    first_bin <- as.character(plotted_seqs$bin_id[which.max(plotted_seqs$y)])
    flip_first <- first_bin %in% flip_bins
    if (sync && flip_first) {
      q <- gggenomes::flip(q, tidyselect::all_of(first_bin))
    }
    if (sync) q <- gggenomes::sync(q, link_track = sync_link_track)
    
    # Without sync, flip every requested bin directly. With sync, apply
    # additional manual flips afterwards; do not undo the reference bin.
    remaining_flips <- if (sync && flip_first) setdiff(flip_bins, first_bin) else flip_bins
    for (b in remaining_flips) {
      q <- gggenomes::flip(q, tidyselect::all_of(b))
    }
    q
  }
  
  # Compare actual sequence orientation before/after sync and explicit flips.
  # Build twice so link display_strand can be set on the input links table;
  # gggenomes still receives the unmodified PAF strand for ribbon geometry.
  original <- gggenomes::gggenomes(seqs = seqs, genes = feats, links = links,
                        adjacent_only = TRUE)
  p <- make_plot(links)
  before <- gggenomes::pull_seqs(original) %>%
    select(seq_id, strand) %>%
    mutate(before_strand = as.character(strand))
  after <- gggenomes::pull_seqs(p) %>%
    select(seq_id, strand) %>%
    mutate(after_strand = as.character(strand))
  orientation <- before %>%
    select(seq_id, before_strand) %>%
    inner_join(after %>% select(seq_id, after_strand), by = "seq_id") %>%
    mutate(flipped = before_strand != after_strand)
  flipped_ids <- orientation$seq_id[!is.na(orientation$flipped) & orientation$flipped]
  links$display_strand <- if_else(
    xor(links$seq_id %in% flipped_ids, links$seq_id2 %in% flipped_ids),
    if_else(links$strand == "+", "-", "+"),
    as.character(links$strand)
  )
  p <- make_plot(links)
  
  # Keep the original strand for geometry; colour by post-flip orientation.
  p <- p + gggenomes::geom_seq()
  if (is.null(colour_map) && is.null(clusters)) {
    p <- p + gggenomes::geom_gene(fill = "white")
  } else if (is.null(colour_map)) {
    p <- p + gggenomes::geom_gene(aes(fill = cluster_id), show.legend = FALSE) +
      ggplot2::scale_fill_discrete(na.value = "white", guide = "none") +
      ggnewscale::new_scale_fill()
  } else {
    p <- p + gggenomes::geom_gene(aes(fill = gene_fill), show.legend = FALSE) +
      ggplot2::scale_fill_identity(guide = "none") +
      ggnewscale::new_scale_fill()
  }
  if (!is.null(annotation_map)) {
    p <- p + ggplot2::geom_text(
      data = gggenomes::genes(!is.na(gene_annotation)),
      aes(x = (x + xend) / 2, y = y + gene_label_nudge_y,
          label = gene_annotation),
      parse = TRUE, angle = 45, hjust = 0, vjust = 0,
      check_overlap = TRUE, show.legend = FALSE, inherit.aes = FALSE
    )
  }
  p <- p +
    ggnewscale::geom_link(
      data = gggenomes::links(1, display_strand == "+"),
      aes(fill = identity),
      colour = "white",
      linewidth = 0.3
    ) +
    scale_fill_gradient(
      name   = "Identity \u2013 same direction",
      low    = "#fee8e8",
      high   = "#b30000",
      limits = c(colour_min, 1),
      breaks = c(colour_min, (colour_min + 1) / 2, 1),
      oob    = scales::squish,
      labels = scales::label_percent(accuracy = 0.1),
      guide  = guide_colourbar(
        order = 1, direction = "horizontal", title.position = "top",
        barwidth = grid::unit(3, "cm"), barheight = grid::unit(0.35, "cm")
      )
    ) +
    ggnewscale::new_scale_fill() +
    gggenomes::geom_link(
      data = gggenomes::links(1, display_strand == "-"),
      aes(fill = identity),
      colour = "white",
      linewidth = 0.3
    ) +
    scale_fill_gradient(
      name   = "Identity \u2013 opposite direction",
      low    = "#e3f2fd",
      high   = "#084594",
      limits = c(colour_min, 1),
      breaks = c(colour_min, (colour_min + 1) / 2, 1),
      oob    = scales::squish,
      labels = scales::label_percent(accuracy = 0.1),
      guide  = guide_colourbar(
        order = 2, direction = "horizontal", title.position = "top",
        barwidth = grid::unit(3, "cm"), barheight = grid::unit(0.35, "cm")
      )
    )
  
  # Reserve a dedicated band under the lowest bin. Expanding only the y range
  # keeps every sequence and link at its original genomic plot coordinates.
  if (is.null(bin_labels)) {
    p <- p + gggenomes::geom_bin_label()
  } else {
    p <- p + gggenomes::geom_bin_label(aes(label = ifelse(
      bin_id %in% names(bin_labels),
      unname(bin_labels[as.character(bin_id)]),
      encodeString(as.character(bin_id), quote = '"')
    )), parse = TRUE)
  }
  p +
    ggplot2::expand_limits(y = min(gggenomes::pull_seqs(p)$y) - legend_space) +
    gggenomes::geom_scalebar(
      x = 0.90, y = scalebar_y, length = 0.10,
      label_position = "above", guide_x = "none"
    ) +
    labs(x = NULL) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.02, legend_y),
      legend.justification = c(0, 0),
      legend.box = "horizontal",
      legend.spacing.x = grid::unit(0.4, "cm")
    )
}