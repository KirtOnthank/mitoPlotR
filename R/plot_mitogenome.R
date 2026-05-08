utils::globalVariables(c(
  "type", "seqnames", "Name", "ID", "strand",
  "molecule", "start", "end", "gene", "pos", "bin"
))

#' Plot a circular mitochondrial genome map
#'
#' Creates a circular mitochondrial genome plot from MITOS2 annotation
#' output and a mitochondrial genome FASTA file. The plot includes
#' forward- and reverse-strand feature tracks, optional sequencing
#' coverage depth, coordinate axes, and center labels.
#'
#' This function is designed primarily for metazoan mitochondrial
#' genomes annotated with MITOS2.
#'
#' @param gff_file Path to a MITOS2 GFF annotation file.
#' @param fasta_file Path to the mitochondrial genome FASTA file.
#' @param depth_file Optional path to a tab-delimited sequencing depth
#'   file containing three columns:
#'   sequence name, position, and depth.
#' @param output_file Output image filename. Supported formats
#'   currently include `.png` and `.svg`. The graphics format is
#'   determined automatically from the file extension.
#' @param species_name Species name to display in the center of the plot.
#'   If `NULL`, the FASTA filename will be used.
#' @param width Width of the output image in pixels.
#' @param height Height of the output image in pixels.
#' @param res Resolution of the output image in DPI.
#' @param bin_size Bin size used for averaging sequencing depth.
#' @param feature_cols Named vector of colors used for feature types.
#' @param label_tRNAs Logical indicating whether tRNA genes should be
#'   labeled on the plot.
#'
#' @details
#' The function expects MITOS2 annotation output in GFF format and a
#' corresponding mitochondrial genome FASTA file.
#'
#' Sequencing depth files should contain three tab-delimited columns:
#'
#' \describe{
#'   \item{Column 1}{Sequence name}
#'   \item{Column 2}{Genome position}
#'   \item{Column 3}{Sequencing depth}
#' }
#'
#' Coverage depth is automatically binned and scaled for plotting.
#'
#' @return
#' Invisibly returns the output filename.
#'
#' @examples
#' \dontrun{
#' plot_mitogenome(
#'   gff_file = "result.gff",
#'   fasta_file = "Octopus_rubescens_mitogenome.fasta",
#'   depth_file = "paired.depth.txt",
#'   output_file = "octopus_mitogenome.png",
#'   species_name = "Octopus rubescens"
#' )
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom Biostrings readDNAStringSet width
#' @importFrom circlize circos.clear circos.par circos.initialize
#' @importFrom circlize circos.trackPlotRegion circos.rect
#' @importFrom circlize circos.text circos.axis circos.lines
#' @importFrom dplyr filter mutate select distinct group_by
#' @importFrom dplyr summarise rename if_else
#' @importFrom grDevices png dev.off
#' @importFrom graphics text
#' @importFrom stats quantile
#' @importFrom utils read.table globalVariables
#'
#' @export

plot_mitogenome <- function(
    gff_file,
    fasta_file,
    depth_file = NULL,
    output_file = "circular.png",
    species_name = NULL,
    width = 2000,
    height = 2000,
    res = 500,
    bin_size = 25,
    feature_cols = c(
      gene = "#4DAF4A",
      tRNA = "#377EB8",
      rRNA = "#E41A1C"
    ),
    label_tRNAs = FALSE
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required.")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }

  file_type <- match.arg(file_type, c("png", "svg"))

  gff <- rtracklayer::import(gff_file)

  seqs <- Biostrings::readDNAStringSet(fasta_file)
  genome_length <- Biostrings::width(seqs)[1]
  molecule_name <- names(seqs)[1]

  if (is.null(species_name)) {
    species_name <- tools::file_path_sans_ext(basename(fasta_file))
  }

  genes_plot <- as.data.frame(gff) |>
    dplyr::filter(type %in% c("gene", "tRNA", "rRNA")) |>
    dplyr::mutate(
      molecule = as.character(seqnames),
      gene = dplyr::if_else(is.na(Name), ID, Name),
      strand = as.character(strand)
    ) |>
    dplyr::select(molecule, start, end, strand, gene, type) |>
    dplyr::distinct(molecule, start, end, strand, gene, type, .keep_all = TRUE)

  depth_binned <- NULL

  if (!is.null(depth_file)) {
    depth <- read.table(
      depth_file,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    colnames(depth) <- c("molecule", "pos", "depth")
    depth$molecule <- molecule_name

    depth_binned <- depth |>
      dplyr::mutate(bin = ((pos - 1) %/% bin_size) * bin_size + 1) |>
      dplyr::group_by(molecule, bin) |>
      dplyr::summarise(depth = mean(depth), .groups = "drop") |>
      dplyr::rename(pos = bin)

    depth_max <- stats::quantile(depth_binned$depth, 0.99, na.rm = TRUE)
    depth_binned$depth_scaled <- pmin(depth_binned$depth, depth_max) / depth_max
  }

  file_type <- tolower(tools::file_ext(output_file))

  if (!file_type %in% c("png", "svg")) {
    stop(
      "Unsupported output format. ",
      "Please use a filename ending in '.png' or '.svg'."
    )
  }

  if (file_type == "png") {

    grDevices::png(
      filename = output_file,
      width = width,
      height = height,
      res = res
    )

  } else if (file_type == "svg") {

    grDevices::svg(
      filename = output_file,
      width = width / res,
      height = height / res
    )

  }

  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 90,
    gap.degree = 0,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.01, 0.01)
  )

  circlize::circos.initialize(
    factors = molecule_name,
    xlim = cbind(1, genome_length)
  )

  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.12,
    panel.fun = function(x, y) {
      sector_name <- circlize::CELL_META$sector.index
      feat <- genes_plot |> dplyr::filter(molecule == sector_name)

      circlize::circos.rect(
        xleft = 0, ybottom = 0.52,
        xright = genome_length, ytop = 0.95,
        col = "white", border = "black", lwd = 0.4
      )

      circlize::circos.rect(
        xleft = 0, ybottom = 0.05,
        xright = genome_length, ytop = 0.48,
        col = "white", border = "black", lwd = 0.4
      )

      for (i in seq_len(nrow(feat))) {
        ybottom <- if (feat$strand[i] == "+") 0.52 else 0.05
        ytop    <- if (feat$strand[i] == "+") 0.95 else 0.48

        circlize::circos.rect(
          xleft = feat$start[i],
          ybottom = ybottom,
          xright = feat$end[i],
          ytop = ytop,
          col = unname(feature_cols[as.character(feat$type[i])]),
          border = "black",
          lwd = 0.4
        )
      }
    }
  )

  label_idx <- if (label_tRNAs) {
    seq_len(nrow(genes_plot))
  } else {
    which(genes_plot$type != "tRNA")
  }

  for (i in label_idx) {
    mid <- (genes_plot$start[i] + genes_plot$end[i]) / 2

    circlize::circos.text(
      x = mid,
      y = 1.5,
      labels = genes_plot$gene[i],
      sector.index = genes_plot$molecule[i],
      facing = "bending.outside",
      niceFacing = TRUE,
      adj = c(0.5, 0),
      cex = 0.5
    )
  }

  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    bg.border = NA,
    track.height = 0.05,
    panel.fun = function(x, y) {
      circlize::circos.axis(
        h = "top",
        major.at = seq(0, genome_length, by = 2000),
        labels.cex = 0.5,
        minor.ticks = 0,
        direction = "inside"
      )
    }
  )

  if (!is.null(depth_binned)) {
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1),
      bg.border = NA,
      bg.col = NA,
      track.height = 0.10,
      track.margin = c(0.04, 0.05),
      panel.fun = function(x, y) {
        sector_name <- circlize::CELL_META$sector.index
        d <- depth_binned[depth_binned$molecule == sector_name, ]

        circlize::circos.rect(
          xleft = d$pos,
          ybottom = 0,
          xright = pmin(d$pos + bin_size - 1, genome_length),
          ytop = d$depth_scaled,
          col = "grey30",
          border = "grey30"
        )

        circlize::circos.lines(
          x = c(1, genome_length),
          y = c(0.5, 0.5),
          col = "red",
          lty = 2,
          lwd = 0.7
        )
      }
    )
  }

  graphics::text(0, 0.15, labels = species_name, cex = 1, font = 3)
  graphics::text(0, 0.05, labels = "Mitochondrial genome", cex = 0.6)
  graphics::text(0, -0.05, labels = paste0(genome_length, " bp"), cex = 0.6)

  grDevices::dev.off()

  invisible(output_file)
}
