
test_that("plot_mitogenome creates a PNG file", {
  gff <- system.file("extdata", "example.gff", package = "mitoPlotR")
  fasta <- system.file("extdata", "example.fasta", package = "mitoPlotR")
  depth <- system.file("extdata", "example.depth.txt", package = "mitoPlotR")

  out <- tempfile(fileext = ".png")

  plot_mitogenome(
    gff_file = gff,
    fasta_file = fasta,
    depth_file = depth,
    output_file = out,
    species_name = "Octopus rubescens"
  )

  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})


test_that("plot_mitogenome creates an SVG file", {
  gff <- system.file("extdata", "example.gff", package = "mitoPlotR")
  fasta <- system.file("extdata", "example.fasta", package = "mitoPlotR")
  depth <- system.file("extdata", "example.depth.txt", package = "mitoPlotR")

  out <- tempfile(fileext = ".svg")

  plot_mitogenome(
    gff_file = gff,
    fasta_file = fasta,
    depth_file = depth,
    output_file = out,
    species_name = "Octopus rubescens"
  )

  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
})


test_that("unsupported output formats produce an error", {
  gff <- system.file("extdata", "example.gff", package = "mitoPlotR")
  fasta <- system.file("extdata", "example.fasta", package = "mitoPlotR")

  out <- tempfile(fileext = ".pdf")

  expect_error(
    plot_mitogenome(
      gff_file = gff,
      fasta_file = fasta,
      output_file = out
    ),
    "Unsupported output format"
  )
})

