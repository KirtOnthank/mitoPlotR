
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mitoPlotR

<!-- badges: start -->

<!-- badges: end -->

mitoPlotR is an R package for visualization and exploratory analysis of
mitochondrial genomes, with a focus on metazoan mitogenomes and MITOS2
annotation workflows

## Installation

You can install mitoPlotR from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("KirtOnthank/mitoPlotR")
```

## Example

``` r
library(mitoPlotR)

plot_mitogenome(
  gff_file = system.file("extdata", "example.gff", package = "mitoPlotR"),
  fasta_file = system.file("extdata", "example.fasta", package = "mitoPlotR"),
  depth_file = system.file("extdata", "example.depth.txt", package = "mitoPlotR"),
  output_file = "example_plot.png",
  species_name = "Octopus rubescens"
)
```

## Example output

![](man/figures/example_plot.png)
