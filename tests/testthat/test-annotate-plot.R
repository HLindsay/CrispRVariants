context("Test annotation and plotting of targets to genes")

suppressPackageStartupMessages(library(GenomicFeatures))
data(gol_clutch1)

txdb_fname <- system.file("extdata", "Danio_rerio.Zv9.73.gol.sqlite",
                   package = "CrispRVariants")
txdb <- loadDb(txdb_fname)

genes <- CrispRVariants:::.getOverlappingGenes(txdb, gol$target, TRUE)
result <- CrispRVariants:::.makeGeneSegments(genes, txdb, gol$target)
all_exons <- result$all_exs

x <- CrispRVariants:::annotateGenePlot(txdb, gol$target)
y <- CrispRVariants:::annotateGenePlot(txdb, 
                             GenomicRanges::GRanges("18", IRanges(1,10)))

test_that("UTR coordinates for gol are corrrect", {
  expect_equal(all_exons[all_exons$type == "utr","start"], 4644226)
  expect_equal(all_exons[all_exons$type == "utr","end"], 4644285)
})


test_that("Annotate gene plot produces a ggplot if genes match",{
  expect_equal(any(class(x) == "ggplot"), TRUE)
  expect_equal(any(class(x) == "grob"), FALSE)
})


test_that("Annotate gene plot produces an empty grob if no genes match",{
  expect_equal(any(class(y) == "ggplot"), FALSE)
  expect_equal(any(class(y) == "grob"), TRUE)
})
