context("CrisprSet methods")


test_that("Variant alleles can be renamed",{
    data("gol_clutch1")
    vc <- variantCounts(gol)
    expect_false("renamed" %in% rownames(vc))
    gol_alns <- alns(gol)
    gol_labels <- mcols(unlist(gol_alns))$allele
    gol_labels[gol_labels == "-3:3D"] <- "Renamed"
    gol$setCigarLabels(labels = relist(gol_labels, gol_alns))
})
