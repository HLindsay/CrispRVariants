context("Text access of attributes from CrisprSet objects")

data("gol_clutch1")

test_that("consensusSeqs returns the correct sequence", {
    expected <- Biostrings::DNAString("GTCTTGGTCTCTCGCAGGATGCTGGAGCCA")
    alleles <- consensusSeqs(gol)
    expect_true(alleles[["-3:3D"]] == expected)
    expect_true(identical(rownames(gol$cigar_freqs), names(alleles)))
})
