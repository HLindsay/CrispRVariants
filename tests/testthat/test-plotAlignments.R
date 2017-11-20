context("Test alignment plot")

# Setup data for testing plotAlignments -----
data("gol_clutch1")

ref_only <- plotAlignments(Biostrings::DNAString("AACCTTGG"), alns = NULL,
                  ins.sites = data.frame())

crispr_plot <- plotAlignments(gol)
# -----

test_that("plotAlignments returns a ggplot object", {
    # plotAlignments should produce a ggplot that has layers
    expect_is(ref_only, "ggplot")
    expect_true(nrow(ref_only$data) > 0)
    
    # The same should be true when plotAlignments is called with a CrisprSet
    expect_is(crispr_plot, "ggplot")
    expect_true(nrow(crispr_plot$data) > 0)
    
    # Original sequence should be present
    expect_identical(as.character(ref_only$data$value), 
                     strsplit("AACCTTGG", "")[[1]])
})
