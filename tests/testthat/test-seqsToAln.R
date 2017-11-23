context("Testing seqsToAln correctly truncates alignments")

# Setup data -----

target <- GenomicRanges::GRanges("1", IRanges::IRanges(10,19), "+")
cigar1 <- "10D20M"
seq1 <- Biostrings::DNAStringSet("AAAAACCCCCTTTTTGGGGG")
# -----

test_that("seqsToAln correctly handles overlapping gaps", {
    expect_equal(seqsToAln(cigar1, seq1, target, aln_start = 5), "-----AAAAA")
    expect_equal(seqsToAln(cigar1, seq1, target, aln_start = 0), "AAAAACCCCC")
    expect_equal(seqsToAln(cigar1, seq1, target, 
                  aln_start = 5, reverse_complement = TRUE), "TTTTT-----")
})


# Data for examples with partial alignments -----
data("gol_clutch1")
gol_alns <- unlist(alns(gol))
aln_seqs <- CrispRVariants:::seqsToAln(cigar(gol_alns),
                                       dnaseq = mcols(gol_alns)$seq,
                                       target = gol$target,
                                       aln_start = start(gol_alns))

keep_wrt_target <- GRanges(18, IRanges(start = c(4647372, 4647396),
                                       end = c(4647392,4647404)))

target <- gol$target

keep <- IRanges(start = c(1,25), end = c(21,33))
# Delete characters at the start
keep_l <- IRanges(start = c(5,25), end = c(21,33))
# Delete characters at the end
keep_r <- IRanges(start = c(1,25), end = c(21,30))
# Delete characters at both ends
keep_b <- IRanges(start = c(5,25), end = c(21,30))

#-----


test_that("selectAlnRegions converts GRanges if possible", {
    expect_identical(.checkRelativeLocs(target, keep_wrt_target), keep)
    
    # Function should fail with a warning if genomic coordinates are given
    # as IRanges
    expect_warning(result <- selectAlnRegions(5, seq1, seq1, target, 
                                              ranges(keep_wrt_target)))
    expect_equal(result, NULL)
})

test_that("Insertions on the left boundary of a gap are retained", {
    original_ins <- gol$insertion_sites$start
    result <- .adjustRelativeInsLocs(gol$target, keep, original_ins, c(0,5))
    expect_equal(length(result), length(original_ins))
    
    # Check insertions at first gap base pos 22 are kept 
    expect_true(! any(is.na(result[original_ins == 22])))
    
    # Check insertions at last gap base 24 are ignored
    expect_true(all(is.na(result[original_ins == 24])))
})

test_that("selectAlnRegions correctly joins multiple segments", {
    select_args <- list(alns = aln_seqs[1:10], reference = gol$ref,
                        target = gol$target, join = " del %s bp ")
  
    run_selectRegions <- function(keep){
       do.call(selectAlnRegions,
               modifyList(select_args, list(keep = keep)))
    }
    
    result <- run_selectRegions(keep)
    result_l <- run_selectRegions(keep_l)
    result_r <- run_selectRegions(keep_r)
    result_b <- run_selectRegions(keep_b)
    # Same thing allowing gaps at borders
    select_args["border.gaps"] <- TRUE
    result_bd <- run_selectRegions(keep)
    result_lbd <- run_selectRegions(keep_l)
    result_rbd <- run_selectRegions(keep_r)
    result_bbd <- run_selectRegions(keep_b)
    
    expect_equal(result$ref, "GTCTTGGTCTCTCGCAGGATG del 3 bp CTGGAGCCA")
    expect_equal(result_l$ref, "TGGTCTCTCGCAGGATG del 3 bp CTGGAGCCA" )
    expect_equal(result_r$ref, "GTCTTGGTCTCTCGCAGGATG del 3 bp CTGGAG")
    expect_equal(result_b$ref, "TGGTCTCTCGCAGGATG del 3 bp CTGGAG")
    expect_equal(result_bd$ref, "GTCTTGGTCTCTCGCAGGATG del 3 bp CTGGAGCCA")
    expect_equal(result_lbd$ref,
                 " del 4 bp TGGTCTCTCGCAGGATG del 3 bp CTGGAGCCA")
    expect_equal(result_rbd$ref,
                 "GTCTTGGTCTCTCGCAGGATG del 3 bp CTGGAG del 3 bp ")
    expect_equal(result_bbd$ref,
                 " del 4 bp TGGTCTCTCGCAGGATG del 3 bp CTGGAG del 3 bp ")
})


test_that("selectAlnRegions works starting from partial alignments", {
  
  
})
