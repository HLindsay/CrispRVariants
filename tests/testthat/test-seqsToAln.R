context("Testing seqsToAln correctly truncates alignments")

# Setup data

target <- GenomicRanges::GRanges("1", IRanges::IRanges(10,19), "+")
cigar1 <- "10D20M"
seq1 <- Biostrings::DNAStringSet("AAAAACCCCCTTTTTGGGGG")

test_that("seqsToAln correctly handles overlapping gaps", {
   expect_equal(seqsToAln(cigar1, seq1, target, aln_start = 5), "-----AAAAA")
   expect_equal(seqsToAln(cigar1, seq1, target, aln_start = 0), "AAAAACCCCC")
   expect_equal(seqsToAln(cigar1, seq1, target, 
                  aln_start = 5, reverse_complement = TRUE), "TTTTT-----")
})


