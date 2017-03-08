context("Testing recreation of reference from alignments")

sq_ln <- c(500)
names(sq_ln) <- "1"

alns <- GenomicAlignments::GAlignments(seqnames=rep("1",6),
                                       pos=as.integer(seq(1,30,5)),
    cigar=c("10M1I11M", "9M1I10M", "10M1D5M", "3M1I7M1D9M", "2M3I2M3I2M","5M1I3M"),
    strand = S4Vectors::Rle(factor(rep("+", 6), levels = c("+","-","*"))),
    seqlengths = sq_ln,
    MD = c("22", "5A4A8", "10^A5", "A1A1A2A1A^A3A2AA1", "AA1A2", "7A"),
    seq = Biostrings::DNAStringSet(c("AAAAAAAAAATAAAAAAAAAAA",
                                     "AAAAATAAACAGAAAAAAAA",
                                     "AAAAAAAAAAAAAAA",
                                     "CAGTAGAACACAAATAAGGA",
                                     "TTCCCAGCCCAA",
                                     "AAAAAGAAT"))
)


test_that("Recreation of reference when location not provided", {
    lf <- Biostrings::letterFrequency(refFromAlns(alns), letters = c("A"), as.prob = TRUE)
    expect_equal(unname(unique(lf)[,"A"]), 1)
})


test_that("Recreation of reference when location provided", {
    test_loc <- GenomicRanges::GRanges("1", IRanges::IRanges(10,15))
    test_outside_loc <- GenomicRanges::GRanges("1", IRanges::IRanges(30,40))
    expect_equal(refFromAlns(alns, location = test_loc), Biostrings::DNAString("AAAAAA"))
    expect_error(refFromAlns(alns, location = test_outside_loc))
})
