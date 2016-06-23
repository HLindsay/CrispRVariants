# Setup data
context("Initialization of CrisprSet objects")



test_that("readsToTargets correctly separates reads by PCR primer",{
  wdths <- c(10,10,5,3)
  gr <- GenomicRanges::GRanges("chr1", IRanges(start = c(5,7,4,13), width = wdths),
                               cigar = sprintf("%sM", wdths), strand = "+", 
                               flag=c(0,0,0,2048))
  names(gr) <- c("A", "B", "C","C")
  seqs <- Biostrings::DNAStringSet(subseq(rep("ACTGACTGAC", 
                                              length(gr)),1, width(gr)))
  gr$seq <- seqs
  galnsl <- GenomicAlignments::GAlignmentsList(list(as(gr, "GAlignments")))
  targets <- GenomicRanges::GRanges("chr1", IRanges(c(10, 12), width = 2))
  primer.ranges <- GenomicRanges::GRanges("chr1", IRanges(c(4, 7), width = c(11,10)))
  references <- Biostrings::DNAStringSet(c("AA","CC"))

  # test that reads can't be separated without primer ranges
  csets <- suppressWarnings(readsToTargets(galnsl, targets, 
             references = references, target.loc = 1, verbose = FALSE))
  # There should be no reads as all are ambiguous
  # and chimeras cannot be distinguished with default tolerance of 5
  expect_equal(length(csets), 0)
  
  csets <- readsToTargets(galnsl, targets, references = references, target.loc = 1,
                          primer.ranges = primer.ranges, chimera.to.target = 0,
                          verbose = FALSE)
  
  expect_equal(names(csets[[1]]$crispr_runs[[1]]$alns), "A")
  expect_equal(names(csets[[2]]$crispr_runs[[1]]$alns), "B")
  # Chimeras can be resolved with zero tolerance
  expect_equal(names(csets[[2]]$crispr_runs[[1]]$chimeras),c("C","C"))
  expect_equal(length(csets[[1]]$crispr_runs[[1]]$chimeras), 0)
})


test_that("Arguments for calling SNVs are passed on",{
  # There is a nucleotide change 19 bases upstream of the cut site
  # in one of these sequences
  bam <- system.file("extdata", "bam/ab1_ptena_wildtype_looking_embryo_1_s.bam",
                      package="CrispRVariants")
  reference <- Biostrings::DNAString("GCCATGGGCTTTCCAGCCGAACGATTGGAAGGT")
  gdl <- GenomicRanges::GRanges("chr17", IRanges(23648469, 23648501),
           strand = "-")
  # With default settings, the snv should not be called
  cset <- readsToTarget(bam, target = gdl, reference = reference,
           target.loc = 22)
  expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 0)

  # Re-initialise increasing detection window.
  # One sequence should now have a SNV
  cset <- readsToTarget(bam, target = gdl, reference = reference,
           target.loc = 22, upstream.snv = 20)
  expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 1)

})


# To do:
# test that reads are not separated when exact match required (allow.partial = FALSE)
# test separateChimeras that chimera away from cut is excluded
# same test via readsToTarget
# test separateChimeras that guide within chimera is not chimeric
# test function of chimera.to.target
# test chimera options "count","exclude","ignore", "merge"

