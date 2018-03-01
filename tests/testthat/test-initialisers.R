# Setup data
context("Initialization of CrisprSet objects")


test_that("readsToTargets correctly separates reads by PCR primer",{    
    wdths <- c(10,10,5,3)
    gr <- GenomicRanges::GRanges("chr1", 
               IRanges::IRanges(start = c(5,7,4,13), width = wdths),
               cigar = sprintf("%sM", wdths), strand = "+", 
               flag=c(0,0,0,2048))
    names(gr) <- c("A", "B", "C","C")
    seqs <- Biostrings::DNAStringSet(subseq(rep("ACTGACTGAC",
                                            length(gr)),1, width(gr)))
    gr$seq <- seqs
    galnsl <- GenomicAlignments::GAlignmentsList(list(as(gr, "GAlignments")))
    targets <- GenomicRanges::GRanges("chr1",
                        IRanges::IRanges(c(10, 12), width = 2), strand = "+")
    primer.ranges <- GenomicRanges::GRanges("chr1", 
                        IRanges::IRanges(c(4, 7), width = c(11,10)))
    references <- Biostrings::DNAStringSet(c("AA","CC"))

    # test that reads can't be separated without primer ranges
    csets <- suppressWarnings(readsToTargets(galnsl, targets,
               references = references, target.loc = 1, verbose = FALSE))
    # There should be no reads as all are ambiguous
    # and chimeras cannot be distinguished with default tolerance of 5
    expect_equal(length(csets), 0)

    csets <- readsToTargets(galnsl, targets,
                            references = references,
                            target.loc = 1, primer.ranges = primer.ranges,
                            chimera.to.target = 0, verbose = FALSE)
  
    expect_equal(names(csets[[1]]$crispr_runs[[1]]$alns), "A")
    expect_equal(names(csets[[2]]$crispr_runs[[1]]$alns), "B")
    # Chimeras can be resolved with zero tolerance
    expect_equal(names(csets[[2]]$crispr_runs[[1]]$chimeras),c("C","C"))
    expect_equal(length(csets[[1]]$crispr_runs[[1]]$chimeras), 0)
})

# Construct data for checking SNVs -----
bam <- system.file("extdata", "bam/ab1_ptena_wildtype_looking_embryo_1_s.bam",
                      package="CrispRVariants")
reference <- Biostrings::DNAString("GCCATGGGCTTTCCAGCCGAACGATTGGAAGGT")
gdl <- GenomicRanges::GRanges("chr17", IRanges(23648469, 23648501),
           strand = "-")

# Constructed alignments for checking snvs
snv <- rep(GRanges("chr17", IRanges(start = 23648469, width = 33),
                   cigar = "33M"),5)
strand(snv) <- c("+","-","+","-","+")
names(snv) <- c("u8","u7","d5","d6","rc")
mcols(snv)$seq <-  Biostrings::DNAStringSet(
                     c("ACCTTCCAATCGTTCGGCTCGAAAGCCCATGGC",
                       "ACCTTCCAATCGTTCGGCCGGAAAGCCCATGGC",
                       "ACCTTGCAATCGTTCGGCTGGAAAGCCCATGGC",
                       "ACCGTCCAATCGTTCGGCTGGAAAGCCCATGGC",
                       "GCCATGGGCTTTCCAGCCGAACGATTGGAAGGT"))
# ----

test_that("Excluding reads by name works correctly", {
    # This bam file contains four non-chimeric and one chimeric read.
    # Here the non-chimeric reads are excluded by name
    
    # To do: exclude.names gives a warning as is.null used with vector
    cset <- suppressWarnings(readsToTarget(bam, gdl, reference = reference,
                exclude.names = c("AB2017","AB2018","AB2021","AB2024")))

    expect_equal(length(cset$crispr_runs[[1]]$alns), 0)
    expect_equal(length(cset$crispr_runs[[1]]$chimeras), 2)
})


test_that("Ambiguous nucleotides are not considered SNVs",{
    # There is an ambiguous nucleotide 19 bases upstream of the cut site
    # in one of these sequences
    cset <- readsToTarget(bam, target = gdl, reference = reference,
                          target.loc = 22, upstream.snv = 18)
    expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 0)
})


test_that("Arguments for calling SNVs are passed on",{
    # There is a nucleotide mismatch at position 19 upstream of target.loc
    
    # Counting 7 bases downstream, no SNVs are present
    cset <- readsToTarget(bam, target = gdl, reference = reference,
                          target.loc = 22, downstream.snv = 7)
    expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 0)

    # Re-initialise increasing detection window.
    # SNV should now be detected
    cset <- readsToTarget(bam, target = gdl, reference = reference,
                          target.loc = 22, upstream.snv = 19)
    expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 1)
    
    # Check that large detection ranges don't cause and error
    cset <- readsToTarget(bam, target = gdl, reference = reference,
                          target.loc = 22, upstream.snv = 30,
                          downstream.snv = 40)
    expect_equal(length(grep("SNV", cset$crispr_runs[[1]]$cigar_labels)), 1)
    
    # More checking with constructed alignments
    snv <- GenomicAlignments::GAlignmentsList(sample1 = as(snv, "GAlignments"))
    cset <- readsToTarget(snv, target = gdl,
                          reference = reference,
                          target.loc = 22)
    
    # Test that variant at position 8 downstream is not counted
    expect_true(all(c("SNV:-8G", "SNV:6C") %in%
                      rownames(variantCounts(cset))))

    expect_equal(variantCounts(cset)["no variant",], 2)
    
    # With changed parameters, neither of these SNVs are counted
    cset2 <- readsToTarget(snv, target = gdl, reference = reference,
                           target.loc = 22, upstream.snv = 7,
                           downstream.snv = 5)
    expect_equal(variantCounts(cset2)["no variant",], 4)
    
})

test_that("Mismatched reference and target detected",{
    snv <- GenomicAlignments::GAlignmentsList(sample1 = as(snv, "GAlignments"))
    cset <- readsToTarget(snv, target = gdl,
                          reference = Biostrings::reverseComplement(reference),
                          target.loc = 22)
    # Expect: there are no indel sequences
    expect_equal(length(variantCounts(cset, include.nonindel = FALSE)), 0)
    
    # Expect: only one non-snv read
    # (test via mutationEfficiency to also test snv recognition)
    expect_equal(mutationEfficiency(cset, snv = "exclude")[["ReadCount"]], 1)
  
})

test_that("readsToTargets (signature character) returns a list of CrisprSets",{
    # Here testing with a single bam, but should still run
    bam <- system.file("extdata",
                       "bam/ab1_ptena_wildtype_looking_embryo_2_s.bam",
                        package="CrispRVariants")

    refs <- Biostrings::DNAStringSet(c("GCCATGGGCTTTCCAGCCGAACGATTGGAAGGT",
                                        "AAAAACCCCCTTTTTGGGGG"))
    # gol plus a dummy guide
    gdl <- GenomicRanges::GRanges(c("chr17", "chr1"),
                IRanges(start = c(23648469, 1), end = c(23648501, 20)),
           strand = c("-", "+"))
    # This will warn that no reads map to the dummy guide
    csets <- suppressWarnings(readsToTargets(bam, targets = gdl,
                references = refs, target.loc = 17, verbose = FALSE))

    expect_true(class(csets[[1]]) == "CrisprSet")
})

test_that("Strand is chosen correctly", {
    expect_equal(rcAlns("+","target"), FALSE)
    expect_equal(rcAlns("-","target"), TRUE)
    expect_equal(rcAlns("+","opposite"), TRUE)
    expect_equal(rcAlns("-","opposite"), FALSE)
    expect_equal(rcAlns("+","positive"), FALSE)
    expect_equal(rcAlns("-","positive"), FALSE)
})




# To do:
# test that reads are not separated when exact match required (allow.partial = FALSE)
# test separateChimeras that chimera away from cut is excluded
# same test via readsToTarget
# test separateChimeras that guide within chimera is not chimeric
# test function of chimera.to.target
# test chimera options "count","exclude","ignore", "merge"
# test snvs with indels

