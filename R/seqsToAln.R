#' @title Creates a text alignment from a set of cigar strings
#'@description Creates a one-to-one text alignment of a set of cigar strings with respect
#'to the reference sequence by collapsing insertions and introducing gaps
#'across deletions.
#'
#'When genomic coordinates for the alignment start and the target region
#'are provided, aligned sequences are cropped to the target region
#'@author Helen Lindsay
#'@param cigar A list of cigar strings to align
#'@param dnaseq The set of sequences corresponding to the cigars, as Biostrings::DNAStrings
#'@param target The target region to return, as GRanges.  Sequences overlapping
#'the target region are trimmed to exactly match it.
#'@param del_char The character to represent deleted bases. Default "-"
#'@param aln_start Genomic start locations of aligned sequences. Should be
#'used in conjunction with target_start and target_end.
#'@param reverse_complement (Default: FALSE)
#'@param allow.partial Are alignments that do not span the target region allowed?
#'(Default: FALSE)
#'@return The sequences with insertions collapsed and deletions padded
#'@rdname seqsToAln
seqsToAln <- function(cigar, dnaseq, target, del_char = "-", 
                      aln_start = NULL, reverse_complement = FALSE,
                      allow_partial = FALSE){
    # Additional trimming is necessary because deletion operations may overhang
    # boundaries of target region
    
    if (as.character(GenomicRanges::strand(target)) == "*") strand <- "+"
  
    result <- GenomicAlignments::sequenceLayer(dnaseq, cigar, D.letter = del_char,
                                               N.letter = "N")
    if (! is.null(aln_start)){
      shifts <- aln_start - GenomicAlignments::start(target)
      result <- Biostrings::stackStrings(result, shift = shifts, from = 1,
                 to = width(target), Lpadding.letter = "+", Rpadding.letter = ".")
      if (! isTRUE(allow_partial) & isTRUE(any(grep("\\+|\\.", result)))){
        stop(paste("When allow_partial is FALSE, dnaseq to be",
                   "trimmed must span the target location"))
      }
    }

    if (isTRUE(reverse_complement)){
      result <- Biostrings::reverseComplement(result)
    }
    result <- as.character(result)
    result <- gsub("\\.", "<", gsub("\\+", ">", result))
    result
}
