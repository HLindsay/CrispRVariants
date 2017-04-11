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
#'@param allow_partial Are alignments that do not span the target region allowed?
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


#'@title seqsToPartialAln
#'@description Given a character vector of pairwise alignments and a
#' region to display, trims alignments to the display region
#'@author Helen Lindsay
#'@param cigars cigar strings corresponding to alignments
#'@param starts Genomic start coordinates of the alignments
#'@param alns Character vector of pairwise alignments, with insertions removed
#'@param reference Reference sequence
#'@param target The complete region spanned by the alignments (GRanges)
#'@param keep Region to display (IRanges or GRanges)
#'@param xbreaks Mapping between coordinates and plot locations.
#'  If supplied, plot locations are adjusted for gaos (Default: NULL)
#'@param xticks Locations for x ticks in plot (Default: NULL)
seqsToPartialAln <- function(cigars, starts, alns, reference, target,
                          keep, xbreaks = NULL, xticks = NULL){

    # Get the region to delete
    if (class(keep) == "GRanges") keep <- ranges(keep)
    to_delete <- Biostrings::gaps(keep, start = start(target), end = end(target))

    # Find the operations that overlap the region to be removed
    # Future - choose a region to display, presently just checking for I
    rr <- cigarRangesAlongReferenceSpace(cigars,
                with.ops = TRUE, ops = c("D","I"))
    rr <- unlist(rr)
    
    # Insertions locations will need to be shifted to match the deleted
    # segments, this isn't implemented yet
    if ("I" %in% names(rr)){
        stop("Displaying partial alignments containing insertions is not
             implemented yet")
    }
    
    qkeep <- shift(keep, 1-start(target))
    result <- extractAt(as(alns, "XStringSet"), qkeep)
    ref_result <- extractAt(as(ref, "XString"), qkeep)    

    ndeleted <- width(Biostrings::gaps(keep))
    if (length(ndeleted) > 1) {
        stop("Displaying more than two segments is not implemented yet")
    }

    ref_paste <- sprintf("/ %s /", ndeleted)    
    aln_paste <- sprintf("/ %s /", paste(rep(".", nchar(ndeleted)), collapse = ""))
    
    aln_result <- lapply(result, paste, collapse = aln_paste)
    ref_result <- paste(ref_result, collapse = ref_paste)
    alns <- as.character(aln_result)
    names(alns) <- names(aln_result)    

    if (! is.null(xbreaks)){
        
    }

    return(list(alns = alns, ref = ref_result))
}

