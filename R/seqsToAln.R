# seqsToAln -----
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
  
    result <- GenomicAlignments::sequenceLayer(dnaseq, cigar,
                                               D.letter = del_char,
                                               N.letter = "N")
    if (! is.null(aln_start)){
      shifts <- aln_start - GenomicAlignments::start(target)
      result <- Biostrings::stackStrings(result, shift = shifts,
                                         from = 1, to = width(target),
                                         Lpadding.letter = "+",
                                         Rpadding.letter = ".")
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
} # ------

# selectAlnRegions ------
#'@title selectAlnRegions
#'@description Given a character vector of pairwise alignments and a
#' region to display, trims alignments to the display regions, joined by
#' a separator "join".  Alignments should be equal length, e.g.
#' created by seqsToAln
#'@author Helen Lindsay
#'@param alns Character vector of pairwise alignments, with insertions removed
#'@param reference Reference sequence
#'@param keep Region to display, relative to the target region, i.e.
#'not genomic coords (IRanges or GRanges)
#'@param join character(1) String used for joining alignment segments.
#'Can accept a placeholder to fill in the number of bases deleted with
#'"%s", e.g. the default "/ %s /" will appear as "/ 3 /" if 3 bases are
#'deleted
#'@param border.gaps (logical(1)) Should bases deleted from the borders
#'be shown?  (Default: FALSE)
#'@return A list of the truncated alignments (alns) and reference (ref)
#'@rdname seqsToAln
selectAlnRegions <- function(alns, reference, target,
                             keep, join = "/ %s /",
                             border.gaps = FALSE){

    # Potential changes:
    # record how many bases deleted from each allele
    # (e.g. to make clear if deleted region was a gap or an insertion)
  
    # make keep relative to target if not already
    keep <- .checkRelativeLocs(target, keep)
    to_delete <- .invertKeepRanges(target, keep)
    if (! class(to_delete) == "IRanges") return()
    
    result <- Biostrings::extractAt(as(alns, "XStringSet"), keep)
    ref_result <- Biostrings::extractAt(as(reference, "XString"), keep)    

    ndeleted <- width(to_delete)
    ref_paste <- sprintf(join, ndeleted)    
    aln_paste <- sprintf(join, strrep(".", nchar(ndeleted)))
    seg_joins <- matrix(c(ref_paste, aln_paste), ncol = 2)
    
    # Adjust to show or hide operations at the borders 
    del_left <- min(start(keep)) > min(start(to_delete))
    del_right <- max(end(keep)) < max(end(to_delete))

    if (isTRUE(del_left) & ! isTRUE(border.gaps)){
      # The first gap is not included in the displayed plot
      seg_joins[1,] <- ""
    } 
    if (isTRUE(del_right) & ! isTRUE(border.gaps)){
      seg_joins[nrow(seg_joins),] <- ""
    }
    if (! isTRUE(del_left)){ seg_joins <- rbind("", seg_joins) }
    if (! isTRUE(del_right)){ seg_joins <- rbind(seg_joins, "") }
    
    # Offsets refer to number of characters added into plot
    offsets <- nchar(seg_joins[1:length(keep),1])
    
    keep_idxs <- coverage(keep, width = width(target))
    offset_idxs <- .offsetIndices(nchar(ref_result), offsets)
    keep_idxs[keep_idxs == 1] <- offset_idxs
    keep_idxs <- as.numeric(keep_idxs)
    
    # Joins include border elements, empty if border.gaps = FALSE
    paste_segs <- function(joins, segs){
       paste(.intersperse(joins, segs), collapse = "")
    }
    
    ref_return <- paste_segs(seg_joins[,1], as.character(ref_result))
    aln_return <- unlist(lapply(result, function(x){ 
      paste_segs(seg_joins[,2], as.character(x))
    })) 
    
    return(list(alns = aln_return, ref = ref_return, keep = keep,
                nchar_join = offsets, indices = keep_idxs))
} # ------

# .invertKeepRanges -----
#'@title Helper functions for selectAlnRegions
#'@description (.invertKeepRanges) Internal CrispRVariants function used
#'by seqsToPartialAlns for checking arguments and getting region to delete.
#'Returns FALSE if no region to delete found, or region to be deleted is
#'entire target.
#'@param target The complete region spanned by the alignments (GRanges)
#'@param keep Region to display, relative to the target region, i.e.
#'not genomic coords (IRanges or GRanges)
#'@return Gaps between keep (IRanges), or FALSE if no gap ranges found
#'@rdname selectAlnRegionsHelpers
.invertKeepRanges <- function(target, keep){
    if (! class(keep) %in% c("IRanges", "GRanges")){
      stop("Ranges to keep should be class IRanges or GRanges")
    }
    
    # Get the region to delete
    keepr <- keep
    if (class(keep) == "GRanges") keepr <- ranges(keep)
    to_delete <- Biostrings::gaps(keepr, start = 1, end = width(target))
    
    # If the region to be deleted is the target region,
    # check if it's in genomic coordinates
    if (identical(to_delete, IRanges(1, width(target)))){
      # Check if all regions to keep are contained within target
      warning("Region to delete is identical to the target region, \n",
              "or was of class IRanges and could not be shifted\n",
              "to relative coordinates. Returning")
      return(FALSE) 
    }
    
    # If there is nothing to delete, warn and return
    if (length(to_delete) == 0){
      warning("Nothing to delete.  Returning")
      return(FALSE)
    } 
    
    return(to_delete)
} # -----

# .checkRelativeLocs -----
#'@description (.checkRelativeLocs) Shift keep to start at 1 if it is
#'within the target
#'@rdname selectAlnRegionsHelpers
.checkRelativeLocs <- function(target, keep){
    # Can only check for contained ranges if keep is a GRanges object
    # (Or both are IRanges)
    if (! (class(keep) == "GRanges" & class(target) == "GRanges")) return(keep)
    ovs <- queryHits(findOverlaps(keep, target, type = "within"))
    if (length(setdiff(ovs, seq_along(keep))) == 0){
      keep <- shift(keep, 1-start(target))
    }
    ranges(keep)
} # -----

# .adjustRelativeInsLocs -----
#'@description (.adjustRelativeInsLocs) Internal CrispRVariants function
#'for shifting insertion locations relative to the target region when
#'removing a segment of the alignments.  Note that this function does not
#'do input checking but assumes this has been done upstream. Insertions at
#'the left border of a gap region are removed.
#'@param starts numeric(n) Insertion locations. When plotting,
#'the insertion symbol appears at the left border of the
#'start location.
#'@param gap_nchars character(n) Number of letters added to when joining
#'segments before each region in keep.  If first base of keep is 1,
#'the first entry of gap_nchars should be 0.
#'@return insertion_sites (data.frame) with modified start column
#'@rdname selectAlnRegionsHelpers
.adjustRelativeInsLocs <- function(target, keep, starts, gap_nchars){
    if ( length(gap_nchars) > length(keep) + 2 ){
      # Can be at most 2 more gaps than kept (boundaries)
      stop("gap_nchars should be the number of characters used\n",
           "to separate segments in keep when creating alignments.\n",
           "length(gap_nchars) should be at most two more than length(keep)") 
    }
    # This is the same as the indices in selectAlnRegions except that
    # left borders of gaps are included. 
  
    # As the insertion is placed at the leftmost border of a nucleotide,
    # it is still possible to display insertions at the left border of
    # the region(s) to remove
    not_rightmost <- ! end(keep) == width(target)
    end(keep[not_rightmost]) <- end(keep)[not_rightmost] + 1
    keep_cov <- as.logical(coverage(keep, width = width(target)))
    new_idxs <- rep(NA, length(keep_cov))
    new_idxs[keep_cov] <- seq_along(which(keep_cov))
    new_starts <- new_idxs[starts]
    
    # Offset the new insertion locations by the number of bases
    # added when joining sequences
    all_offsets <- rep(gap_nchars, width(keep))
    to_display <- ! is.na(new_starts)
    new_starts[to_display] <- new_starts[to_display] + all_offsets[new_starts[to_display]] 
    new_starts
    
} # -----

# .offsetIndices -----
#' @description (.offsetIndices) Get indices of a vector grouped by "x",
#' cumulatively adding offsets to each group according to "offset"
#' @param x  A vector of group lengths
#' @param offset A vector of offset lengths matching x
#' @author Helen Lindsay
#' @examples 
#' CrispRVariants:::.offsetIndices(rep(2,5), c(0:4)*10)
#' @rdname selectAlnRegionsHelpers
.offsetIndices <- function(x, offset){
  stopifnot(length(x) == length(offset))
  indices <- seq_len(sum(x))
  indices <- indices + rep(cumsum(offset), x)
  indices
} # -----

