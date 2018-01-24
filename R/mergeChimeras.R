#'@title mergeChimeras
#'@description Merges chimeric alignments where the individual segments
#' border an unmapped region (a long deletion).  If bases of the read are
#' mapped to both ends of the gap, the multimapped reads are only included
#' in the leftmost genomic segment. If there are more than max_unmapped
#' unmapped bases between the mapped bases, the read is not considered
#' mergeable.  Currently experimental and only tested with reads mapped by
#' bwa mem.
#'@param bam A GenomicAlignments::GAlignments object
#'@param chimera_idxs Indices of chimeric reads within bam
#'@param verbose  Should information about the number of mergeable
#' alignments be printed?  (Default: TRUE)
#'@param max_read_overlap  Maximum number of bases in a mergeable read
#'that are aligned to two genomic locations (Default: 10)
#'@param max_unmapped  Maximum number of bases in a mergeable read
#' that are unmapped and located between two mapped segments (Default: 4)
#'@param name Name of the sample, used when reporting verbose output.
#'@author Helen Lindsay
#'@return A list of the merged and unmerged chimeric alignments
#'@rdname mergeChimeras
mergeChimeras <- function(bam, chimera_idxs = NULL, verbose = TRUE,
                        max_read_overlap = 10, max_unmapped = 4, name = NULL){
    
    # max_unmapped refers to unmapped sections between chimeric segments
    
    # Case: Aligned regions overlap
    #  1-2-3-4-5
    #      3-4-5-6-7
    #
    # Case: Inversion:
    # 1-2-3-4-5
    #            9-8-7-6

    message(paste0("Caution: mergeChimeras assumes a sorted bam file\n",
                   "and has only been tested with bwa mem alignments!\n"))

    if (is.null(chimera_idxs)){
        chimera_idxs <- findChimeras(bam)
    }

    if (length(chimera_idxs) == 0){                          
        return(list(merged = bam, unmerged = GenomicAlignments::GAlignments()))
    }

    # Do all reads within a chimera map to the same chromosome?
    nms <- rle(names(bam)[chimera_idxs])
    nms_codes <- rep(1:length(nms$lengths), nms$lengths)
    sqs <- seqnames(bam)[chimera_idxs]
    sqs_codes <- rep(1:length(sqs@lengths), sqs@lengths)
    codes <- rle(paste(nms_codes, sqs_codes, sep = "."))
    one_chr <- rep(codes$lengths, codes$lengths) == rep(nms$lengths, nms$lengths)
    
    # And onto the same strand? (i.e. not inversion)
    strds <- strand(bam)[chimera_idxs]
    strd_rle <- rle(paste0(nms_codes, strds))
    same_strd <- rep(strd_rle$lengths, strd_rle$lengths) ==
                    rep(nms$lengths, nms$lengths)
  
    # Are single chr chimeras gaps? (start(n+1) > end(n))
    del_lns <- start(bam)[chimera_idxs[-1]] -
                           end(bam)[chimera_idxs[-length(chimera_idxs)]]
    is_after <- c(TRUE, del_lns > 0)
    change_pts <- cumsum(nms$lengths) + 1 # note starts from second and includes last
    change_pts <- c(1, change_pts[1:length(change_pts) -1])
    is_after[change_pts] <- TRUE
    codes <- rle(paste(nms_codes, is_after, sep = "."))
    has_genome_gap <- rep(codes$lengths, codes$lengths) ==
                          rep(nms$lengths, nms$lengths)

    # Is the same read segment used in multiple sections of a chimera?
    # For merge-able alignments, the sum of the widths of the aligned regions of read n-1
    # should be less than or equal to the first aligned base of read n wrt the full seq
    cigars <- cigar(bam)[chimera_idxs]
    qrng <- GenomicAlignments::cigarRangesAlongQuerySpace(cigars,
               before.hard.clipping = TRUE, ops = c("M", "I"))
    first_aligned <- min(start(qrng))
    last_aligned <- max(end(qrng))

    # Gaps not correct: if the first part of the read maps after the
    # second part of the read gap - however, here only care about +ve / -ve
    gaps <- c(1,  first_aligned[-1] - last_aligned[-length(last_aligned)])
    gaps[change_pts] <- 1
    read_mapped_uniq <- gaps > -1 * max_read_overlap # 
    read_covered <- gaps <= max_unmapped
       
    gap_codes <- rle(paste(read_mapped_uniq & read_covered, nms_codes, sep = "."))
    has_read_gap <- rep(gap_codes$lengths, gap_codes$lengths) ==
                    rep(nms$lengths, nms$lengths)

    # Rearrangements: where the first aligned base of the second segment is
    # earlier than the first 
    rearr <- first_aligned[-1] - first_aligned[-length(first_aligned)]
    rearr <- c(1, rearr)
    rearr[change_pts] <- 1
    rearr <- rearr < 0
    
    # Remove sections of a read that map to two genomic locations
    # - identify reads with negative gap = overlap between segments
    read_gaps <- c(0, first_aligned[-1] - last_aligned[-length(last_aligned)] - 1)
    read_gaps[change_pts] <- 0
    to_cut <- read_gaps < 0
    
    # Start constructing new cigars and test for gaps in cut sections
    ch_ends <- cumsum(nms$lengths)
    new_cigars <- GenomicAlignments::cigarNarrow(cigars, start = 1,
                                                 end = width(bam[chimera_idxs]))
    new_cigars <- as.character(new_cigars)
    
    first_range <- as.numeric(gsub('M.*', "", new_cigars[to_cut]))
    # Cut overlap off rightmost read, adjust genomic coordinates
    first_range <- first_range + read_gaps[to_cut]
      
    # Remove reads with indels in the segment to exclude.
    rm_nms <- inverse.rle(nms)[which(to_cut)[which(first_range < 0)]]
    indel_in_cut <- names(bam[chimera_idxs]) %in% rm_nms
  
    mergeable <- one_chr & same_strd & has_genome_gap &
                 has_read_gap & !rearr & ! indel_in_cut
    
    chimeras <- bam[chimera_idxs][mergeable]
    chimeras <- chimeras[!duplicated(names(chimeras))] # select the leftmost
    remaining <- bam[chimera_idxs][!mergeable]
    
    if (!any(mergeable)){                          
      return(list(merged = chimeras, unmerged = remaining))
    }
    
    if (isTRUE(verbose)){ 
      format_zero <- function(x) ifelse(is.nan(x), 0, x)
      nchm <- length(chimera_idxs)
      noc <- sum(one_chr == "TRUE")
      nss <- sum(one_chr & !same_strd == "TRUE")
      nrearr <- sum(one_chr & same_strd & rearr == "TRUE")
      ngdup <- sum(!has_genome_gap & one_chr & same_strd & !rearr == "TRUE")
      nrdup <- sum(has_genome_gap & one_chr & same_strd &
                   !read_mapped_uniq & !rearr == "TRUE")
      nrumap <- sum(has_genome_gap & one_chr & same_strd &
                   !read_covered & !rearr == "TRUE")
      mrg <- sum(one_chr & same_strd & has_genome_gap & has_read_gap & ! rearr == "TRUE")
      if (! is.null(name)) cat(sprintf("Chimera statistics for %s:\n", name))  
      cat(sprintf(paste0("%s (%.2f%%) chimeras in %s reads\n",
      "  %s (%.2f%%) map to the same chromosome\n", 
      "    %s (%.2f%%) map to different strands (inversions)\n",
      "    %s (%.2f%%) rearrangements (end of read maps before start)\n",
      "    %s (%.2f%%) genomic duplications (different read locs mapped",
           " to same genomic loc)\n",
      "    %s (%.2f%%) read duplications (different genomic locs mapped",
           " to same read loc)\n",
      "    %s (%.2f%%) unmapped segments between mapped segments\n",
      "    %s (%.2f%%) are long gaps\n\n"),
      nchm, format_zero(nchm/length(bam)*100), length(bam),
      noc, format_zero(noc/nchm*100),
      nss, format_zero(nss/noc*100),
      nrearr, format_zero(nrearr/noc*100),
      ngdup, format_zero(ngdup/noc*100),
      nrdup, format_zero(nrdup/noc*100),
      nrumap, format_zero(nrumap/noc*100),
      mrg, format_zero(mrg/noc*100))) 
    }
    
    # For first member of chimera: get everything up to and including the last M
    new_cigars[change_pts] <- gsub("(^.*M)[0-9]+[HS]","\\1", cigars[change_pts])
    # For last member of chimera: get everything except for clipping at the start
    new_cigars[ch_ends] <- gsub("H","S",gsub("^[0-9]+[HS](.*)", "\\1", cigars[ch_ends]))
    
    new_cigars[to_cut] <- paste0(first_range, gsub('[0-9]+(M.*)', "\\1",
                                                   new_cigars[to_cut]))
    
    # Stick cigars together padding segments with deletions
    
    ch_gstarts <- start(bam)[chimera_idxs]
    ch_gends <- end(bam)[chimera_idxs]
    ### NEED TO ADJUST ENDS FOR READS TRIMMED
    ch_gstarts[to_cut] <- ch_gstarts[to_cut] - read_gaps[to_cut]

    ggaps <- c(sprintf("%sD", -1*(ch_gends[-length(chimera_idxs)] - ch_gstarts[-1] +1)),0)
    ggaps[ch_ends] <- ""
    new_cigars <- paste0(new_cigars, ggaps)

    new_cigars <- aggregate(new_cigars, list(nms_codes), FUN = paste0, collapse = "")$x

    # WARNING: assumption:  The primary alignment is not hard clipped
    new_cigars <- gsub("H", "S", new_cigars)    

    not_suppl <- !bitwAnd(mcols(bam[chimera_idxs])$flag, 2048)
    sqs <- mcols(bam[chimera_idxs])$seq[not_suppl]
    names(sqs) <- names(bam[chimera_idxs])[not_suppl]
    sqs <- unname(sqs[names(chimeras)])
    keep_chs <- mergeable[!duplicated(nms_codes)]
   
    # WARNING: if there is an unmapped section of the read < max_unmapped,
    # it is deleted 
    if (any(read_gaps > 0 )){
      keep_chs <- mergeable[!duplicated(nms_codes)]
      ch_qr <- sum(relist(elementNROWS(qrng), PartitioningByWidth(nms$lengths)))
      qrng <- relist(unlist(qrng), PartitioningByWidth(ch_qr))
      #qrng <- relist(unlist(qrng), IRanges::PartitioningByWidth(nms$lengths))
      gps <- Biostrings::gaps(qrng)[keep_chs]
      sqs <- Biostrings::replaceAt(sqs, gps, "")
    }
    
    galns <- GAlignments(names = names(chimeras),
                    seqnames = as.factor(seqnames(chimeras)),
                    pos = start(chimeras), cigar = new_cigars[keep_chs],
                    strand = strand(chimeras),
                    seqlengths = seqlengths(chimeras),
                    seq = sqs, flag = 0)

    return(list(merged = galns, unmerged = remaining))

}
