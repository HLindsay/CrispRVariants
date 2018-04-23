# .formatVarLabels -----
#'@title formatVarLabels
#'@description Internal CrispRVariants function for creating allele
#'labels given variants and positions
#'@param grl  (GRangesList) A GRangesList of variants.  The position
#'of the variants is used in labels 
#'@param labels (character(n)) A vector of labels for each variant.
#'In CrispRVariants, this is the size and type of the variant, e.g. "9D"
#'for a 9 bp deletion.
#'@param position One of "start" and "end".  Determines whether the
#'start or the end coordinate is used when labeling variants.
#'@param genome.to.pos Optional named vector for transforming 
#'variant coordinates into another coordinate system (Default: NULL)
#'@param pos.to.lab (character(1)) Character to join positions and labels
#'(Default: ":", e.g. -3:9D)
#'@param as.string Should individual variant labels be pasted into
#'a single comma separated string when one alignment has multiple variants?
#' (Default: TRUE)
#'@rdname alleleLabelsHelpers
.formatVarLabels <- function(grl, labels, position = c("start", "end"),
                             genome.to.pos = NULL, pos.to.lab = ":",
                             as.string = TRUE){
  
    # Select either the start or end of the operation as the op position
    position <- match.arg(position)
    
    if (position == "start"){
      glocs <- start(grl)
    } else {
      glocs <- end(grl)
    } 
    glocs <- as.character(unlist(glocs))
    
    # Map position to another value if a mapping is given
    if (! is.null(genome.to.pos)){
      glocs <- genome.to.pos[glocs]
    }
    
    # Add location to cigar operations
    lns <- lengths(grl)
    op_labs <- paste(glocs, unlist(labels), sep = pos.to.lab)
    op_labs <- relist(op_labs, IRanges::PartitioningByWidth(lns))
    
    # Reverse each member of the list if using the end as the position 
    if (position == "end"){
      op_labs <- rev(relist(rev(unlist(op_labs)), rev(op_labs)))
    }
    
    # Collapse into a single string
    if (isTRUE(as.string)){
      op_labs <- paste(op_labs, collapse = ",")
    }
    
    op_labs
} # -----

# indelLabels -----
#'@title indelLabels
#'@description Makes allele labels for insertion / deletion variants
#'@author Helen Lindsay
#'@param alns (GAligments) aligned reads for finding variants
#'@param rc Should the variants be displayed with respect to the 
#'negative strand? (Default: FALSE)
#'@param genome.to.pos A vector with names being genomic locations and
#'values being positions to use in labels (Default: NULL)
#'@param keep.ops CIGAR operations to remain in the variant label
#'(usually indels)
#'@param regions  IRanges(k) Regions for counting insertions and
#'deletions.  Insertions on the right border are not counted.
#'@param ... extra formatting arguments
#'@return A vector of labels for alns
indelLabels = function(alns, rc = FALSE, genome.to.pos = NULL,
                       keep.ops = c("I","D","N"), regions = NULL, 
                       ...){
    
    # MAKE CONSISTENT: USE GENOMIC OR RELATIVE COORDS
    cigs <- GenomicAlignments::cigar(alns)
    op_locs <- selectOps(cigs, ops = keep.ops,
                         op.regions = regions,
                         pos = start(alns))
    
    position <- ifelse(isTRUE(rc), "end", "start")
    op_labs <- .formatVarLabels(op_locs$op_rngs, op_locs$op_labels,
                                position = position,
                                genome.to.pos = genome.to.pos, ...)
    
    op_labs
} # -----

# .resolveRefStrand -----
# If the reference provided maps to the negative strand, use
# the reverse complement when finding mismatches using the bam file
.resolveRefStrand <- function(strand, reference){
    if (strand == "*") strand <- "+"
    if (! strand %in% c("+","-","*")){
      stop('Strand should be either "+" or "-"')
    }
    # If reference is negative, use the positive to match bam format
    if (strand == "-"){
      reference <- Biostrings::reverseComplement(reference)
    }
    reference
} # -----

# findMismatches -----
# As only a single reference is provided, assume that alns start and end 
# at the same position
#'@title findMismatches
#'@description Assume that the reference may be on the negative strand and
#'regions are given with respect to the reference sequence.
#'@param alns A GAlignments object, where the aligned sequences should span the
#'reference sequence
#'@param ref.start (numeric(1)) The genomic start position of the reference sequence
#'@param ref.seq A DNAString object, the sequence for comparison when checking
#'for mismatches.  The sequence does not necessarily have to match the mapping
#'reference sequence.  Must span all regions if regions are provided.
#'@param regions A GRanges object, regions to check for mismatches with coordinates
#'relative to the reference sequence 
#'@param strand One of "+", "-"
#'@param min.pct (numeric(1), between 0 and 100) Only return SNVs that occur
#'at in least min.pct % of reads.  This tests the occurrence of the exact base
#'change, not any change at a position. 
#'@return A data frame of sequence indices, genomic position of mismatch
#'and mismatch base
#'@rdname alleleLabelsHelpers
.findMismatches <- function(alns, ref.seq, ref.start, regions = NULL,
                            strand = "+", min.pct = 0){
  
    # Using a single reference, rather than a reference per region,
    # as CrispRVariants considers variants within individual reads 
    # Sequences with deletions do not count towards min.pct for each column
 
    rwdth = nchar(ref.seq)
    
    # Check that regions are within the reference range
    if (! is.null(regions)){
      if (! all(regions %within% IRanges(1, rwdth))){
        stop("Region for counting mismatches longer than reference sequence")
      }
      # Convert regions to positive strand if they are negative
      if (as.character(strand) == "-"){
        temp <- structure(rev(1:nchar(ref.seq)), names = 1:nchar(ref.seq))
        strt <- temp[as.character(start(regions))]
        end <- temp[as.character(end(regions))]
        regions <- IRanges(end, strt)
      }
    }
    
    # As the alignments are to the +ve strand, get the +ve strand reference
    strand <- as.character(strand)
    ref.seq <- .resolveRefStrand(strand, ref.seq)
    
    # Get sequences aligned to the reference region, i.e. insertions ignored
    sqs <- sequenceLayer(mcols(alns)$seq, cigar(alns))
    if (! all(start(alns) == ref.start)){
      sqs <- stackStrings(sqs, from = 1,
                          shift = start(alns) - ref.start,
                          to = rwdth,
                          Lpadding.letter = "N",
                          Rpadding.letter = "N")
    }

    # Positions will be used for the labels
    posns <- ref.start:(ref.start + rwdth - 1)
    names(posns) <- 1:rwdth
    
    # HERE TO LOOP OVER REGIONS
    # unlist(unname(GAlignmentsList(alns[[1]])))
    
    
    
    # Note: using replaceAt as looking at single base positions
    if (! is.null(regions)){
      del_regions <- gaps(regions, start = 1, end = rwdth)
      sqs <- Biostrings::replaceAt(sqs, del_regions, "")
      posns <- posns[as.logical(coverage(regions, width = rwdth))]
      ref.seq <- replaceAt(ref.seq, del_regions, "")
    } 
    
    sqsm <- as.matrix(sqs)
    ref <- as.matrix(ref.seq)
    eq_ref <- t(t(sqsm) == as.vector(ref))
    
    is_base <- matrix(sqsm %in% c("A","C","T","G"), nrow = nrow(sqsm))
    snv <- is_base & ! eq_ref
    wh.snv <- which(snv, arr.ind = TRUE)
    if (strand == "-"){
      rc_base <- as.character(reverseComplement(DNAStringSet(sqsm[snv])))
      result <- data.frame(seq = wh.snv[,1],
                           pos = posns[wh.snv[,2]],
                           base = rc_base)
    } else {
      result <- data.frame(seq = wh.snv[,1],
                           pos = posns[wh.snv[,2]],
                           base = sqsm[snv])      
    }

    result <- result[order(result$seq), ]
    
    if (min.pct > 0){
      pcts <- prop.table(consensusMatrix(sqs)[c("A","C","T","G"),], 2) * 100
      # Note: keep includes reference but will not appear in result table
      keep <- which(pcts >= min.pct, arr.ind = TRUE)
      keep <- paste(posns[keep[,2]], rownames(keep))
      pos_base <- paste(result$pos, result$base)
      result <- result[pos_base %in% keep,]
    }
    
    result
  
} # -----

# nonindelLabels -----
#'@title nonindelLabels
#'@description Make variant labels for variants without an insertion or deletion
#'@param alns A GAlignments object, where the aligned sequences should span the
#'reference sequence
#'@param target (GRanges(1)) The region for counting mismatches 
#'@param ref.seq A DNAString object, the sequence for comparison when checking
#'for mismatches.  The sequence does not necessarily have to match the mapping
#'reference sequence.  Must span all regions if regions are provided.
#'@param regions A GRanges object, regions to check for mismatches with coordinates
#'relative to the reference sequence 
#'@param min.pct (numeric(1), between 0 and 100) Only return SNVs that occur
#'at in least min.pct % of reads.  This tests the occurrence of the exact base
#'change, not any change at a position. 
#'@param mismatch.label (character(1)) Label to append to the start of mismatch
#'strings, if returning as a single string (Default: "SNV:")
#'@param genome.to.pos Optional named vector for transforming 
#'variant coordinates into another coordinate system (Default: NULL)
#'@param as.string Should individual variant labels be pasted into
#'a single comma separated string when one alignment has multiple variants?
#' (Default: TRUE)
#'@return A data frame of sequence indices, genomic position of mismatch
#'and mismatch base
mismatchLabels <- function(alns, target, ref.seq, 
                           regions = NULL, min.pct = 0,
                           mismatch.label = "SNV",
                           genome.to.pos = NULL,
                           as.string = TRUE){
    # Multiple counting rules
    # Check rules don't overlap
    # Each rule may have from bases?  and to bases
  

    mm <- .findMismatches(alns, ref.seq = ref.seq,
                          ref.start = start(target),
                          regions = regions,
                          strand = as.character(strand(target)),
                          min.pct = min.pct)
    
    allele_labels <- rep("", length(alns))
    if (nrow(mm) == 0){
      return(allele_labels)
    }
    
    if (! is.null(genome.to.pos)){
      mm$pos <- genome.to.pos[as.character(mm$pos)]
    } 
    
    mm_labs <- paste0(mm$pos, mm$base)
    temp <- split(mm_labs, mm$seq)
    temp <- relist(unlist(unname(temp)),
                   PartitioningByWidth(lengths(temp)))
    
    if (isTRUE(as.string)){
      temp <- paste(temp, collapse = ",")
      tnms <- names(temp)
      temp <- paste(mismatch.label, temp, sep = ":")
      names(temp) <- tnms
    }
    
    allele_labels[as.numeric(names(temp))] <- unname(temp)
    allele_labels
  
} # -----


matchLabels <- function(alns, target, match.label = "No variant"){
    # Label partial alignments differently to perfect matches
    # RENAMING PARTIAL, NONMATCHING ALNS NOT FINISHED

    spans <- start(alns) <= start(target) & end(alns) >= end(target)
    nsp_cigar <- cigar(alns)[! spans] 
    rngs <- unname(ranges(alns))[! spans]
    splt <- split(rngs, cigar(alns)[! spans])
    nunq <- splt(lengths(splt) > 1)

    aln_labels <- rep(match.label, length(alns))
    aln_labels[! spans] <- cigar(alns)[! spans]
}


# .defaultCigarLabels -----
# temporary function until CrisprRun class is removed
.defaultCigarLabels = function(cset, renumbered = TRUE,
                               match_label = "no variant",
                               mismatch_label = "SNV",
                               split.snv = TRUE,
                               upstream.snv = 8, downstream.snv = 6,
                               bpparam = BiocParallel::SerialParam()){
  
  # Are the match and mismatch labels already set in the cset?
  if (all(lengths(alns(cset))) == 0) return(NULL)
  
  g_to_t <- NULL
  target.loc = cset$pars$target.loc
  target_start = start(cset$target)
  target_end = end(cset$target)
  rc = cset$pars$rc
  ref = cset$ref
  
  if (isTRUE(renumbered)){
    if (any(is.na(c(target_start, target_end, rc)))){
      stop("Must specify target.loc (cut site), target_start,
           target_end and rc for renumbering")
    }

    g_to_t <- cset$.genomeToTargetLocs(target.loc, target_start, target_end, 
                     as.character(strand(cset$target)) %in% c("+", "*"))
  }
  
  # rc means "display on negative strand"
  # Reversing ranges is dealt with in .findMismatches

  # This section is slow
  cig_by_run <- BiocParallel::bplapply(cset$crispr_runs,
                         function(crun)  crun$getCigarLabels(
                           target = cset$target,
                           target.loc = target.loc, #cut.site,
                           genome_to_target = g_to_t,
                           ref = ref,
                           separate.snv = split.snv,
                           match.label = cset$pars$match_label,
                           mismatch.label = cset$pars$mismatch_label,
                           rc = rc, upstream = upstream.snv,
                           downstream = downstream.snv), BPPARAM = bpparam)
  
  cig_by_run
} # -----
