#'@title CrisprRun class
#'@description A ReferenceClass container for a single sample of alignments narrowed
#'to a target region.  Typically CrisprRun objects will not be accessed directly,
#'but if necessary via a CrisprSet class which contains a list of CrisprRun objects.  
#'Note that the CrispRVariants plotting functions don't work on CrisprRun objects.
#'@param bam a GAlignments object containing (narrowed) alignments to the target region.
#' Filtering of the bam should generally be done before initialising a CrisprRun object
#'@param target The target location, a GRanges object
#'@param genome.ranges A GRangesList of genomic coordinates for the cigar operations.
#' If bam is a standard GAlignments object, this is equivalent to
#' cigarRangesAlongReferenceSpace + start(bam)
#'@param rc (reverse complement)  Should the alignments be reverse complemented,
#'i.e. displayed with respect to the negative strand?  (Default: FALSE)
#'@param name A name for this set of reads, used in plots if present (Default: NULL)
#'@param chimeras Off-target chimeric alignments not in bam.  (Default: empty)
#'@param verbose Print information about initialisation progress (Default: TRUE)
#'@field alns A GAlignments object containing the narrowed reads.  Note that if the alignments
#'are represented with respect to the reverse strand, the "start" remains with repect to the
#'forward strand, whilst the cigar and the sequence are reverse complemented.
#'@field name The name of the sample
#'@field cigar_labels A vector of labels for the reads, based on the cigar strings,
#'optionally renumbered with respect to a new zero point (e.g. the cut site) and
#'shortened to only insertion and deletion locations.
#'Set at initialisation of a CrisprSet object, but not at
#'initialisation of a CrisprRun object.
#'@field chimeras Chimeric, off-target alignments corresponding to alignments in alns
#'@seealso \code{\link[CrispRVariants]{CrisprSet}}
#'@author Helen Lindsay
#'@examples
#'# readsToTarget with signature("GAlignments", "GRanges") returns a 
#'# CrisprRun object
#'
#'bam_fname <- system.file("extdata", "gol_F1_clutch_1_embryo_1_s.bam",
#'                         package = "CrispRVariants")
#'param <- Rsamtools::ScanBamParam(what = c("seq", "flag"))
#'alns <- GenomicAlignments::readGAlignments(bam_fname, param = param, 
#'          use.names = TRUE)
#'
#'reference <- Biostrings::DNAString("GGTCTCTCGCAGGATGTTGCTGG")
#'gd <- GenomicRanges::GRanges("18", IRanges::IRanges(4647377, 4647399), strand = "+")
#'
#'crispr_run <- readsToTarget(alns, target = gd, reference = reference,
#'                            name = "Sample name", target.loc = 17)
#'
#'# Alternatively, CrisprRun objects can be accessed from a CrisprSet object
#'# e.g. crispr_set$crispr_runs[[1]]
#'@export CrisprRun
#'@exportClass CrisprRun
CrisprRun = setRefClass(
  Class = "CrisprRun",
  fields = c(alns = "GAlignments",
             name = "character",
             cigar_labels = "character",
             chimeras = "GAlignments",
             chimera_combs = "data.frame")
)

CrisprRun$methods(
  initialize = function(bam, target, rc = FALSE, name = NULL,
                        chimeras = GenomicAlignments::GAlignments(), verbose = TRUE){
    
    # If no name is provided, use the coordinates
    if (is.null(name)){
        name <<- sprintf("%s:%s-%s", seqnames(target), start(target), end(target))
    } else { name <<- name}

    if (isTRUE(verbose)) message(sprintf("\nInitialising CrisprRun %s\n", .self$name))

    alns <<- bam
    chimeras <<- chimeras
    chimera_combs <<- .self$.splitChimeras()

    if (length(bam) == 0) { return() }

    # recalculate genome.ranges in case of keep_unpaired = FALSE, not tested
    .self$getInsertionSeqs(target = target)
  },

  show = function(){
    print(c(class(.self), sprintf("CrisprRun object named %s, with %s on target alignments.",
                                  .self$name, length(.self$alns)), .self$alns))
  },

#  removeSeqs = function(idxs){
#'
#Description:
#  Remove sequences from a CrisprRun object and from the internal CrisprRun
#  fields that store insertion locations for plotting.
#
#Input parameters:
#  idxs:     Indexes of reads to remove'
#
#    # note insertions table is not updated
#    ins_key_idxs <- which(names(.self$ins_key) %in% idxs)
#
#    if (length(ins_key_idxs) > 0){
#      .self$field("ins_key", .self$ins_key[-ins_key_idxs])
#    }
#
#    # Insertion key refers to indexs with in the alignments, these have shifted after
#    # filtering.  Note - would need to shift insertions idxs if insertions is also updated
#    subtract <- rep(0, length(.self$alns))
#    subtract[idxs] <- 1
#    subtract <- cumsum(subtract)
#
#    temp <- .self$ins_key
#    nm_as_num <- as.numeric(names(temp))
#    names(temp) <- nm_as_num - subtract[nm_as_num]
#    .self$field("ins_key", temp)
#
#    # Remove the extra sequences from the chimeras
#    rm_by_nm <- names(alns)[idxs]
#    ch_to_keep <- !(names(.self$chimeras) %in% rm_by_nm)
#    .self$field("chimeras", .self$chimeras[ch_to_keep])
#
#    .self$field("alns", .self$alns[-idxs])
#    .self$field("cigar_labels", .self$cigar_labels[-idxs])
#  },

  getInsertionSeqs = function(target){
' 
Description:
  Set the "insertions" field - a table of the locations of insertions,
  and the "ins_key" field which relates sequences indices to the insertions 
  they contain
Input parameters:
'
    # Note that the start of a ref_ranges insertion is its genomic end (rightmost base)

    genome_ranges <- cigarRangesAlongReferenceSpace(cigar(.self$alns),
                                     pos = start(.self$alns), ops = "I")
    ref_ranges <- cigarRangesAlongReferenceSpace(cigar(.self$alns), ops = "I")
    idxs <- rep(seq_along(ref_ranges), lengths(ref_ranges))
    tseqs <- as.character(mcols(.self$alns)$seq)[idxs]

    if (length(tseqs) == 0) {
      return(NULL)
    }

    query_ranges <- cigarRangesAlongQuerySpace(cigar(.self$alns), ops = "I")
    qranges <- unlist(query_ranges)

    ins_seqs <- as.character(subseq(tseqs, start(qranges), end(qranges)))
    ins_starts <- start(unlist(shift(genome_ranges, 1 - start(target))))
    genomic_starts <- unlist(start(genome_ranges)) -1 # -1 for leftmost base
    
    df <- data.frame(start = ins_starts, seq = ins_seqs,
            genomic_start = genomic_starts,
            label = .self$cigar_labels[idxs], idxs)
    df$seq <- as.character(df$seq)
    agg <- aggregate(df$idx, by = as.list(df[,c(1:4)]), c)
    agg$count <- lengths(agg$x)
    colnames(agg) <- c("start", "seq", "genomic_start",
                       "cigar_label", "idxs", "count")
    return(agg)
  },

  .checkNonempty = function(){
    if (length(.self$alns) == 0){
      message("No on target alignments")
      return(FALSE)
    }
    TRUE
  },

  .splitChimeras = function(){
    splits <- split(cigar(.self$chimeras), names(.self$chimeras))
    if (length(splits) == 0) return(data.frame())
    combination <- sapply(splits, base::paste, collapse=";")
    tt <- as.data.frame(table(combination))
    tt <- tt[order(tt$Freq, decreasing = TRUE),]
    tt
  },


  getCigarLabels = function(target.loc, genome_to_target, ref,
                             separate.snv = TRUE, match.label = "no variant",
                             mismatch.label = "SNV",
                             rc = FALSE, keep.ops = c("I","D","N"), upstream = 8,
                            downstream = min(5, width(ref) - cut_site)){
    '
Description:
  Sets the "cig_labels" field, returns the cigar labels.

Input parameters:
  target.loc:       The location of the cut site with respect to the target
  genome_to_target: A vector with names being genomic locations and values
                    being locations with respect to the cut site
  separate.snv:     Should single nucleotide variants be called?  
                    (Default: TRUE)
  match.label:      Label for non-variant reads (Default: no variant)
  mismatch.label:   Label for single nucleotide variants (Default: SNV)
  rc:               Should the variants be displayed with respect to the 
                    negative strand? (Default: FALSE)
  keep.ops:         CIGAR operations to remain in the variant label
                    (usually indels)
  upstream:         distance upstream of the cut site to call SNVs
  downstream:       distance downstream of the cut site to call SNVs'
    
    cigs <- cigar(.self$alns)
    wdths <- explodeCigarOpLengths(cigs)
    ops <- explodeCigarOps(cigs)
    temp <- CharacterList(relist(paste0(unlist(wdths), unlist(ops)), wdths))
    ops <- CharacterList(ops)
    keep <- ops %in% keep.ops

    rranges <- cigarRangesAlongReferenceSpace(cigs, pos = start(.self$alns))[keep]
    # here get snvs, add snv ops

    if (isTRUE(rc)){
      glocs <- end(rranges)
    } else {
      glocs <- start(rranges)
    }
    
    temp <- paste(genome_to_target[as.character(unlist(glocs))],
                  unlist(temp[keep]), sep = ":")
        
    temp <- as.list(relist(temp, IRanges::PartitioningByEnd(cumsum(sum(keep)))))
          
    complex <- sum(keep) > 1

    if (isTRUE(rc)){
      temp[complex] <- sapply(temp[complex], function(x) paste(rev(x), collapse = ","))
    } else {
      temp[complex] <- sapply(temp[complex], base::paste, collapse = ",")
    }
    
    spans <- width(.self$alns) == nchar(ref) 
    temp[sum(keep) == 0 & spans] <- match.label
    temp[sum(keep) == 0 & ! spans] <- cigs[sum(keep) == 0 & ! spans]

    renamed <- as.character(temp)
    
    if (isTRUE(separate.snv)){
      renamed <- .self$.splitNonIndel(ref, renamed, rc, match_label = match.label,
                                    mismatch_label = mismatch.label,
                                    cut_site = target.loc, upstream = upstream,
                                    downstream = downstream)
      
    }

    .self$field("cigar_labels", renamed)
    
    #__________
    # Check if partial alns starting from different places have same labels
    strts <- start(.self$alns)[sum(keep) == 0 & ! spans]
    temp <- split(strts, renamed[sum(keep) == 0 & ! spans])
    if (! all(lengths(lapply(temp, unique)) == 1) ){
       rn <- unlist(renamed[sum(keep) == 0 & ! spans]) 
       gen_strt <- genome_to_target[as.character(unlist(temp))]
        
    }
    
    renamed
  },

  .splitNonIndel = function(ref, cig_labels, rc, match_label = "no variant",
                          mismatch_label = "SNV", cut_site = 17,
                          upstream = 8, downstream = 5){

  # NOTE: SNVS not called in partial alignments
  # Only consider mismatches up to (upstream) to the left of the cut and
  # (downstream) to the right of the cut
  # The cut site is between cut_site and cut_site + 1
  upstream = min(cut_site, upstream)
  downstream = min(downstream, nchar(ref) - cut_site)
  no_var <- which(cig_labels == match_label & mcols(.self$alns)$seq != ref)
  if (length(no_var) == 0) return(cig_labels)

  if ((cut_site-upstream + 1) < 0 | (cut_site + downstream) > length(ref)){
    stop("Specified range for detecting SNVs is greater than target range")
  }

  snv_range <- c((cut_site-upstream + 1):(cut_site + downstream))
  sqs <- mcols(.self$alns)$seq[no_var]
  if (isTRUE(rc)) sqs <- reverseComplement(sqs)

  no_var_seqs <- as.matrix(sqs)
  no_var_seqs <- no_var_seqs[,snv_range, drop = FALSE]

  rr <- strsplit(as.character(ref[snv_range]), "")[[1]]
  result <- apply(no_var_seqs, 1, function(x){
    snvs <- which((x != rr & x != "N")) - upstream - 1
    snvs[snvs >= 0] <- snvs[snvs >= 0] + 1
    sprintf("%s:%s", mismatch_label, paste(snvs, collapse = ","))
  })
  result[result == sprintf("%s:", mismatch_label)] <- match_label
  cig_labels[no_var] <- result
  cig_labels
  }

)
