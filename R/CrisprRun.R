# CrisprRun class definition -----
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
#'@param verbose Print information about initialisation progress (Default: FALSE)
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
) # -----

# CrisprRun initializer -----
CrisprRun$methods(
  initialize = function(bam, target, rc = FALSE, name = NULL,
                        chimeras = GenomicAlignments::GAlignments(),
                        verbose = FALSE){
    
    # If no name is provided, use the coordinates
    if (is.null(name)){
        name <<- sprintf("%s:%s-%s", seqnames(target), start(target), end(target))
    } else { name <<- name}

    if (isTRUE(verbose)){
      message(sprintf("\nInitialising CrisprRun %s\n", .self$name))
    }
    
    alns <<- bam
    chimeras <<- chimeras
    chimera_combs <<- .self$.splitChimeras()

    if (length(bam) == 0) { return() }

    # recalculate genome.ranges in case of keep_unpaired = FALSE, not tested
    .self$getInsertionSeqs(target = target)
  }, # -----

  # show -----
  show = function(){
    print(c(class(.self), sprintf("CrisprRun object named %s, with %s on target alignments.",
                                  .self$name, length(.self$alns)), .self$alns))
  }, # ------

  # getInsertionSeqs -----
  getInsertionSeqs = function(target){
' 
Description:
  Return a table relating insertion sequences to alignment indices
Input parameters:
'
    df <- getInsertionsTable(.self$alns, pos = start(target))
    if (is.null(df)) {return(NULL)}
    df$label <- .self$cigar_labels[df$idx]
    
    agg <- aggregate(seq_len(nrow(df)), by = df[,c(1:3,5)], c,
                     simplify=FALSE)
    agg$count <- lengths(agg$x)
    colnames(agg) <- c("start", "seq", "genomic_start",
                       "cigar_label", "idxs", "count")
    return(agg)
  }, # -----

  # .checkNonempty -----
  .checkNonempty = function(){
    if (length(.self$alns) == 0){
      message("No on target alignments")
      return(FALSE)
    }
    TRUE
  }, # -----

  # splitChimeras -----
  .splitChimeras = function(){
    splits <- split(cigar(.self$chimeras), names(.self$chimeras))
    if (length(splits) == 0) return(data.frame())
    combination <- sapply(splits, base::paste, collapse=";")
    tt <- as.data.frame(table(combination))
    tt <- tt[order(tt$Freq, decreasing = TRUE),]
    tt
  }, # -----

  # getCigarLabels -----
  getCigarLabels = function(target, target.loc, genome_to_target, ref,
                            separate.snv, rc, match.label, mismatch.label,
                            keep.ops = c("I","D","N"),
                            upstream = min(target.loc, 8),
                            downstream = min(6, width(ref) - target.loc + 1),
                            regions = NULL, snv.regions = NULL){
    '
Description:
  Sets the "cig_labels" field, returns the cigar labels.

Input parameters:
  target:           (GRanges) the counting region.
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
  downstream:       distance downstream of the cut site to call SNVs
  regions:          IRanges(k) Regions for counting insertions and
                    deletions.  Insertions on the right border are not
                    counted.
  snv.regions       Regions for calling SNVS'
    
    
    # target.loc is wrt reference sequence, -ve if ref is -ve
    # upstream and downstream defined wrt the targt.loc
    
    # Find indels
    
    labels <- indelLabels(.self$alns, rc,
                          genome.to.pos = genome_to_target)
    
    relative_start <- max(1, target.loc - upstream + 1)
    relative_end <- min(target.loc + downstream, nchar(ref))
    regions <- IRanges(relative_start, relative_end)
    
    # Find snvs
    
    if (any(labels == "") & isTRUE(separate.snv)){
      mm_labs <- mismatchLabels(.self$alns[labels == ""], target, ref, 
                       regions = IRanges(relative_start, relative_end), 
                       genome.to.pos = genome_to_target)
    
      labels[labels == ""] <- mm_labs  
    }
    
    # Remaining are either no variant or partial
    
    # Adjust cigar for partial alignments
    spans <- width(.self$alns) == nchar(ref) 
    labels[labels == "" & spans] <- match.label
    labels[labels == "" & !spans]  <- cigar(.self$alns)[labels == "" & ! spans]
    
    .self$field("cigar_labels", labels)
    mcols(.self$alns)$"allele" = labels
    labels

    ## Adjust cigar for partial alignments
    #spans <- width(.self$alns) == nchar(ref) 
    #temp[ops_per_aln == 0 & spans] <- match.label
    #temp[ops_per_aln == 0 & ! spans] <- cigs[ops_per_aln == 0 & ! spans]
#
     #__________
    # Check if partial alns starting from different places have same labels
    #strts <- start(.self$alns)[ops_per_aln == 0 & ! spans]
    #temp <- split(strts, renamed[ops_per_aln == 0 & ! spans])
    #if (! all(lengths(lapply(temp, unique)) == 1) ){
    #   rn <- unlist(renamed[ops_per_aln == 0 & ! spans]) 
    #   gen_strt <- genome_to_target[as.character(unlist(temp))]
    #    
    #}
    
    #renamed
  } # -----

  # .splitNonIndel -----
#  .splitNonIndel = function(ref, cig_labels, rc, match_label = "no variant",
#                          mismatch_label = "SNV", cut_site = 17,
#                          upstream = 8, downstream = 6){
#
#    # NOTE: SNVS not called in partial alignments
#    # Only consider mismatches up to (upstream) to the left of the cut and
#    # (downstream) to the right of the cut
#    # The cut site is between cut_site and cut_site + 1
#    upstream = min(cut_site, upstream)
#    downstream = min(downstream, nchar(ref) - cut_site)
#    
#    is_match <- cig_labels == match_label
#    
#    # If negative strand, rc reference for checking identity
#    # Sequences are later reverse complemented to match the ref
#    test_ref <- ref 
#    if(isTRUE(rc)){ test_ref <- Biostrings::reverseComplement(test_ref)}
#    no_var <- which(is_match & mcols(.self$alns)$seq != test_ref)
#    
#    if (length(no_var) == 0) return(cig_labels)
#  
#    if ((cut_site-upstream + 1) < 0 | (cut_site + downstream) > length(ref)){
#      stop("Specified range for detecting SNVs is greater than target range")
#    }
#  
#    snv_range <- c((cut_site-upstream + 1):(cut_site + downstream))
#    sqs <- mcols(.self$alns)$seq[no_var]
#    if (isTRUE(rc)) sqs <- reverseComplement(sqs)
#  
#    no_var_seqs <- as.matrix(sqs)
#    no_var_seqs <- no_var_seqs[,snv_range, drop = FALSE]
#  
#    rr <- strsplit(as.character(ref[snv_range]), "")[[1]]
#    result <- apply(no_var_seqs, 1, function(x){
#      snvs <- which((x != rr & x != "N")) - upstream - 1
#      snvs[snvs >= 0] <- snvs[snvs >= 0] + 1
#      sprintf("%s:%s", mismatch_label, paste(snvs, collapse = ","))
#    })
#    result[result == sprintf("%s:", mismatch_label)] <- match_label
#    cig_labels[no_var] <- result
#    cig_labels
#    } # -----

)
