# narrowAlignments -----
#'@title Narrow a set of aligned reads to a target region
#'@description Aligned reads are narrowed to the target region.  In
#'the case of reads with deletions spanning the boundaries of the target,
#'reads are narrowed to the start of the deletion,
#'@author Helen Lindsay
#'@rdname narrowAlignments
#'@export
setGeneric("narrowAlignments", function(alns, target, ...) {
  standardGeneric("narrowAlignments")})


#'@param alns A GAlignments object including a metadata column "seq"
#'containing the sequence
#'@param target A GRanges object
#'@param reverse.complement Should the aligned reads be reverse complemented?
#'@param verbose (Default: FALSE)
#'@param minoverlap Minimum overlapping region between alignments and target.
#'If not specified, alignments must span the entire target region.
#'(Default: NULL)
#'@param clipping.ops CIGAR operations corresponding to clipping
#'(Default: c("S","H"))
#'@param match.ops CIGAR operations corresponding to a match, i.e.
#'a non-indel position (Default: c("M","X","="))
#'@param ... additional arguments
#'@return The narrowed alignments (GAlignments)
#'@rdname narrowAlignments
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "CrispRVariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'target <- GenomicRanges::GRanges("18", IRanges::IRanges(4647377, 4647399),
#'           strand = "+")
#'narrowAlignments(bam, target, reverse.complement = FALSE)
setMethod("narrowAlignments", signature("GAlignments", "GRanges"),
    function(alns, target, ..., reverse.complement,
             minoverlap = NULL, verbose = FALSE,
             clipping.ops = c("S","H"), match.ops = c("M","X","=")){
       
    # Narrowing example:
    # 3-4-5-6-7-8-9-10 Read
    #     5-6-7-8      Target sequence
    # target_start 5 - (read_start 3 - 1) = index 3
    # target_end 8 - target_start 5 + cigstart 3 = index
    
    # Notes:
    # Not using cigarNarrow as the seqs aren't narrowed
    # and want to keep indel operations bordering target range
    
    if (is.null(minoverlap)){
      # alns must span target
      alns <- alns[start(alns) <= start(target) & end(alns) >= end(target) &
                     seqnames(alns) == as.character(seqnames(target))]
    } else {
      alns <- alns[queryHits(findOverlaps(alns, target,
                                          minoverlap = minoverlap,
                                          ignore.strand = TRUE))]
    }
    if (length(alns) == 0){ return(GenomicAlignments::GAlignments()) }
    
    m_cols <- as.list(mcols(alns))
    if ("qual" %in% names(m_cols) & ! "seq" %in% names(m_cols)){
      stop("Metadata col 'seq' must be present if 'qual' is present")
    }
    
    if (isTRUE(verbose)) message("narrowing alignments\n")
    
    ref_ranges <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar(alns))
    genomic <- GenomicRanges::shift(ref_ranges, start(alns)-1)
    
    # Find the on target operations
    clipped <- unlist(GenomicAlignments::explodeCigarOps(cigar(alns)))
    clipped <- clipped %in% clipping.ops
    on_target <- unlist(start(genomic) <= end(target) & end(genomic) >= start(target))
    on_target <- on_target & ! clipped
    on_target <- relist(on_target, ref_ranges)
    
    # Make masks to pull out first and last on-target ranges
    first_on_tg <- ! duplicated(on_target) & on_target
    last_on_tg <- ! duplicated(on_target, fromLast = TRUE) & on_target
    
    # Find cases where the first or last on target operation is a match ("M"),
    # operations overlap boundaries (partial alignments don't need trimming)
    ops <- unlist(explodeCigarOps(cigar(alns)))
    first_op_m <- ops[unlist(first_on_tg)] %in% match.ops &
      start(alns) < start(target)
    last_op_m <- ops[unlist(last_on_tg)] %in% match.ops &
      end(alns) > end(target)
    
    # Narrow border match operations to match target
    # need to know where the target start (genomic) is wrt
    # sequence (query) coordinates
    cig <- GenomicAlignments::cigar(alns)
    q_ranges <- GenomicAlignments::cigarRangesAlongQuerySpace(cig)
    sq_starts <- start(q_ranges)[first_on_tg]
    
    genomic_offset_fom <- start(target) - start(genomic)[first_on_tg][first_op_m]
    sq_starts[first_op_m] <- sq_starts[first_op_m] + genomic_offset_fom
    
    # Adjust the genomic starting location for the alignments that will be narrowed
    genomic_starts <- start(genomic)[first_on_tg]
    genomic_starts[first_op_m] <- start(target)
    
    sq_ends <- end(q_ranges)[last_on_tg]
    genomic_offset_lom <- end(target) - start(genomic)[last_on_tg][last_op_m]
    sq_ends[last_op_m] <- start(q_ranges)[last_on_tg][last_op_m] + genomic_offset_lom
    
    # Adjust reference ranges to match
    ref_ranges <- unlist(ref_ranges)
    is_first_m <- which(unlist(first_on_tg))[first_op_m]
    
    # Calculate new ends before changing starts
    is_last_m <- which(unlist(last_on_tg))[last_op_m]
    lom_ends <- as.integer(start(ref_ranges[is_last_m]) + genomic_offset_lom)
    
    # Now adjust starts and ends
    start(ref_ranges[is_first_m]) <- as.integer(start(ref_ranges[is_first_m]) +
                                                  genomic_offset_fom)
    end(ref_ranges[is_last_m]) <- lom_ends
    
    # Get the width of on target reference ranges, recreate cigars
    # The width is the reference range unless an insertion, in which case the query range
    ref_ranges <- relist(ref_ranges, genomic)
    ref_ranges <- ref_ranges[on_target]
    qwdths <- unlist(width(q_ranges)[on_target])
    wdths <- width(ref_ranges)
    wdthl <- unlist(wdths)
    opsl <- ops[unlist(on_target)]
    wdthl[opsl == "I"] <- qwdths[opsl == "I"]
    new_cigs <- sapply(relist(paste0(wdthl, opsl), wdths), paste0, collapse = "")
    
    ga_params <- list(seqnames = seqnames(alns), pos = as.integer(genomic_starts),
                      cigar = new_cigs, names = names(alns), strand = strand(alns),
                      seqlengths = GenomeInfoDb::seqlengths(alns))
    
    if ("seq" %in% names(m_cols)){
      m_cols$seq <- subseq(m_cols$seq, start = as.numeric(sq_starts),
                           end = as.numeric(sq_ends))
    }
    if ("qual" %in% names(m_cols)){
      m_cols$qual <- subseq(m_cols$qual, start = as.numeric(sq_starts),
                            end = as.numeric(sq_ends))
    }
    
    ga_params <- c(ga_params, m_cols)
    new_alns <- do.call(GAlignments, ga_params)
    
    new_alns
}) # -----
