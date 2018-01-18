#'@title getInsertionsTable
#'@description Returns a table of insertion sequences present in
#'a GAlignments object.  This table is aggregated and used by
#'plotAlignments.
#'@author Helen Lindsay
#'@param alns A GAlignments object
#'@param pos (Integer(1)) The amount by which to shift genomic
#'coordinates upstream to get coordinates relative to a display region
#'@return A data frame of insertion sequences, genomic and
#'relative locations
getInsertionsTable <- function(alns, pos = 1L){
    ins_on_qry <- cigarRangesAlongQuerySpace(cigar(alns), ops = "I")
    if (all(lengths(ins_on_qry) == 0)){ return(NULL) }
    ins_on_genome <- cigarRangesAlongReferenceSpace(cigar(alns),
                                                    pos = start(alns),
                                                    ops = "I") 
    
    ins_seqs <- extractAt(mcols(alns)$seq, ins_on_qry)
    # Subtract 1 to get coordinates of the leftmost mapped base
    genome_starts <- start(ins_on_genome) - 1
    relative_starts <- start(unlist(shift(ins_on_genome, 1 - pos)))
    
    df <- data.frame(start = relative_starts,
                     seq = as.character(unlist(ins_seqs)),
                     genomic_start = unlist(genome_starts),
                     idx = rep(seq_along(ins_seqs), lengths(ins_seqs)),
                     stringsAsFactors = FALSE)
    df
}