#'@title refFromAlns
#'@description Reconstruct the reference sequence from alignments reads using the CIGAR
# string and MD tags.  The reference sequence is the sequence the reads were aligned
# to.  If a location is supplied, at least one alignment must span the location.
#'@param alns  Alignments to use for inferring the reference sequence
#'@param location  The location to infer the reference for.
# containing a single range (Optional).
#'@param ... additional arguments
#'@return The reference sequences corresponding to the provided alignments
#'@author Helen Lindsay
#'@rdname refFromAlns
#'@example 
#'bam <- system.file("extdata", "bam/ab1_ptena_wildtype_looking_embryo_1_s.bam",
#'                   package = "CrispRVariants")
#'alns <- GenomicAlignments::readGAlignments(bam,
#'           param = Rsamtools::ScanBamParam(tag = "MD", what = "seq"))
#' To get the reference sequence from the a given location:
#'location <- GenomicRanges::GRanges("chr17", IRanges::IRanges(23648420,23648430))
#'refFromAlns(alns, location = location)
#'@export
setGeneric("refFromAlns", function(alns, location, ...) {
  standardGeneric("refFromAlns")})


#'@param keep.names Should read names be added to the result if present?
#'(Default: FALSE)
#'@return A DNAStringSet (signature = c("GAlignments", "ANY"))
#'@rdname refFromAlns
setMethod("refFromAlns", signature(alns = "GAlignments", location = "ANY"),
    function(alns, location, ..., keep.names = FALSE){
      if (! "MD" %in% names(mcols(alns))){
        stop("Alignments must contain a metadata column with the MD tag.")
      }

      if (! "seq" %in% names(mcols(alns))){
          stop("Alignments must include the read sequences")
      }
      
      # Using cigar, get alignments with respect to the reference
      onref <- GenomicAlignments::sequenceLayer(mcols(alns)$seq, cigar(alns))
      mds <- mcols(alns)$MD
      if (any(is.na(mds))){
          stop("Cannot infer reference from MD tag when MD tag is NA")
      }
      
      # Split in two steps to include the split character in the results
      all.ops <- gsub("([\\^]*[ACGT]+)", ";\\1;", mds)
      
      # Remove leading ";" if present
      all.ops <- ifelse(startsWith(all.ops, ";"),
                          substr(all.ops, 2, nchar(all.ops)),
                          all.ops)

      split.ops <- strsplit(all.ops, "[;]+")
      so <- gsub("\\^", "", unlist(split.ops))
      is.chr <- grepl("[ACTG]+", so)
      idxs <- so
      idxs[is.chr] <- nchar(so[is.chr])
      idxs <- as.numeric(idxs)


      # Find ranges of characters to substitute
      start_locs <- unlist(lapply(relist(idxs, split.ops), function(x){
          cumsum(c(1,x))[1:length(x)]
      }))
      
      rngs <- IRanges(start_locs[is.chr], width = idxs[is.chr])
      ops.per.aln <- unlist(lapply(relist(is.chr, split.ops), sum))
      rngs <- relist(rngs, IRanges::PartitioningByWidth(ops.per.aln))
      replace.vals <- relist(so[is.chr], IRanges::PartitioningByWidth(ops.per.aln))
      refs <- replaceAt(onref, rngs, replace.vals)

      if (isTRUE(keep.names) & ! is.null(names(alns))){
          names(refs) <- names(alns)
      }

      refs
})


#'@return A DNAString (signature = c("GAlignments", "GRanges"))
#'@rdname refFromAlns
setMethod("refFromAlns", signature(alns = "GAlignments", location = "GRanges"),
          function(alns, location, ...){

    if (length(location) != 1){
        stop("A single genomic location is required")
    }
    
    kalns <- alns[seqnames(alns) == as.character(seqnames(location)) & 
                  start(alns) <= start(location) & end(alns) >= end(location)]
    
    if (length(kalns) == 0){
        stop("No alignment spans this location")
    }
    
    # As the reference is the same for all reads, only one alignment is needed
    kalns <- kalns[1]
    refs <- refFromAlns(kalns, keep.names = FALSE)

    # Take the substring of the location of interest
    ir <-  IRanges::IRanges(start(location)-start(kalns) + 1,
                             end(location)-start(kalns)+ 1)
    refs <- DNAString(substr(refs, start(ir), end(ir)))
    if (as.character(strand(location)) == "-"){
       refs <- Biostrings::reverseComplement(refs)
    }
    refs
})

