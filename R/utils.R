# mergeCrisprSets -----
#'@title Merge two CrisprSets
#'@description Merge two CrisprSet objects sharing a reference and target location
#'@author Helen Lindsay
#'@rdname mergeCrisprSets
setGeneric("mergeCrisprSets", function(x,y, ...) {
  standardGeneric("mergeCrisprSets")})


#'@rdname mergeCrisprSets
#'@param x A CrisprSet object
#'@param y A second CrisprSet object
#'@param order A list of sample names, matching the names in x and y,
#'specifying the order of the samples in the new CrisprSet. (Not implemented yet)
#'@param x.samples A subset of column names or indices to keep from CrispRSet x
#'(Default: NULL, i.e. keep all)
#'@param y.samples A subset of column names or indices to keep from CrispRSet y
#'(Default: NULL, i.e. keep all)
#'@param names New names for the merged CrisprSet object (Default: NULL)
#'@param ... extra arguments
#'@return A merged CrisprSet object
#'@examples
#'# Load the metadata table
#'md_fname <- system.file("extdata", "gol_F1_metadata_small.txt", package = "CrispRVariants")
#'md <- read.table(md_fname, sep = "\t", stringsAsFactors = FALSE)
#'
#'# Get bam filenames and their full paths
#'bam_fnames <- sapply(md$bam.filename, function(fn){
#'  system.file("extdata", fn, package = "CrispRVariants")})
#'
#'reference <- Biostrings::DNAString("GGTCTCTCGCAGGATGTTGCTGG")
#'gd <- GenomicRanges::GRanges("18", IRanges::IRanges(4647377, 4647399),
#'       strand = "+")
#'
#'crispr_set1 <- readsToTarget(bam_fnames[c(1:4)], target = gd,
#'       reference = reference, names = md$experiment.name[1:4], target.loc = 17)
#'crispr_set2 <- readsToTarget(bam_fnames[c(5:8)], target = gd,
#'       reference = reference, names = md$experiment.name[5:8], target.loc = 17)
#'mergeCrisprSets(crispr_set1,crispr_set2)
#'@export
setMethod("mergeCrisprSets", signature(x = "CrisprSet", y = "CrisprSet"),
          function(x,y, ..., x.samples = NULL, y.samples = NULL, names = NULL, 
                   order = NULL){
            # To do: add order vector

            warn_msg <- "CrisprSets must have the same %s for merging"
            ref = x$ref
            if (ref != y$ref){
              stop(sprintf(warn_msg, "reference sequence"))
            }
            target = x$target
            if (target != y$target){
              stop(sprintf(warn_msg, "guide sequence"))
            }
            t.loc = x$pars$target.loc
            if (t.loc != y$pars$target.loc){
              stop(sprintf(warn_msg, "target location (zero point)"))
            }
            
            if (is.null(x.samples)) x.samples <- seq_along(x$crispr_runs)
            if (is.null(y.samples)) y.samples <- seq_along(y$crispr_runs)
            cruns <- c(x$crispr_runs[x.samples], y$crispr_runs[y.samples])
            new_names <- sapply(cruns, function(x) x$name)
            
            if (! is.null(names)){
              if (! length(cruns) == length(names)){
                stop("Length of 'names' must equal the number of samples")
              }
              new_names <- names
            }
            #if (! is.null(order)){
            #  if (! length(order) == length(cruns)){
            #    stop("Length of 'order' must equal the number of samples")
            #  }
            #}
            cset <- CrisprSet(cruns, reference = ref, target = target,
                      names = new_names, target.loc = t.loc)
            
            # Recreate variant counts table
            temp <- alns(cset)
            cig_by_run <- relist(mcols(unlist(temp))$allele, temp)
            cset$.countCigars(cig_by_run)
            cset$.getInsertions()
            
            cset
            
          }) # -----

# indelCounts -----
#'@title Count the number of reads containing an insertion or deletion
#'@description Counts the number of reads containing a deletion or insertion
#'(indel) of any size in a set of aligned reads.
#'For countDeletions and countInsertions Reads may be filtered according to
#'whether they contain more than one indel of the same or different types.
#'@author Helen Lindsay
#'@param alns The aligned reads
#'@param multi.del  If TRUE, returns the exact number of deletions,
#'i.e., if one read contains 2 deletions, it contributes 2 to the
#'total count (default: FALSE)
#'@param del.and.ins If TRUE, counts deletions regardless of whether
#'reads also contain insertions.  If FALSE, counts reads that contain
#'deletions but not insertions (default: FALSE)
#'@param del.ops Cigar operations counted as deletions.  Default: c("D")
#'@param ... extra arguments
#'@rdname indelCounts
#'@examples
#'bam_fname <- system.file("extdata", "gol_F1_clutch_2_embryo_4_s.bam",
#'                          package = "CrispRVariants")
#'bam <- GenomicAlignments::readGAlignments(bam_fname, use.names = TRUE)
#'countDeletions(bam)
#'countInsertions(bam)
#'countIndels(bam)
#'indelPercent(bam)
#'@export
setGeneric("countDeletions", function(alns, ...) {
  standardGeneric("countDeletions")})


#'@return countDeletions: The number of reads containing a deletion (integer)
#'@rdname indelCounts
setMethod("countDeletions", signature("GAlignments"),
          function(alns, ..., multi.del = FALSE, del.and.ins = FALSE,
                   del.ops=c("D")){

  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))

  if (isTRUE(multi.del)){
    has_del <- any(cigar_ops %in% del.ops)
  } else{
    has_single_del <- sum(cigar_ops %in% del.ops) == 1
  }

  if (del.and.ins){
    if (multi.del)  return(sum(has_del))
    return(sum(has_single_del))
  }

  has_ins <- any(cigar_ops == "I")

  if (multi.del){
    return(sum(has_del & ! has_ins))
  }

  sum(has_single_del & ! has_ins)
})



#'@rdname indelCounts
#'@export
setGeneric("countInsertions", function(alns, ...) {
  standardGeneric("countInsertions")})


#'@param multi.ins  If TRUE, returns the exact number of insertions,
#'i.e., if one read contains 2 insertions, it contributes 2 to the
#'total count (default: FALSE)
#'@param ins.and.del If TRUE, counts insertions regardless of whether
#'reads also contain deletions  If FALSE, counts reads that contain
#'insertions but not deletions (default: FALSE)
#'@return countInsertions: The number of reads containing an insertion (integer)
#'@rdname indelCounts
setMethod("countInsertions", signature("GAlignments"),
          function(alns, ..., ins.and.del = FALSE, multi.ins = FALSE, del.ops = c("D")){

  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))

  if (isTRUE(multi.ins)){
    has_ins <- any(cigar_ops == "I")
  } else{
    has_single_ins <- sum(cigar_ops == "I") == 1
  }

  if (ins.and.del){
    if (multi.ins)  return(sum(has_ins))
    return(sum(has_single_ins))
  }

  has_del <- any(cigar_ops %in% del.ops)

  if (multi.ins) return(sum(has_ins &! has_del))
  result <- sum(has_single_ins &! has_del)
  result
})


#'@rdname indelCounts
#'@export
setGeneric("countIndels", function(alns) {
  standardGeneric("countIndels")})

#'@return countIndels: The number of reads containing at least one insertion
#'@rdname indelCounts
setMethod("countIndels", signature("GAlignments"),
          function(alns){
  cigar_ops <- CharacterList(explodeCigarOps(cigar(alns)))
  return(sum(any(cigar_ops %in% c("I", "D", "N"))))
})

#'@rdname indelCounts
setGeneric("indelPercent", function(alns) {
  standardGeneric("indelPercent")})


#'@rdname indelCounts
#'@return indelPercent: The percentage of reads containing an insertion or
#'deletion (numeric)
#'@export
setMethod("indelPercent", signature("GAlignments"),
          function(alns){
  return((countIndels(alns) / length(alns))*100)
}) # -----

# dispatchDots -----
#'@title dispatchDots
#'@description Update default values for func with values from dot args
#'@author Helen Lindsay
#'@param func Function to call
#'@param call If TRUE, call the function with the argument list and
#'return this result (Default: FALSE)
#'@param ... dot args to pass to function
#'@return A list of arguments to pass to func, or if call is TRUE, the result
#'of calling func with these arguments.
#'@rdname dispatchDots
#'@examples
#'# Set up a function to dispatch dot arguments to:
#'f <- function(a=1, b=2, c=3){
#'  print(c(a,b,c))
#'}
#'# Set up a function for passing dots:
#'g <- function(...){
#'  CrispRVariants:::dispatchDots(f, ...)
#'}
#'
#'g(a = 5)
#'g(a = 5, call = TRUE)
#'# Unrelated arguments will not be passed on
#'g(a = 5, d = 6)
dispatchDots <- function(func, ..., call = FALSE){
    if (isTRUE(isS4(func))){
      warning("dispatchDots may not work as expected with S4 functions")
    }
    func_defaults <- formals(func)
    #dots <- as.list(...)
    #print(dots)
    
    result <- utils::modifyList(func_defaults, list(...))[names(func_defaults)]
    if (isTRUE(call)) return(do.call(func, result))
    result
} # ------


# .getAxisCoords -----
#'@title .getAxisCoords
#'@description Manually specify x-tick locations and labels,
#'as sometimes ggplot defaults are too dense. Used internally
#'by CrispRVariants for creating alignment plot with plotAlignments.
#'@param locations character(n) Actual x coordinates, or the desired
#'range of the x coordinates.  If labels are provided, all tick locations
#'must be in locations and have a matching label.
#'@param labels character(n) labels for the x axis ticks.  Should be
#'the same length as locations if provided.  Note that if not all
#'tick locations are included in locations, it must be possible to
#'extrapolate labels from locations  (Default: NULL)
#'@param loc.boundaries numeric(i) Locations that must be included.
#'(Default: NULL)
#'@param lab.boundaries numeric(j) Labels that must be included.
#'(Default: c(-1,1), for showing the cut sites).  Boundaries must
#'be in labels and have a matching tick location.
#'@param label.at numeric(1) Add ticks when label modulo label.at
#'is zero (Default = 5)
#'@param min.tick.sep numeric(1) Minimum distance between ticks, excluding
#'boundary ticks.  (Default: 1)
#'@return A list containing vectors named tick_locs and tick_labs 
#'@author Helen Lindsay
#'@rdname getAxisCoords
.getAxisCoords <- function(locations, labels = NULL, loc.boundaries = NULL,
                           lab.boundaries = c(-1,1), label.at = 5,
                           min.tick.sep = 1){
    
    # Potential improvement:
    # Parameter min.tick.sep for preventing tick labels becoming too close
    # Parameter loc.at similar to label.at?
    # label.at if locations are character?
  
    input_labels <- labels
    
    # Input checks ----
    # If provided, labels should be the same length as locations
    if (! is.null(labels)){
      if (! length(locations) == length(labels)){
          stop("Length of locations should equal length of labels")
      }
    } else{
      labels <- locations
    }
    
    if (! is.numeric(label.at) | !
        class(loc.boundaries)  %in% c("numeric", "NULL", "integer")){
      stop("label.at and loc.boundaries should be numeric")
    }
    if (length(label.at) != 1){
      stop("label.at should be a single, numeric value")
    }
    if (! label.at > 0){
      stop("label.at can not be zero or less than zero")
    } # -----
  
    boundaries <- c(locations[locations %in% loc.boundaries],
                    locations[labels %in% lab.boundaries])

    tick_locations <- locations
    close_to_bdry <- outer(locations, boundaries, function(x,y){abs(x-y)})
    close_to_bdry <- apply(close_to_bdry, 1, min) < min.tick.sep
       
    # Select tick locations according to interval, remove ticks
    # too close to boundaries
    tick_locations <- tick_locations[labels %% label.at == 0 & ! close_to_bdry]
    tick_locations <- unique(c(tick_locations, boundaries))
    tick_locations <- tick_locations[order(tick_locations)]
  
    tick_locs_to_labels <- match(tick_locations, locations)
    tick_labels <- labels[tick_locs_to_labels]
    
    return(list(tick_locs = tick_locations, tick_labs = tick_labels))
  
} # -----


# .intersperse -----
#'@title .intersperse
#'@description create a vector of elements in outer interspersed with
#'elements in inner.  Similar to python zip. No element checking.
#'@param outer vector that will be the first and last elements
#'@param inner vector that will join elements of outer
#'@author Helen Lindsay
#'@examples CrispRVariants:::.intersperse(c(1:10), c(1:9)*10)
#'@return A vector interspersing elements of outer and inner.  If
#'outer is c(a,b,c) and inner is c(d,e), returns c(a,d,b,e,c)
#'@rdname intersperse
.intersperse <- function(outer, inner){
    result <- vector(class(outer), length(c(outer,inner)))
    result[2 * seq_along(outer) - 1] <- outer
    result[2 * seq_along(inner)] <- inner
    result
} # -----


# .explodeCigarOpCombs -----
#'@title .explodeCigarOpCombs
#'@description Breaks cigar strings into individual operations
#'@param cigar character(m) A vector of cigar strings
#'@param ops character(n) Which operations should be kept?
#'@return The operations, as a CharacterList
#'@author Helen Lindsay
#'@rdname explodeCigarOpCombs
#'@return Exploded cigar operations with operation widths
.explodeCigarOpCombs <- function(cigar, ops = GenomicAlignments::CIGAR_OPS){
    wdths <- GenomicAlignments::explodeCigarOpLengths(cigar, ops = ops)
    keep.ops <- GenomicAlignments::explodeCigarOps(cigar, ops = ops)
    result <- IRanges::CharacterList(relist(paste0(unlist(wdths),
                                                   unlist(keep.ops)),
                                            wdths))
} # -----


# selectOps -----
#'@title selectOps
#'@description select CIGAR operations in a region of interest.
#'@author Helen Lindsay
#'@rdname selectOps
setGeneric("selectOps", function(cigar, ...) {
  standardGeneric("selectOps")})


# selectOps -----
#'@param cigar CIGAR strings
#'@param ops  CIGAR operations to consider (Default: all)
#'@param op.regions (GRanges) Return operations only in these regions
#'@param pos An offset for the cigar ranges
#'@param ... Extra arguments (Not currently used)
#'@return A GRanges list of opertion locations in reference space
#'with a metadata column for the operation width in query space. 
#'@rdname selectOps
setMethod("selectOps", signature("character"),
          function(cigar, ..., ops = GenomicAlignments::CIGAR_OPS,
                   op.regions = NULL, pos = 1L){
            
    op_rngs <- cigarRangesAlongReferenceSpace(cigar, with.ops = TRUE,
                                              ops = ops, pos = pos)
    wdths <- explodeCigarOpLengths(cigar, ops = ops)
    op_labs <- .explodeCigarOpCombs(cigar, ops)
    
    temp <- unlist(op_rngs)
    mcols(temp)$qwidth <- unlist(wdths)
    op_rngs <- relist(temp, op_rngs)
    
    # If counting regions are specified, return operations in these regions
    if (! is.null(op.regions)){
       # Check for consistency - should insertions on the right border be kept?
       keep_rngs <- relist(overlapsAny(unlist(op_rngs), op.regions), op_rngs)
       op_rngs <- op_rngs[keep_rngs]
       op_labs <- op_labs[keep_rngs]
    }
    list(op_rngs = op_rngs, op_labels = op_labs)
}) # -----


isFALSE <- function(x){
  is.logical(x) && ! is.na(x) & ! isTRUE(x)
}
