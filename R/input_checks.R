#__________________________________________________________________
# Input checks used during initialisation
#__________________________________________________________________
.checkReadsToTarget <- function(target, reference, target.loc,
                                reverse.complement, orientation,
                                chimeras){
  
  if (! class(target.loc) %in% c("numeric", "integer", "NULL")){
    stop("target.loc should be a single integer")
  }
  
  # Check lengths match
  
  ln_msg <- "readsToTarget accepts a single %s, as a%s object.\n"
  ln_msg2 <- "Use readsToTargets for processing multiple targets at once"
  if (length(target) > 1){
    stop(sprintf(ln_msg, "target", " GRanges"), ln_msg2)
  }
  if (! class(reference) %in% c("DNAString", "NULL")){
    if (length(reference) != 1){
      stop(sprintf(ln_msg, "reference sequence", " DNAString"), ln_msg2)
    }
  }  
  
  if (length(target.loc) > 1){
    stop(sprintf(ln_msg, "target.loc", "n integer"))
  }
  
  if (! is.null(reference)){
    # readsToTarget signature GAlignments does not take a reference
    if (! width(target) == nchar(reference)){
      stop("The width of the target should equal the",
           "number of characters in the reference sequence")
    }  
  }
  
  # Warn if target.loc is larger than width of target
  if (! is.null(target.loc)){
    if (target.loc > width(target)){
      warning(sprintf("target.loc (%s) > target width (%s).", 
                      target.loc, width(target)),
              "Is this intentional?")
    }
  }
  
  # Check reverse.complement and target are compatible
  if (reverse.complement == FALSE & orientation == "target"){
    warning(paste0("Conflicting choice of 'reverse.complement'",
                   " and 'orientation'.\nOrienting to the target",
                   " strand.\nSpecify 'orientation' for additional",
                   " options."))
  }
  
  if (class(chimeras) == "character"){
    # Chimeras can either be a GAlignments object or a tag
    # indicating how they are to be treated
    if (chimeras == "merge"){
      message(paste0("Caution: mergeChimeras assumes a sorted bam file\n",
                     "and has only been tested with bwa mem alignments!\n"))
    }
  }
}

.checkFnamesExist <- function(reads){
  if (! all(file.exists(reads))){
    stop("reads should be a vector of file names. ",
         "No file found matching some entries in reads")
  }
}

.checkReadsToTargets <- function(targets, primer.ranges, references){
  if (length(primer.ranges) > 0 & length(primer.ranges) != length(targets)){
    stop("primer.ranges should be the amplified regions, one per target")
  }
  if (length(targets) != length(references)){
    stop("A reference sequence for each target must be provided")
  }
  TRUE
}

.checkForPaired <- function(bams){
  is_first <- function(bams){ bitwAnd(mcols(bams)$flag, 64) }
  if ("flag" %in% names(mcols(bams))){
    if (any(is_first(bams))){
      warning("collapse.pairs set to FALSE but flags indicate paired reads",
              immediate. = TRUE)
      return(FALSE)
    }
  }
  TRUE
} 