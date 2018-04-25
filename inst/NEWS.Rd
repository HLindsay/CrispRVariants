\name{NEWS}
\title{News for Package \pkg{CrispRVariants}}


\section{Changes in version 1.7.12}{
  \itemize{
      \item Change CrisprSe$insertion_sites "idxs" column to read indices
      for reads with multiple insertions
  }
}


\section{Changes in version 1.7.10}{
  \itemize{
      \item \code{CrisprSet} method \code{setCigarLabels} now allows renaming labels
  }
}


\section{Changes in version 1.7.9}{
  \itemize{
      \item Updates mergeCrisprSets in accordance with new order
      of counting operations
  }
}


\section{Changes in version 1.7.7}{
  \itemize{
      \item Adds min and max for guide bounding box in plot
      \item Fixes bug in \code{plotFreqHeatmap} caused when
      using "group" with a single row count matrix 
  }
}

\section{Changes in version 1.7.5}{
   \itemize{
      \item Update tests after changes to allele counting
      \item Added \code{create.plot} argument for \code{plotFreqHeatmap}
      with signature \code{CrisprSet}.
      \item Adds an option "style" to \code{plotAlignments} for colouring
      only mismatch nucleotides
      \item Changes to \code{narrowAlignments} for PacBio cigar format
      \item Bug fix in \code{collapsePairs}.  Only occurred
      when running outside of \code{readsToTarget}.
      \item Adds "alleles" accessor for relating variant labels
      to the truncated cigar strings
      \item Return unmergeable alignments instead of raising an error
      \item Minor code changes to make it easier to run a
      non-standard counting pipeline
      \item Code from initialisers split into separate files for 
      easier readability
  }
}


\section{Changes in version 1.7.4}{
   \itemize{
      \item Added \code{create.plot} argument for \code{plotFreqHeatmap}
      with signature \code{CrisprSet}.
  }
}


\section{Changes in version 1.7.2}{
   \itemize{
      \item Allowing plotting arbitrarily many subsets of
      aligned regions
  }
}


\section{Changes in version 1.7.1}{
   \itemize{
      \item Reorganising \code{plotAlignments} code for allowing
      plotting subsets of the aligned regions.  
      \item Minor change to transcript plot \code{plotVariants} to
      make background white not transparent.
  }
}


\section{Changes in version 1.5.9}{
   \itemize{
      \item More comprehensive input checking in \code{readsToTarget},
      removed redundant checking from CrisprSet initializer.
      \item Fix bug where sequences falsely called no variant if
      target is on the negative strand but positive strand reference
      given.
      \item Changed default SNV calling to 6 bases downstream instead
      of 5 to cover PAM
      \item Added tests for mismatched reference and target
  }
}

\section{Changes in version 1.5.8}{
   \itemize{
      \item Adds an option to filter variants by name when counting or
      plotting
  }
}

\section{Changes in version 1.5.7}{
   \itemize{
      \item Fixes major bug preventing filtering in \code{plotFreqHeatmap}
      \item Fixes bug in \code{mergeChimeras} if no chimeras mergeable
  }
}

\section{Changes in version 1.5.6}{
   \itemize{
      \item Autogenerate bam index for 
      \code{readsToTarget} option chimeras = "ignore"
  }
}

\section{Changes in version 1.5.3}{
  \itemize{
    \item \code{plotAlignments} now accepts the same filtering arguments
    as \code{variantCounts}
  }
}


\section{Changes in version 1.5.1}{
  \itemize{
    \item New argument \code{alleles} in \code{plotAlignments} and
    \code{plotFreqHeatmap} for selecting which alleles to display or specifiying
    a plotting order
    \item Removed unnecessary fields from CrisprRun class
    \item New accessor function \code{alns} to get alignments from a
     \code{CrisprSet}
    \item Improvements to plotAlignments to avoid unnecessary symbols in legend 
    \item Fix to header of plotFreqHeatmap when using type = "proportions"
    and providing sample order
    \item \code{consensusSeqs} now returns cigars as metadata by default
    \item Started indenting with four spaces at the start
  }
}


\section{Changes in version 1.3.7}{
  \itemize{
    \item Updates to vignette
    \item Fix bug removing variants by name in \code{variantCounts}
    \item Fixed argument \code{legend.symbol.size} being ignored in
    plotAlignmenta,DNAString-method. 
  }
}

\section{Changes in version 1.3.6}{
  \itemize{
    \item Fix in new function mergeChimeras when no chimeras present
  }
}

\section{Changes in version 1.3.5}{
  \itemize{
    \item New option "minoverlap" in \code{readsToTarget} allows reads that
    do not span the target region to be considered
    \item \code{plotAlignments} now works with character as well as DNAString
    objects
    \item Merging of long gaps mapped as chimeras now possible    
  }
}

\section{Changes in version 1.3.4}{
  \itemize{
      \item New function \code{refFromAlns} infers the reference sequence from aligned reads
      \item Fixed bug causing an empty plot when plotting a single alignment
      with a large deletion
      \item Changed \code{annotateGenePlot} from panel.margin to panel.spacing in accordance
      with recent ggplot2 versions
      \item Added "create.plot" argument to \code{plotAlignments} for signature 
      \code{CrisprSet} to make plot customisation easier.
      \item Fixed bug in argument names when all alignments are chimeric
      \item CrisprRun name now defaults to the coordinates when no name is provided
  }
}


\section{Changes in version 1.3.3}{
  \itemize{
      \item Fixed bug causing incorrect x-axis position in plotAlignments when
      strand unspecified
  }
}


\section{Changes in version 1.1.6}{
  \itemize{
  	\item \code{plotAlignments} can now mark codon boundaries if codon frame is specified.
  	\item Added citation
  }
}

\section{Changes in version 1.1.5}{
  \itemize{
  	\item More flexible specification of strand with new \code{readsToTarget} parameter 'orientation'
    \item Fixed warning caused by implicit embedding of S4 objects
    \item Added tests for 'plotAlignments' and 'annotateGenePlot'
    \item Minor speedup and internal restructuring of 'annotateGenePlot'
    \item Added CRISPR biocView
    \item Changed NEWS to rd format
  }
}

\section{Changes in version 1.1.4}{
  \itemize{
    \item Fixed a bug that prevented SNV settings being used in some circumstances
  }
}

\section{Changes in version 1.1.2}{
  \itemize{
    \item new function \code{consensusSeqs} returns the consensus sequences
    of the variant alleles
  }
}
