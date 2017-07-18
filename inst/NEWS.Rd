\name{NEWS}
\title{News for Package \pkg{CrispRVariants}}

\section{Changes in version 1.5.1}{
  \itemize{
    \item Fix to header of plotFreqHeatmap when using type = "propotions"
    and providing sample order
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
