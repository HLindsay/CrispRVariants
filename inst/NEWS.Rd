\name{NEWS}
\title{News for Package \pkg{CrispRVariants}}

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
