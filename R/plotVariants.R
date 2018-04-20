# plotVariants generic -----
#'@title Plot alignments, frequencies and location of target sequence
#'@description Combines a plot of transcript structure, alleles
#' aligned with respect to a reference genome and a heatmap of
#' counts or proportions of each allele in a set of data.
#'@rdname plotVariants
#'@param obj The object to be plotted
#'@return A ggplot2 plot of the variants
#'@export
setGeneric("plotVariants", function(obj, ...) {
  standardGeneric("plotVariants")})
# -----

# plotVariants CrisprSet -----
#'@rdname plotVariants
#'@param txdb GenomicFeatures:TxDb object (default: NULL)
#'@param plotAlignments.args Extra arguments for plotAlignments
#'@param plotFreqHeatmap.args Extra arguments for plotFreqHeatmap
#'@param add.chr If target chromosome does not start with "chr", e.g.
#'"chr5", add the "chr" prefix.  (Default:TRUE)
#'@param ... extra arguments for plot layout
#'@seealso \code{\link{arrangePlots}} for general layout options
#' and \code{\link{annotateGenePlot}} for options relating
#' to the transcript plot.
#'@examples
#'#Load a CrisprSet object for plotting
#'data("gol_clutch1")
#'
#'#Load the transcript db.  This is a subset of the Ensembl Danio Rerio v73 gtf
#'# for the region 18:4640000-4650000 which includes the targeted gol gene
#'
#'library(GenomicFeatures)
#'fn <- system.file("extdata", "Danio_rerio.Zv9.73.gol.sqlite",
#'                  package = "CrispRVariants")
#'txdb <- loadDb(fn)
#'
#'# Plot the variants
#'p <- plotVariants(gol, txdb = txdb)
#'
#'#In the above plot, the bottom margin is too large, the legend is
#'#cut off, and the text within the plots should be larger.
#'#These issues can be fixed with some adjustments:
#'p <- plotVariants(gol, txdb = txdb,
#'                  plotAlignments.args = list(plot.text.size = 4, legend.cols = 2),
#'                  plotFreqHeatmap.args = list(plot.text.size = 4),
#'                  left.plot.margin = grid::unit(c(0.1,0,0.5,1), "lines"))
#'
setMethod("plotVariants", signature("CrisprSet"),
          function(obj, ..., txdb = NULL, add.chr = TRUE,
                   plotAlignments.args = list(),
                   plotFreqHeatmap.args = list()){

    include_txs <- TRUE
    if(!(class(txdb) == "TxDb" | class(txdb) == "TranscriptDb") ){
      if (is.null(txdb)){
        include_txs <- FALSE
      } else{
        stop("txdb should be a (GenomicFeatures) transcript database object")
      }
    }
    if (isFALSE(plotAlignments.args$create.plot) |
        isFALSE(plotFreqHeatmap.args$create.plot) ){
      stop("Argument create.plot can not be FALSE for plotVariants.\n",
        "Call plotAlignments or plotFreqHeatmap directly to get plot data.")
    }
    
    dots <- list(...)
    annotate_nms <- c("target.colour", "target.size",
                      "gene.text.size", "panel.spacing")

    annotate_args <- dots[names(dots) %in% annotate_nms]
    dots[annotate_nms] <- NULL

    arrange_nms <-  c("fig.height","col.wdth.ratio",
                      "row.ht.ratio", "left.plot.margin")
    arrange_args <- dots[names(dots) %in% arrange_nms]

    target <- obj$target
    if (isTRUE(add.chr) & isTRUE(include_txs)){
      # If adding "chr" to target chromosomes matches txdb chromosomes, do so
      target_levels <- GenomeInfoDb::seqlevels(target)
      txdb_levels <- GenomeInfoDb::seqlevels(txdb)
      wchr <- paste0("chr", target_levels)
      idxs <- wchr %in% txdb_levels
      target_levels[idxs] <- wchr[idxs]
      target <- GenomeInfoDb::renameSeqlevels(target, target_levels)
    }

    if (isTRUE(include_txs)){
      annotate_args <- modifyList(list(txdb = txdb, target = target),
                                  annotate_args)
      gene_p <- do.call(annotateGenePlot, annotate_args)
      gene_p <- ggplot2::ggplotGrob(gene_p)
    } else {
      arrange_args[["row.ht.ratio"]] <- c(0,1)
      gene_p <- grid::grid.rect(gp=grid::gpar(col="white"), draw = FALSE)
      no_ins <- nrow(obj$insertion_sites) == 0
      if (no_ins & ! "left.plot.margin" %in% names(arrange_args)){
        # Names tend to clip if there are no insertions, increase default
        arrange_args[["left.plot.margin"]] <- grid::unit(c(0.2,0,3,0.5), "lines")
      }
    }

    plotAlignments.args[["obj"]] <-  obj
    aln_p <- do.call(plotAlignments, plotAlignments.args)
    aln_p <- aln_p + theme(legend.spacing = grid::unit(0.2,"cm"))
  
    plotFreqHeatmap.args[["obj"]] <- obj
    heat_p <- do.call(plotFreqHeatmap, plotFreqHeatmap.args)
    heat_p <- heat_p + theme(plot.background=element_rect(fill = "transparent",
                                                        colour = NA),
                           plot.margin = grid::unit(c(1, 0.25, 0.5, 0), "lines"))

    arrange_args = modifyList(list("top.plot" = gene_p, "left.plot" = aln_p,
                                   "right.plot" = heat_p), arrange_args)
    result <- do.call(arrangePlots, arrange_args)

    return(result)
}) # -----

# arrangePlots ------ 
#'@title Arrange plots for plotVariants:CrisprSet
#'@description Arranges 3 plots in two rows.  The vertical margins of the
#'left.plot and right.plot constrained to be equal
#'@param top.plot ggplot grob, placed on top of the figure, spanning the figure
#'width
#'@param left.plot ggplot, placed in the second row on the left
#'@param right.plot ggplot, placed in the second row on the right.
#'y-axis labels are removed.
#'@param fig.height Actual height for the figure. If not provided,
#'figure height is the sum of the row.ht.ratio (Default: NULL)
#'@param row.ht.ratio Vector specifying row height ratio (Default: c(1,6))
#'@param col.wdth.ratio  Vector specifying column width ratio (Default: c(2, 1))
#'@param left.plot.margin Unit object specifying margins of left.plot.
#'Margins of right.plot are constrained by the left.plot.
#'@return The arranged plots
arrangePlots <- function(top.plot, left.plot, right.plot, fig.height = NULL,
                      col.wdth.ratio  = c(2, 1), row.ht.ratio = c(1,6),
                      left.plot.margin = grid::unit(c(0.1,0,3,0.2), "lines")){

  # Set the size ratio of the top and bottom rows
  plot_hts <- if (is.null(fig.height)){ row.ht.ratio
  }else { fig.height/sum(row.ht.ratio)*row.ht.ratio }

  # Remove y-axis labels from right plot
  right.plot <- right.plot + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank())

  # Adjust margins of left.plot
  left.plot <- left.plot + theme(plot.margin = left.plot.margin)

  # Convert plots to grobs, lock plot heights
  p2 <- ggplot2::ggplotGrob(left.plot)
  p3 <- ggplot2::ggplotGrob(right.plot)
  p3$heights <- p2$heights

  # Return arranged plots
  return(gridExtra::grid.arrange(top.plot,
         gridExtra::arrangeGrob(p2, p3, ncol = 2, widths = col.wdth.ratio),
         nrow = 2, heights = plot_hts, newpage = FALSE))
} # -----

# getOverlappingGenes -----
.getOverlappingGenes <- function(txdb, target, all.transcripts = TRUE){
    genomicfeatures <- requireNamespace("GenomicFeatures")
    stopifnot(isTRUE(genomicfeatures))

    columns <- c("GENEID","TXNAME", "EXONSTART", "EXONEND", "TXSTRAND",
                 "CDSCHROM", "CDSSTART", "CDSEND")

    trns <- GenomicFeatures::transcripts(txdb)
    exs <- findOverlaps(target, trns, ignore.strand = TRUE)

    if (length(exs) == 0){
        return(grid::grid.rect(gp=grid::gpar(col="white"), draw = FALSE))
    }
    # Find the genes that overlap
    genes <- AnnotationDbi::select(txdb,
              keys = as.character(mcols(trns[subjectHits(exs)])$tx_id),
              keytype = "TXID", columns = columns)
              

    # Find all (possibly non-overlapping) transcripts of overlapping genes
    if (isTRUE(all.transcripts)){
        genes <- suppressWarnings(AnnotationDbi::select(txdb,
              keys = unique(genes$GENEID),
              keytype = "GENEID",
              columns = c("TXID", columns) ))
    }
    
    genes
} # -----

# makeGeneSegments -----
.makeGeneSegments <- function(genes, txdb, target){
    gene_gr <- GRanges(seqnames(target)[1],
                     IRanges(genes$EXONSTART, genes$EXONEND),
                     strand = genes$TXSTRAND, txid = genes$TXID)
 
    # Assumption: if no end then also no start
    #no_cds <- is.na(genes$CDSSTART)
    #strts <- genes$CDSSTART
    #ends <- genes$CDSEND
    #strts[no_cds] <- genes$EXONSTART[no_cds]
    #ends[no_cds] <- genes$EXONEND[no_cds]
 
    has_cds <- ! is.na(genes$CDSSTART)
    cds_gr <- GRanges(seqnames(target)[1], 
                      IRanges(genes$CDSSTART[has_cds], genes$CDSEND[has_cds]),
                      strand = genes$TXSTRAND[has_cds], txid = genes$TXID[has_cds])
    #cds_gr <- GRanges(seqnames(target)[1], 
    #                 IRanges(strts, ends),
    #                 strand = genes$TXSTRAND, txid = genes$TXID)
 
    txid <- as.character(unique(genes$TXID))
    all_rngs <- c(cds_gr, gene_gr)
  

    # Get UTRs matching each transcript
    all_sections <- lapply(txid, function(tid) {
        exons <- disjoin(all_rngs[all_rngs$txid == tid]) 
        exons$type <- "exon"
        utr <- ! exons %in% cds_gr
        if (any(utr)) exons[utr]$type <- "utr"
        exons
    })

    gene_spans <- do.call(c,unname(lapply(all_sections, range)))
    all_exs <- do.call(c, all_sections)

    all_exs <- data.frame(start = start(all_exs),
                end = end(all_exs),
                ts = rep(1:length(all_sections), lapply(all_sections, length)),
                type = all_exs$type)
  
    colnames(all_exs) <- c("start", "end", "ts", "type")
  
    list("all_exs" = all_exs, "gene_spans" = gene_spans)
} # -----

# annotateGenePlot -----
#'@title Plots and annotates transcripts
#'@description Plots the gene structure, annotates this with the
#'target location
#'@param txdb A GenomicFeatures:TxDb object
#'@param target Location of target (GRanges)
#'@param target.colour Colour of box indicating targt region
#'@param target.size Thickness of box indicating target region
#'@param gene.text.size Size for figure label
#'@param panel.spacing Unit object, margin size
#'@param plot.title A title for the plot.  If no plot.title is supplied,
#'the title is the list of gene ids shown (default).
#'If plot.title == FALSE, the plot will not have a title.
#'@param all.transcripts If  TRUE (default), all transcripts of genes overlapping
#'the target are shown, including transcripts that do not themselves overlap the target.
#'If FALSE, only the transcripts that overlap the target are shown.
#'@return A ggplot2 plot of the transcript structures
annotateGenePlot <- function(txdb, target, target.colour = "red",
                        target.size = 1, gene.text.size = 10,
                        panel.spacing = grid::unit(c(0.1,0.1,0.1,0.1), "lines"),
                        plot.title = NULL, all.transcripts = TRUE){

    genomicfeatures <- requireNamespace("GenomicFeatures")
    stopifnot(isTRUE(genomicfeatures))

    genes <- .getOverlappingGenes(txdb, target, all.transcripts = all.transcripts)
    if ("grob" %in% class(genes)) return(genes)

    result <- .makeGeneSegments(genes, txdb, target)
    all_exs <- result$all_exs
    gene_spans <- result$gene_spans

    min_st <- min(all_exs$start)
    max_end <- max(all_exs$end)

    tcks <- unname(quantile(min_st:max_end, seq(1,100, by = 2)*0.01))
    tcks <- lapply(as(gene_spans, "GRangesList"), function(sp){
      tcks[tcks > start(sp) & tcks < end(sp)]
    })

    tck_lns <- lapply(tcks, length)
    tcks <- data.frame("tloc" = unlist(tcks),
                       "ys" = rep(1:length(tcks), lapply(tcks, length)))
    lns <- data.frame(tloc = c(start(gene_spans),end(gene_spans)),
                      ys = rep(seq_along(gene_spans),2))

    all_exs$ymax <- all_exs$ts + 0.3
    all_exs$ymin <- all_exs$ts - 0.3
    is_utr <- all_exs$type == "utr"
    all_exs$ymax[is_utr] <- all_exs$ts[is_utr] + 0.2
    all_exs$ymin[is_utr] <- all_exs$ts[is_utr] - 0.2

    target_df <- data.frame(xmin = start(target), xmax = end(target),
                            ymin = 0, ymax = ceiling(max(all_exs$ymax)))

    if (is.null(plot.title)){ plot.title <- paste(unique(genes$GENEID), sep = ";")}

    # Choose either right or left pointing arrows for the transcript plots
    #strands <- unlist(lapply(all_sections, function(x) as.character(strand(x[1]))))

    strands <- rep(as.character(strand(gene_spans)), tck_lns)
    strands[strands == "-"] <- 60
    strands[strands == "+"] <- 62
    tcks$shp <- as.integer(strands)

    p <- ggplot2::ggplot(tcks) +
      geom_point(aes_(x = quote(tloc), y = quote(ys), group = quote(ys),
                       shape = quote(shp)), size = 2) +
      geom_line(data = lns, aes_(x = quote(tloc), y = quote(ys),
                                  group = quote(ys))) +
      scale_shape_identity()

    p <- p + geom_rect(data = all_exs, fill = "black", color = "black",
                       aes_(xmin = quote(start), xmax = quote(end),
                             ymin = quote(ymin), ymax = quote(ymax)))

    p <- p + geom_rect(data = target_df,
                       aes_(xmin = quote(xmin), xmax = quote(xmax),
                             ymin = quote(ymin), ymax = quote(ymax)),
                       colour = target.colour, fill = NA, size = target.size)

    if (! isFALSE(plot.title)){
      p <- p + ggtitle(plot.title)
    }
    p <- p + theme_minimal() +
         theme(axis.text.x = element_text(size = gene.text.size),
            axis.text.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "white", colour=  NA),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.spacing = panel.spacing,
            text = element_text(size = gene.text.size),
            axis.ticks.y = element_blank()) +
         ylab(NULL) + xlab(NULL)

    return(p)
} # -----
