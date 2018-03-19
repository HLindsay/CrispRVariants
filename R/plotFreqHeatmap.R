# plotFreqHeatmap generic -----
#'@title Plot a table of counts with colours indicating frequency
#'@description Creates a heatmap from a matrix of counts or proportions,
#'where tiles are coloured by the proportion and labeled with the value.
#'@rdname plotFreqHeatmap
#'@export
setGeneric("plotFreqHeatmap", function(obj, ...) {
  standardGeneric("plotFreqHeatmap")})
# -----

# plotFreqHeatmap matrix -----
#'@rdname plotFreqHeatmap
#'@param obj A matrix of counts with rows = feature, columns = sample
#'@param col.sums Alternative column sums to be used for calculating
#'the tile colours if as.percent = TRUE, e.g. if "obj" is a subset of
#'a larger data set.  If "NULL" (default), the column sums of "obj" are used.
#'@param header Alternative column titles, e.g. column sums for the
#'unfiltered data set when obj is a subset.  If set to "NA", column sums
#'of obj are displayed. If "NULL", no header is displayed (Default: NA).
#'@param header.name Label for the header row (Default: "Total")
#'@param group Grouping factor for columns.  If supplied, columns are
#'ordered to match the levels  (Default: NULL)
#'@param group.colours Colours for column groups, should match levels of "group".
#'If "NULL", groups are coloured differently (Default: NULL)
#'@param as.percent  Should colours represent the percentage of reads
#'per sample (TRUE) or the actual counts (FALSE)?  (Default: TRUE)
#'@param x.axis.title A title for the x-axis.  (Default: NULL)
#'@param x.size Font size for x-labels (Default: 16)
#'@param y.size Font size for y-labels (Default: 16)
#'@param x.angle  Angle for x-labels (Default: 90, i.e. vertical)
#'@param legend.text.size Font size for legend (Default: 16)
#'@param plot.text.size Font size counts within plot (Default: 3)
#'@param line.width Line thickness of title box'
#'@param x.hjust Horizontal justification of x axis labels (Default: 1)
#'@param legend.position The position of the legend (Default: right)
#'@param x.labels X-axis labels (Default: NULL, column.names of the matrix,
#'doesn't do anything at the moment)
#'@param legend.key.height The height of the legend key, as a "unit" object.
#'(See \code{\link[grid]{unit}}).
#'@param ... additional arguments
#'@return The ggplot2 plot of the variant frequencies
setMethod("plotFreqHeatmap", signature("matrix"),
          function(obj, ..., col.sums = NULL, header = NA, header.name = "Total",
                   group = NULL, group.colours = NULL, as.percent = TRUE,
                   x.axis.title = NULL, x.size = 6, y.size = 8, x.angle = 90,
                   legend.text.size = 6, plot.text.size = 3, line.width = 1,
                   x.hjust = 1, legend.position = "right", x.labels = NULL,
                   legend.key.height = grid::unit(1, "lines")) {

  # Potential improvements:
  # Allow a separate object for colours
  # param colour.vals A matrix of the same dimensions as obj containing
  # numbers that will be used to colour the heatmap.
  
  # If there should be a header, create it and make space in plot
  add_header <- ifelse(is.null(header), FALSE, TRUE)
  
  if (isTRUE(add_header)){
    # If header is not specified, it should be colSums
    if (is.na(header[1])) header <- colSums(obj)
    
    # If header is specified, it should have enough values
    if (!length(header) == ncol(obj)){
      stop("Header length should equal to the number of columns of obj")
    }
    temp <- rownames(obj)
    obj <- rbind(Total = rep(NA, ncol(obj)), obj)
    rownames(obj) <- c(header.name, temp)
  }
  
  # If a sample group is supplied, reorder the columns of counts
  if (! is.null(group)){
    if (! length(group) == ncol(obj)){
      stop("Group should be a vector or factor the same length as ncol(obj)")
    }
    if (! class(group) == "factor") group <- factor(group, levels = unique(group))

    obj <- obj[,order(group), drop = FALSE]
    if (length(header) > 1){
      header <- header[order(group)]
    }
    group <- group[order(group)]

    # if group.colours are not provided, use defaults
    if (is.null(group.colours)){
      # default is a colourblind safe palette
      clrs <- c("#000000","#0072B2","#E69F00","#009E73",
                "#56B4E9","#D55E00","#CC79A7")

      clrs <- clrs[group]
    } else {
      clrs <- group.colours[group]
    }
  }

  counts <- reshape2::melt(obj)
  colnames(counts) <- c("Feature", "Sample","Count")
  counts$Feature <- factor(counts$Feature, levels = rev(levels(counts$Feature)))

  # Create coloured tile background
  if (isTRUE(as.percent)){
    
    # If col.sums are provided, use these, otherwise use the column sums
    if (is.null(col.sums)){
      col.sums <- colSums(obj, na.rm = TRUE)
    } else if (! length(col.sums) == ncol(obj)){
      stop("Number of col.sums should equal to the number of columns of obj")
    }
    totals <- as.numeric(col.sums)
    
    m <- t(t(obj)/totals) * 100
    m <- reshape2::melt(m)
    colnames(m) <- c("Feature", "Sample","Percentage")
    m$Feature <- factor(m$Feature, levels = rev(levels(m$Feature)))
    
    g <- ggplot(m) + geom_tile(aes_(x = quote(Sample),
                                     y = quote(Feature),
                                     fill = quote(Percentage)))
  }
  else{
    g <- ggplot(counts) + geom_tile(aes_(x = quote(Sample),
                                          y = quote(Feature),
                                          fill = quote(Count)))
  }

  # Add the count numbers to the boxes,
  # if totals are included, add bold boxes to indicate these
  counts$ff <- "plain"
  box_coords <- data.frame(xmin = integer(), xmax=integer(),
                           ymin=integer(), ymax=integer())
  box_row <- 1
  panel_ranges <- ggplot_build(g)$layout$panel_ranges[[1]]
  xranges <- panel_ranges$x.range
  yranges <- panel_ranges$y.range

  if (isTRUE(add_header)){
    idxs <- which(is.na(counts$Count))
    counts$Count[idxs] <- header
    counts$ff[idxs] <- "bold"
    box_coords[box_row,] <- c(min(xranges), max(xranges),
                              nrow(obj) -0.5, max(yranges))
  }

  if (nrow(box_coords) > 0){
    g <- g + geom_rect(data = box_coords,
                       aes_(xmin = quote(xmin), xmax = quote(xmax),
                           ymin = quote(ymin), ymax = quote(ymax)),
                    color = "black", size = line.width, fill = "transparent")
  }
  # Plot the counts
  g <- g + geom_text(data = counts,
                     aes_(x = quote(Sample), y = quote(Feature),
                           label = quote(Count), fontface = quote(ff)),
                     size = plot.text.size)

  # Colour the boxes - white for 0, darkred for highest
  hmcols <- grDevices::colorRampPalette(c("white","gold","orange",
                                          "orangered","red", "darkred"))(50)

  # No colour bar legend if the max count is 1
  if (max(counts$Count) > 1){
    if (isTRUE(as.percent)){
      g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white",
                                    guide = "legend", limits = c(0, 100))
    } else {
      g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white",
                                    guide = "legend")
    }
  } else {
    g <- g + scale_fill_gradientn(colours = hmcols, na.value = "white",
                                  guide = "none")
  }

  # Set plot labels
   if (is.null(x.labels)) x.labels <- colnames(obj)
  g <- g + ylab(NULL) + xlab(x.axis.title) + theme_bw() +
  theme(axis.text.x = element_text(size = x.size, angle = x.angle, 
                                     hjust = x.hjust, vjust = 1),
          axis.text.y = element_text(size = y.size),
          legend.text = element_text(size = legend.text.size),
          legend.title = element_text(size = legend.text.size),
          legend.position = legend.position,
          legend.key.height = legend.key.height)

  # Colour xlabels by group, make sure not to replace original args
  if (! is.null(group)){
    g <- g + theme(axis.text.x=element_text(colour= clrs,
                   size = x.size, angle = x.angle, hjust = x.hjust))
  }

  return(g)
}) # ------


# plotFreqHeatmap CrisprSet ----- 
#'@rdname plotFreqHeatmap
#'@param top.n  Show the n top ranked variants.  Note that if the nth and n+1th
#'variants have equal rank, they will not be shown.   (Default: 50)
#'@param min.freq i (%) only plot variants with frequency >= i% in at least
#' one sample (Default: 0, i.e no frequency cutoff)
#'@param min.count i (integer) only plot variants with count >= i in at least
#' one sample (default: 0, i.e no count cutoff)
#'@param type Plot either "counts" or "proportions"
#'@param order A list of column names or indices specifying the order of the
#'columns in the plot
#'@param alleles A list of alleles to include.  Can be used to display
#'only alleles of interest or to order the alleles.
#'The default value NULL means all alleles passing the frequency cut offs 
#'will be included.
#'@examples
#'#Load a CrisprSet object for plotting
#'data("gol_clutch1")
#'
#'# Plot the frequency heatmap
#'plotFreqHeatmap(gol)
setMethod("plotFreqHeatmap", signature("CrisprSet"),
          function(obj, ..., top.n = 50, min.freq = 0, min.count = 1,
                   type = c("counts", "proportions"), order = NULL,
                   alleles = NULL) {

  result <- obj$heatmapCigarFreqs(top.n = top.n, min.freq = min.freq,
                     min.count = min.count, type = type, order = order,
                     alleles = alleles, ...)
  return(result)
}) # ------