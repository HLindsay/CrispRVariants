## ---- eval = FALSE-------------------------------------------------------
#  crispr_set <- readsToTarget(reads, target = target, reference = reference,
#                              target.loc = target.loc)
#  plotVariants(crispr_set)
#  # or use plotVariants(crispr_set, txdb) to additionally show the target
#  # location with respect to the transcripts if a Transcript Database
#  # txdb is available

## ---- message=FALSE, warning=FALSE---------------------------------------
library(CrispRVariants)
library(sangerseqR)

# List AB1 filenames, get sequence names,  make names for the fastq files
# Note that we only include one ab1 file with CrispRVariants because
# of space constraints.  All bam files are included

data_dir <- system.file(package="CrispRVariants", "extdata/ab1/ptena")
fq_dir <- tempdir()
ab1_fnames <- dir(data_dir, "ab1$", recursive=TRUE, full=TRUE)
sq_nms <- gsub(".ab1","",basename(ab1_fnames))

# Replace spaces and slashes in filename with underscores
fq_fnames  <- paste0(gsub("[\ |\\/]", "_", dirname(ab1_fnames)), ".fastq")

# abifToFastq to read AB1 files and write to FASTQ
dummy <- mapply( function(u,v,w) {
        abifToFastq(u,v,file.path(fq_dir,w))
}, sq_nms, ab1_fnames, fq_fnames)

## ---- message=FALSE, warning = FALSE-------------------------------------
length(unique(ab1_fnames))
length(unique(fq_fnames))

## ---- message = FALSE, warning=FALSE, eval=FALSE-------------------------
#  library("Rsamtools")
#  
#  # BWA indices were generated using bwa version 0.7.10
#  bwa_index <- "GRCHz10.fa.gz"
#  bam_dir <- system.file(package="CrispRVariants", "extdata/bam")
#  fq_fnames <- file.path(fq_dir,unique(fq_fnames))
#  bm_fnames <- gsub(".fastq$",".bam",basename(fq_fnames))
#  srt_bm_fnames <- file.path(bam_dir, gsub(".bam","_s",bm_fnames))
#  
#  # Map, sort and index the bam files, remove the unsorted bams
#  for(i in 1:length(fq_fnames)) {
#    cmd <- paste0("bwa mem ", bwa_index, " ", fq_fnames[i],
#                  " | samtools view -Sb - > ", bm_fnames[i])
#    message(cmd, "\n"); system(cmd)
#    indexBam(sortBam(bm_fnames[i],srt_bm_fnames[i]))
#    unlink(bm_fnames[i])
#  }

## ---- message=FALSE------------------------------------------------------
# The metadata and bam files for this experiment are included with CrispRVariants
library("gdata")
md_fname <- system.file(package="CrispRVariants", "extdata/metadata/metadata.xls")
md <- gdata::read.xls(md_fname, 1)
md

# Get the bam filenames from the metadata table
bam_dir <- system.file(package="CrispRVariants", "extdata/bam")
bam_fnames <- file.path(bam_dir, md$bamfile)

# check that all files exist
all( file.exists(bam_fnames) )

## ---- message=FALSE------------------------------------------------------
library(rtracklayer)
# Represent the guide as a GenomicRanges::GRanges object
gd_fname <- system.file(package="CrispRVariants", "extdata/bed/guide.bed")
gd <- rtracklayer::import(gd_fname)
gd

## ---- message=FALSE------------------------------------------------------
gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")

## ---- eval=FALSE---------------------------------------------------------
#  system("samtools faidx GRCHz10.fa.gz")
#  
#  reference=system(sprintf("samtools faidx GRCHz10.fa.gz %s:%s-%s",
#                           seqnames(gdl)[1], start(gdl)[1], end(gdl)[1]),
#                   intern = TRUE)[[2]]
#  
#  # The guide is on the negative strand, so the reference needs to be reverse complemented
#  reference=Biostrings::reverseComplement(Biostrings::DNAString(reference))
#  save(reference, file = "ptena_GRCHz10_ref.rda")

## ------------------------------------------------------------------------
ref_fname <- system.file(package="CrispRVariants", "extdata/ptena_GRCHz10_ref.rda")
load(ref_fname)
reference

## ---- tidy = FALSE-------------------------------------------------------
# First read the alignments into R.  The alignments must include
# the read sequences and the MD tag
alns <- GenomicAlignments::readGAlignments(bam_fnames[[1]], 
          param = Rsamtools::ScanBamParam(tag = "MD", what = c("seq", "flag")),
          use.names = TRUE)

# Then reconstruct the reference for the target region.
# If no target region is given, this function will reconstruct
# the complete reference sequence for all reads.
rfa <- refFromAlns(alns, gdl)

# The reconstructed reference sequence is identical to the sequence
# extracted from the reference above
print(rfa == reference)

## ---- message=FALSE------------------------------------------------------
# Note that the zero point (target.loc parameter) is 22
crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 22)
crispr_set

# The counts table can be accessed with the "variantCounts" function
vc <- variantCounts(crispr_set)
print(class(vc))

## ---- eval = FALSE-------------------------------------------------------
#  # In R
#  library(GenomicFeatures)
#  gtf_fname <- "Danio_rerio.GRCz10.81_chr17.gtf"
#  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gtf")
#  saveDb(txdb, file= "GRCz10_81_chr17_txdb.sqlite")

## ---- echo=FALSE, message=FALSE------------------------------------------
library(GenomicFeatures)
txdb_fname <- system.file("extdata/GRCz10_81_ptena_txdb.sqlite", 
                          package="CrispRVariants")
txdb <- loadDb(txdb_fname)

## ---- message = FALSE----------------------------------------------------
# The gridExtra package is required to specify the legend.key.height 
# as a "unit" object.  It is not needed to call plotVariants() with defaults
library(gridExtra)

# Match the clutch id to the column names of the variants
group <- md$Group

## ----ptena_plot, fig.width = 8.5, fig.height = 7.5, message = FALSE, fig.cap = "(Top) schematic of gene structure showing guide location (left) consensus sequences for variants (right) variant counts in each embryo."----
p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 8, 
    row.ht.ratio = c(1,8), col.wdth.ratio = c(4,2),
    plotAlignments.args = list(line.weight = 0.5, ins.size = 2, 
                               legend.symbol.size = 4),
    plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, group = group, 
                                legend.text.size = 8, 
                                legend.key.height = grid::unit(0.5, "lines"))) 

## ------------------------------------------------------------------------
# Calculate the mutation efficiency, excluding indels that occur in the "control" sample
# and further excluding the "control" sample from the efficiency calculation
eff <- mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control")
eff

# Suppose we just wanted to filter particular variants, not an entire sample.
# This can be done using the "filter.vars" argument
eff2 <- mutationEfficiency(crispr_set, filter.vars = "6:1D", exclude.cols = "control")

# The results are the same in this case as only one variant was filtered from the control
identical(eff,eff2)

## ------------------------------------------------------------------------
sqs <- consensusSeqs(crispr_set)
sqs

# The ptena guide is on the negative strand.
# Confirm that the reverse complement of the "no variant" allele
# matches the reference sequence:
Biostrings::reverseComplement(sqs[["no variant"]]) == reference

## ------------------------------------------------------------------------
ch <- getChimeras(crispr_set, sample = "ptena 4")

# Confirm that all chimeric alignments are part of the same read
length(unique(names(ch))) == 1

# Set up points to annotate on the plot
annotations <- c(resize(gd, 1, fix = "start"), resize(gd, 1, fix = "end"))
annotations$name <- c("ptena_start", "ptena_end")

plotChimeras(ch, annotations = annotations)


## ------------------------------------------------------------------------
mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control",
                   include.chimeras = FALSE)

## ---- fig.width = 8.5, fig.height = 7.5, message = FALSE, warning = FALSE----

crispr_set_rev <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                                names = md$Short.name, target.loc = 22, 
                                orientation = "opposite")
plotVariants(crispr_set_rev)

## ---- warning = FALSE----------------------------------------------------
# We create a longer region to use as the "target"
# and the corresponding reference sequence
gdl <- GenomicRanges::resize(gd, width(gd) + 20, fix = "center")
reference <- Biostrings::DNAString("TCATTGCCATGGGCTTTCCAGCCGAACGATTGGAAGGTGTTTA")

# At this stage, target should be the entire region to display and target.loc should
# be the zero point with respect to this region
crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 10,
                            verbose = FALSE)

# Multiple guides are added at the stage of plotting
# The boundaries of the guide regions must be specified with respect to the
# given target region
p <- plotVariants(crispr_set, 
       plotAlignments.args = list(pam.start = c(6,35),
                              target.loc = c(10, 32),
                              guide.loc = IRanges::IRanges(c(6, 25),c(20, 37))))
p

## ---- message = FALSE----------------------------------------------------
# Setup for ptena data set
library("CrispRVariants")
library("rtracklayer")
library("GenomicFeatures")
library("gdata")

# Load the guide location
gd_fname <- system.file(package="CrispRVariants", "extdata/bed/guide.bed")
gd <- rtracklayer::import(gd_fname)
gdl <- resize(gd, width(gd) + 10, fix = "center")

# The saved reference sequence corresponds to the guide 
# plus 5 bases on either side, i.e. gdl
ref_fname <- system.file(package="CrispRVariants", 
                         "extdata/ptena_GRCHz10_ref.rda")
load(ref_fname)

# Load the metadata table, which gives the sample names
md_fname <- system.file(package="CrispRVariants",
                        "extdata/metadata/metadata.xls")
md <- gdata::read.xls(md_fname, 1)

# Get the list of bam files
bam_dir <- system.file(package="CrispRVariants", "extdata/bam")
bam_fnames <- file.path(bam_dir, md$bamfile)

# Check that all files were found
all(file.exists(bam_fnames))

crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 22,
                            verbose = FALSE)

# Load the transcript database
txdb_fname <- system.file("extdata/GRCz10_81_ptena_txdb.sqlite", 
                          package="CrispRVariants")
txdb <- AnnotationDbi::loadDb(txdb_fname)


## ---- fig.height = 5, warning = FALSE------------------------------------
p <- plotVariants(crispr_set, txdb = txdb)

## ---- fig.height = 5, warning = FALSE------------------------------------
p <- plotVariants(crispr_set, txdb = txdb, row.ht.ratio = c(1,3))

## ---- fig.height = 5, message = FALSE, warning = FALSE-------------------
p <- plotVariants(crispr_set, txdb = txdb, col.wdth.ratio = c(4,1))

## ------------------------------------------------------------------------
# Load gol data set
library("CrispRVariants")
data("gol_clutch1")

## ---- fig.height = 2.5, message = FALSE, warning = FALSE-----------------
library(GenomicFeatures)
p <- plotVariants(gol, plotAlignments.args = list(top.n = 3),
             plotFreqHeatmap.args = list(top.n = 3),
             left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"))

## ---- fig.height = 2.5, message = FALSE, warning = FALSE-----------------
plotVariants(gol, plotAlignments.args = list(top.n = 3),
             plotFreqHeatmap.args = list(top.n = 3, order = c(1,5,3)),
             left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"))

## ---- fig.height = 2.5, warning = FALSE----------------------------------
plotAlignments(gol, top.n = 3, ins.size = 6)

## ---- fig.height = 2.5---------------------------------------------------
plotAlignments(gol, top.n = 3, legend.symbol.size = 6)

## ---- fig.height = 3, warning = FALSE------------------------------------
plotAlignments(gol, top.n = 5, max.insertion.size = 25)

## ---- fig.height = 3, warning = FALSE------------------------------------
# Here we set a fairly high value of 50% for min.insertion.freq
# As ambiguous nucleotides occur frequently in this data set,
# there are no alleles passing this cutoff.
plotAlignments(gol, top.n = 5, min.insertion.freq = 50)

## ---- fig.height = 3, warning = FALSE------------------------------------
plotAlignments(gol, top.n = 5, max.insertion.size = 25, min.insertion.freq = 50)

## ---- fig.height = 2.5, warning = FALSE----------------------------------
# No white space between rows
plotAlignments(gol, top.n = 3, tile.height = 1)

## ---- fig.height = 3, warning = FALSE------------------------------------
# More white space between rows
plotAlignments(gol, top.n = 3, tile.height = 0.3)

## ---- fig.height = 2.5, warning = FALSE----------------------------------
plotAlignments(gol, top.n = 3, highlight.guide = FALSE)

## ---- fig.height = 3, message = FALSE------------------------------------
library(IRanges)
guide <- IRanges::IRanges(15,28)
plotAlignments(gol, top.n = 3, guide.loc = guide)

## ---- fig.height = 2.5---------------------------------------------------
# Here we increase the size of the axis labels and make
# two columns for the legend
plotAlignments(gol, top.n = 5, axis.text.size = 12, 
               legend.text.size = 12, legend.cols = 2)


## ---- fig.height = 3, warning = FALSE------------------------------------
# Don't highlight the PAM sequence
plotAlignments(gol, top.n = 3, highlight.pam = FALSE)

## ---- fig.height = 3, warning = FALSE------------------------------------

# Highlight 3 bases upstream to 3 bases downstream of the target.loc
plotAlignments(gol, top.n = 3, pam.start = 19, pam.end = 25)

## ---- fig.height = 3, warning = FALSE------------------------------------

plotAlignments(gol, top.n = 3, guide.loc = IRanges(5,10),
               pam.start = 8, pam.end = 13)


## ---- fig.height = 3, warning = FALSE------------------------------------
plotAlignments(gol, top.n = 3, line.weight = 3)

## ---- fig.height = 3, warning = FALSE------------------------------------
plotAlignments(gol, top.n = 3, codon.frame = 1)

## ---- eval = FALSE-------------------------------------------------------
#  plot_data <- plotAlignments(gol, top.n = 3, create.plot = FALSE)
#  names(plot_data)
#  # This data can be modified as required, then replotted using:
#  do.call(plotAlignments, plot_data)

## ----hmap_default, fig.height = 3, fig.width = 4, fig.align='center', fig.cap = "plotFreqHeatmap with default options"----
# Save the plot to a variable then add a title using ggplot2 syntax.
# If the plot is not saved to a variable the unmodified plot is displayed.
p <- plotFreqHeatmap(gol, top.n = 3)
p + labs(title = "A. plotFreqHeatmap with default options")

## ---- fig.height = 2.5, fig.width = 5, fig.align='center', fig.cap = "plotFreqHeatmap showing allele proportions"----
    
plotFreqHeatmap(gol, top.n = 3, type = "proportions")

## ---- fig.height = 2.5, fig.width = 4, fig.align='center', fig.cap = "plotFreqHeatmap with X-axis labels coloured by experimental group and tiles coloured by count instead of proportion"----
ncolumns <- ncol(variantCounts(gol))
ncolumns
grp <- rep(c(1,2), each = ncolumns/2)
p <- plotFreqHeatmap(gol, top.n = 3, group = grp, as.percent = FALSE)
p + labs(title = "B. coloured X labels with tiles coloured by count")

## ---- fig.height = 2.5, fig.width = 5, fig.align='center', fig.cap = "plotFreqHeatmap with labels showing allele proportions, header showing counts per sample and modified legend position."----
grp_clrs <- c("red", "purple")
p <- plotFreqHeatmap(gol, top.n = 3, group = grp, group.colours = grp_clrs,
                type = "proportions", header = "counts",
                legend.position = "bottom")
p <- p + labs(title = "C. Modified plotFreqHeatmap")
p

## ---- fig.height = 2.5, fig.width = 4, fig.align='center'----------------
plotFreqHeatmap(gol, top.n = 3, 
                legend.key.height = ggplot2::unit(1.5, "lines"))

## ---- eval = FALSE-------------------------------------------------------
#  var_counts <- variantCounts(gol, top.n = 3)
#  # (additional modifications to var_counts can be added here)
#  plotFreqHeatmap(var_counts)

## ---- fig.height = 2.5---------------------------------------------------
barplotAlleleFreqs(crispr_set, txdb = txdb)

## ---- fig.height = 2.5, message = FALSE----------------------------------
barplotAlleleFreqs(crispr_set, txdb = txdb, palette = "bluered")

## ---- fig.height = 2.5, message = FALSE----------------------------------
barplotAlleleFreqs(crispr_set, txdb = txdb, include.table = FALSE)

## ---- fig.height = 2.5---------------------------------------------------
var_counts <- variantCounts(crispr_set)
barplotAlleleFreqs(var_counts)

## ---- fig.height = 2.5---------------------------------------------------
rainbowPal9 <- c("#781C81","#3F4EA1","#4683C1",
                 "#57A3AD","#6DB388","#B1BE4E",
                 "#DFA53A","#E7742F","#D92120")

barplotAlleleFreqs(var_counts, classify = FALSE, bar.colours = rainbowPal9)


## ---- fig.height = 2.5---------------------------------------------------
# Classify variants as insertion/deletion/mixed
byType <- crispr_set$classifyVariantsByType()
byType

# Classify variants by their location, without considering size
byLoc <- crispr_set$classifyVariantsByLoc(txdb=txdb)
byLoc
# Coding variants can then be classified by setting a size cutoff
byLoc <- crispr_set$classifyCodingBySize(byLoc, cutoff = 6)
byLoc

# Combine filtering and variant classification, using barplotAlleleFreqs.matrix
vc <- variantCounts(crispr_set)

# Select variants that occur in at least two samples
keep <- names(which(rowSums(vc > 0) > 1))
keep

# Use this classification and the selected variants
barplotAlleleFreqs(vc[keep,], category.labels = byLoc[keep])

## ---- fig.height = 2.5---------------------------------------------------
p <- plotAlignments(gol, top.n = 3)
p + theme(legend.margin = ggplot2::unit(0, "cm"))

## ---- fig.height = 1-----------------------------------------------------

# Get a reference sequence
library("CrispRVariants")
data("gol_clutch1")
ref <- gol$ref

#Then to make the plot:
plotAlignments(ref, alns = NULL, target.loc = 22, ins.sites = data.frame())


## ---- message = FALSE, warning = FALSE-----------------------------------
library(Biostrings)
library(CrispRVariants)
library(rtracklayer)

## ---- warning = FALSE----------------------------------------------------
# This is a small, manually generated data set with a variety of different mutations 
bam_fname <- system.file("extdata", "cntnap2b_test_data_s.bam", 
                         package = "CrispRVariants")
guide_fname <- system.file("extdata", "cntnap2b_test_data_guide.bed",
                           package = "CrispRVariants")
guide <- rtracklayer::import(guide_fname)
guide <- guide + 5
reference <- Biostrings::DNAString("TAGGCGAATGAAGTCGGGGTTGCCCAGGTTCTC")

cset <- readsToTarget(bam_fname, guide, reference = reference, verbose = FALSE, 
                      name = "Default")
cset2 <- readsToTarget(bam_fname, guide, reference = reference, verbose = FALSE, 
                       chimera.to.target = 100, name = "Including long dels")

default_var_counts <- variantCounts(cset)
print(default_var_counts)
print(c("Total number of reads: ", colSums(default_var_counts)))

# With chimera.to.target = 100, an additional read representing a large deletion is 
# reported in the "Other" category.
var_counts_inc_long_dels <- variantCounts(cset2)
print(var_counts_inc_long_dels)
print(c("Total number of reads: ", colSums(var_counts_inc_long_dels)))

# This alignment can be viewed using `plotChimeras` 
ch <- getChimeras(cset2, sample = 1)
plotChimeras(ch, annotations = cset2$target)

