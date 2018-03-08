## ---- eval = FALSE---------------------------------------------------------
#  crispr_set <- readsToTarget(reads, target = target, reference = reference,
#                              target.loc = target.loc)
#  plotVariants(crispr_set)
#  # or use plotVariants(crispr_set, txdb) to additionally show the target
#  # location with respect to the transcripts if a Transcript Database
#  # txdb is available

## ---- message=FALSE, warning=FALSE-----------------------------------------
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

## ---- message=FALSE, warning = FALSE---------------------------------------
length(unique(ab1_fnames))
length(unique(fq_fnames))

## ---- message = FALSE, warning=FALSE, eval=FALSE---------------------------
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

## ---- message=FALSE--------------------------------------------------------
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

## ---- message=FALSE--------------------------------------------------------
library(rtracklayer)
# Represent the guide as a GenomicRanges::GRanges object
gd_fname <- system.file(package="CrispRVariants", "extdata/bed/guide.bed")
gd <- rtracklayer::import(gd_fname)
gd

## ---- message=FALSE--------------------------------------------------------
gdl <- GenomicRanges::resize(gd, width(gd) + 10, fix = "center")

## ---- eval=FALSE-----------------------------------------------------------
#  system("samtools faidx GRCHz10.fa.gz")
#  
#  reference=system(sprintf("samtools faidx GRCHz10.fa.gz %s:%s-%s",
#                           seqnames(gdl)[1], start(gdl)[1], end(gdl)[1]),
#                   intern = TRUE)[[2]]
#  
#  # The guide is on the negative strand, so the reference needs to be reverse complemented
#  reference=Biostrings::reverseComplement(Biostrings::DNAString(reference))
#  save(reference, file = "ptena_GRCHz10_ref.rda")

## --------------------------------------------------------------------------
ref_fname <- system.file(package="CrispRVariants", "extdata/ptena_GRCHz10_ref.rda")
load(ref_fname)
reference

## ---- tidy = FALSE---------------------------------------------------------
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

## ---- message=FALSE--------------------------------------------------------
# Note that the zero point (target.loc parameter) is 22
crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 22)
crispr_set

# The counts table can be accessed with the "variantCounts" function
vc <- variantCounts(crispr_set)
print(class(vc))

## ---- eval = FALSE---------------------------------------------------------
#  # In R
#  library(GenomicFeatures)
#  gtf_fname <- "Danio_rerio.GRCz10.81_chr17.gtf"
#  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_fname, format = "gtf")
#  saveDb(txdb, file= "GRCz10_81_chr17_txdb.sqlite")

## ---- echo=FALSE, message=FALSE--------------------------------------------
library(GenomicFeatures)
txdb_fname <- system.file("extdata/GRCz10_81_ptena_txdb.sqlite", 
                          package="CrispRVariants")
txdb <- loadDb(txdb_fname)

## ---- message = FALSE------------------------------------------------------
# The gridExtra package is required to specify the legend.key.height 
# as a "unit" object.  It is not needed to call plotVariants() with defaults
library(gridExtra)

# Match the clutch id to the column names of the variants
group <- md$Group

## ----ptena-plot, fig.width = 8.5, fig.height = 7.5, message = FALSE, fig.cap = "(Top) schematic of gene structure showing guide location (left) consensus sequences for variants (right) variant counts in each embryo."----
p <- plotVariants(crispr_set, txdb = txdb, gene.text.size = 8, 
    row.ht.ratio = c(1,8), col.wdth.ratio = c(4,2),
    plotAlignments.args = list(line.weight = 0.5, ins.size = 2, 
                               legend.symbol.size = 4),
    plotFreqHeatmap.args = list(plot.text.size = 3, x.size = 8, group = group, 
                                legend.text.size = 8, 
                                legend.key.height = grid::unit(0.5, "lines"))) 

## --------------------------------------------------------------------------
# Calculate the mutation efficiency, excluding indels that occur in the "control" sample
# and further excluding the "control" sample from the efficiency calculation
eff <- mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control")
eff

# Suppose we just wanted to filter particular variants, not an entire sample.
# This can be done using the "filter.vars" argument
eff2 <- mutationEfficiency(crispr_set, filter.vars = "6:1D", exclude.cols = "control")

# The results are the same in this case as only one variant was filtered from the control
identical(eff,eff2)

## --------------------------------------------------------------------------
sqs <- consensusSeqs(crispr_set)
sqs

# The ptena guide is on the negative strand.
# Confirm that the reverse complement of the "no variant" allele
# matches the reference sequence:
Biostrings::reverseComplement(sqs[["no variant"]]) == reference

## --------------------------------------------------------------------------
ch <- getChimeras(crispr_set, sample = "ptena 4")

# Confirm that all chimeric alignments are part of the same read
length(unique(names(ch))) == 1

# Set up points to annotate on the plot
annotations <- c(resize(gd, 1, fix = "start"), resize(gd, 1, fix = "end"))
annotations$name <- c("ptena_start", "ptena_end")

plotChimeras(ch, annotations = annotations)


## --------------------------------------------------------------------------
mutationEfficiency(crispr_set, filter.cols = "control", exclude.cols = "control",
                   include.chimeras = FALSE)

## ---- fig.width = 8.5, fig.height = 7.5, message = FALSE, warning = FALSE----

crispr_set_rev <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                                names = md$Short.name, target.loc = 22, 
                                orientation = "opposite")
plotVariants(crispr_set_rev)

## ---- warning = FALSE------------------------------------------------------
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

## ---- message = FALSE------------------------------------------------------
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


## ---- fig.height = 5, warning = FALSE--------------------------------------
p <- plotVariants(crispr_set, txdb = txdb)

## ---- fig.height = 5, warning = FALSE--------------------------------------
p <- plotVariants(crispr_set, txdb = txdb, row.ht.ratio = c(1,3))

## ---- fig.height = 5, message = FALSE, warning = FALSE---------------------
p <- plotVariants(crispr_set, txdb = txdb, col.wdth.ratio = c(4,1))

## --------------------------------------------------------------------------
# Load gol data set
library("CrispRVariants")
data("gol_clutch1")

