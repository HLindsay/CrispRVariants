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

