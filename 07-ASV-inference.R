## ----knitr_options, echo=FALSE--------------------------------------------------------------------------------------
knitr::opts_chunk$set(results = 'hide',
                      eval = FALSE,
                      #knitr::opts_knit$set(root.dir = "/fs/project/PAS0471/workshops/2020-12_micro/$USER"),
                      class.source = "r_code",
                      class.output = "r_output",
                      class.warning = "r_warning",
                      class.message = "r_warning",
                      class.error = "r_error")


## ----libload--------------------------------------------------------------------------------------------------------
.libPaths(new = '/fs/project/PAS0471/.R/4.0/')


## -------------------------------------------------------------------------------------------------------------------
packages <- c("tidyverse", "gridExtra", "dada2",
              "phyloseq", "DECIPHER", "phangorn")
pacman::p_load(char = packages)

# If you wanted to install and load these packages yourself,
# just make sure you have the pacman package installed:
## install.packages('pacman')
# Then, the code above would work, as it will install pacakes as needed.


## -------------------------------------------------------------------------------------------------------------------
# Dir with input fastq files:
indir <- 'data/processed/fastq_trimmed/second_trim'

# Dirs for output:
qc_dir <- 'analysis/QC'
filter_dir <- 'data/processed/fastq_filtered'
outdir <- 'analysis/ASV_inference'

# Fasta file with training data:
# (Check for an up-to-date version at <https://benjjneb.github.io/dada2/training.html>)
tax_key <- 'data/ref/silva_nr99_v138_train_set.fa' 

# File with sample metadata:
metadata_file <- 'metadata/sample_meta.txt'


## -------------------------------------------------------------------------------------------------------------------
fastqs_raw_F <- sort(list.files(indir, pattern = '_R1_001.fastq.gz', full.names = TRUE))
fastqs_raw_R <- sort(list.files(indir, pattern = '_R2_001.fastq.gz', full.names = TRUE))

head(fastqs_raw_F)


## -------------------------------------------------------------------------------------------------------------------
# Load and prepare sample metadata:
metadata_df <- read.table(file = metadata_file, sep = "\t", header = TRUE)

colnames(metadata_df)[1] <- 'sample_ID'
rownames(metadata_df) <- metadata_df$sample_ID

head(metadata_df)


## -------------------------------------------------------------------------------------------------------------------
metadata_df$sample_ID

head(basename(fastqs_raw_F))    # basename() strips the dir name from the filename


## -------------------------------------------------------------------------------------------------------------------
# sub() arguments: sub(pattern, replacement, vector)
# replace with "" -> replace with nothing
sample_IDs <- sub("-V4-V5_.*", "", basename(fastqs_raw_F))


## -------------------------------------------------------------------------------------------------------------------
identical(sort(metadata_df$sample_ID), sample_IDs)

setdiff(sort(metadata_df$sample_ID), sample_IDs) # REMOVE MISSING SAMPLES?


## -------------------------------------------------------------------------------------------------------------------
plotQualityProfile(c(fastqs_raw_F[1], fastqs_raw_R[1]))


## -------------------------------------------------------------------------------------------------------------------
pdf(file.path(qc_dir, "error_profiles.pdf"))  # Open a pdf file
for (sample_index in 1:3) {                   # Loop through first three file pairs
  print(plotQualityProfile(                   # Print plots into pdf
    c(fastqs_raw_F[sample_index], fastqs_raw_R[sample_index])) # F and R together
    )
}
dev.off()                                    # Close the pdf file


## -------------------------------------------------------------------------------------------------------------------
fastqs_filt_F <- file.path(filter_dir, paste0(sample_IDs, '_F_filt.fastq'))
fastqs_filt_R <- file.path(filter_dir, paste0(sample_IDs, '_R_filt.fastq'))


## -------------------------------------------------------------------------------------------------------------------
filter_dir('Filtering and Trimming...')
Sys.time()
filter_results <-
  filterAndTrim(fastqs_raw_F, fastqs_filt_F,
                fastqs_raw_R, fastqs_filt_R,
                truncLen = c(250,210),
                trimLeft = 10,
                maxN = 0,
                maxEE = c(2,2),
                truncQ = 2,
                rm.phix = FALSE,
                compress = FALSE, multithread = TRUE, verbose = TRUE) 
filter_dir('...Done!')
Sys.time()

head(filter_results)


## -------------------------------------------------------------------------------------------------------------------
fastqs_derep_F <- derepFastq(fastqs_filt_F, verbose = FALSE)
fastqs_derep_R <- derepFastq(fastqs_filt_R, verbose = FALSE)

names(fastqs_derep_F) <- sample_IDs
names(fastqs_derep_R) <- sample_IDs


## -------------------------------------------------------------------------------------------------------------------
print('Learning errors...')
Sys.time()

errF <- learnErrors(fastqs_derep_F, multithread = TRUE, verbose = TRUE)
errR <- learnErrors(fastqs_derep_R, multithread = TRUE, verbose = TRUE)

print('...Done!')
Sys.time()


## -------------------------------------------------------------------------------------------------------------------
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


## -------------------------------------------------------------------------------------------------------------------
print('Inferring ASVs (running the dada algorithm)...')
Sys.time()

dadaFs <- dada(fastqs_derep_F, err = errF, pool = FALSE, multithread = TRUE)
dadaRs <- dada(fastqs_derep_R, err = errR, pool = FALSE, multithread = TRUE)

print('...Done.')
Sys.time()


## -------------------------------------------------------------------------------------------------------------------
dadaFs[[1]]


## -------------------------------------------------------------------------------------------------------------------
mergers <- mergePairs(dadaFs, fastqs_derep_F,
                      dadaRs, fastqs_derep_R,
                      verbose = TRUE)


## -------------------------------------------------------------------------------------------------------------------
saveRDS(fastqs_derep_F, file = file.path(outdir, "fastqs_derep_F.rds"))
saveRDS(fastqs_derep_R, file = file.path(outdir, "fastqs_derep_R.rds"))

rm(fastqs_derep_F, fastqs_derep_R)


## -------------------------------------------------------------------------------------------------------------------
seqtab_all <- makeSequenceTable(mergers)

# The dimensions of the object are the number of samples and number of ASVs:
dim(seqtab_all)


## -------------------------------------------------------------------------------------------------------------------
table(nchar(getSequences(seqtab_all)))

# If you need to remove sequences of a particular length (e.g. too long):
# seqtab2 <- seqtab[, nchar(colnames(seqtab_all)) %in% seq(250,256)]


## -------------------------------------------------------------------------------------------------------------------
seqtab <- removeBimeraDenovo(seqtab_all,
                             method = "consensus",
                             multithread = TRUE,
                             verbose = TRUE)
ncol(seqtab)

# Proportion of retained sequences:
sum(seqtab) / sum(seqtab_all)


## -------------------------------------------------------------------------------------------------------------------
saveRDS(seqtab, file = file.path(outdir, "seqtab_V4.rds"))


## -------------------------------------------------------------------------------------------------------------------
getN <- function(x) sum(getUniques(x))

nreads_summary <- cbind(filter_results,
                        sapply(dadaFs, getN),
                        sapply(dadaRs, getN),
                        sapply(mergers, getN),
                        rowSums(seqtab))

colnames(nreads_summary) <- c('input', 'filtered', 'denoisedF',
                              'denoisedR', 'merged', 'nonchim')
rownames(nreads_summary) <- sample_IDs
head(nreads_summary)

outfile_nreads <- file.path(outdir, 'nreads_summary.txt')
write.table(nreads_summary, file = outfile_nreads,
            sep = "\t", quote = FALSE, row.names = TRUE)


## -------------------------------------------------------------------------------------------------------------------
print('Assigning taxa to ASVs...')
Sys.time()

taxa <- assignTaxonomy(seqtab, tax_key, multithread = TRUE)
colnames(taxa) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')

print('...Done.')
Sys.time()


## -------------------------------------------------------------------------------------------------------------------
# Rename ASVs: 
asv_seqs <- colnames(seqtab)
asv_headers <- paste('>ASV', 1:ncol(seqtab), sep = '_')


## -------------------------------------------------------------------------------------------------------------------
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(outdir, 'ASVs.fa'))


## -------------------------------------------------------------------------------------------------------------------
# Create a count (OTU) table:
otu_df <- t(seqtab)
row.names(otu_df) <- sub('>', '', asv_headers)

# Create a taxon table:
tax_df <- taxa
row.names(tax_df) <- sub(">", "", asv_headers)


## -------------------------------------------------------------------------------------------------------------------
ps <- phyloseq(otu_table(otu_df, taxa_are_rows = TRUE),
               sample_data(metadata_df),
               tax_table(tax_df))

# Write tables and phyloseq object to file:
write.table(otu_df, file.path(outdir, 'ASVs_counts.tsv'),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(tax_df, file.path(outdir, 'ASVs_taxonomy.tsv'),
            sep = "\t", quote = FALSE, col.names = NA)
saveRDS(ps, file = file.path(outdir, 'ps_V4.rds'))


## -------------------------------------------------------------------------------------------------------------------
print('Done with ASV inference.')


## ---- eval = FALSE--------------------------------------------------------------------------------------------------
## #seqtab<- readRDS('seqtab_V4.rds')
## 
## seqs <- getSequences(seqtab)
## 
## # This propagates to the tip labels of the tree.
## # At this stage ASV labels are full ASV sequence
## names(seqs) <- seqs
## alignment <- AlignSeqs(DNAStringSet(seqs),
##                        anchor = NA,
##                        iterations = 5,
##                        refinements = 5)
## 
## print('Computing pairwise distances from ASVs...')
## Sys.time()
## phang.align <- phyDat(as(alignment, 'matrix'), type = 'DNA')
## dm <- dist.ml(phang.align)
## treeNJ <- NJ(dm) # Note, tip order is not sequence order
## fit = pml(treeNJ, data = phang.align)
## print('...Done.')
## Sys.time()
## 
## print('....')
## Sys.time()
## fitGTR <- update(fit, k=4, inv=0.2)
## print('...Done.')
## Sys.time()
## 
## print('Computing likelihood of tree...')
## Sys.time()
## fitGTR <- optim.pml(fitGTR, model='GTR', optInv=TRUE, optGamma=TRUE,
##                       rearrangement = 'stochastic', control = pml.control(trace = 0))
## print('...Done'.)
## Sys.time()


## -------------------------------------------------------------------------------------------------------------------
knitr::purl(input = 'markdown/07-ASV-inference.Rmd',
            output = 'scripts/02-dada2_V4.R')


## #!/bin/bash

## #SBATCH --nodes=1

## #SBATCH --ntasks-per-node=28

## #SBATCH --time=100:00:00

## #SBATCH --account=PAS0471

## 
## module load gnu/9.1.0

## module load mkl/2019.0.3

## module load R/4.02

## 
## Rscript scripts/02-dada2_V4.R


## sbatch scripts/02-dada2.sh

