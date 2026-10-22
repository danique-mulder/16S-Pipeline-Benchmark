#### Install and load packages ############################################################
# Activate bioconductor repos
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(c('BiocManager'), quiet = TRUE, update = FALSE, ask = FALSE)
}

require(BiocManager)

# Packages to install/load
pkg <- c("optparse", "dada2", "ggplot2", "parallel", "Biostrings", "dplyr") 

# install uninstalled packages
pkg_already_installed <- pkg[pkg %in% installed.packages()]
packages_to_install <- setdiff(pkg, pkg_already_installed)

if (length(packages_to_install)) {
  install(packages_to_install)
}

loaded <- FALSE
loaded <- sapply(pkg, require, character.only = TRUE)
print(as.matrix(loaded))

#### Parse command line arguments 
# Define the options
option_list <- list(
  make_option(c("-p", "--path"), type = "character", default = NULL, help = "Input file"),
  make_option(c("-d", "--temporary"), type = "character", default = NULL, help = "Temporary file path"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output file path"),
  make_option(c("-c", "--cores"), type = "numeric", default = 96, help = "Number of threads available"),
  make_option(c("-v", "--verbose"), type = "logical", default = TRUE, help = "Verbose"),
  make_option(c("-a", "--pattern"), type = "character", default = "_R1_001.fastq", help = "Pattern to parse sample names"),
  make_option(c("-t", "--truncf"), type = "numeric", default = 190, help = "Position at which to trunctate forward reads"),
  make_option(c("-u", "--truncr"), type = "numeric", default = 135, help = "Position at which to trunctate reverse reads")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Access the parsed arguments
path <- opt$path
tmp_path <- opt$temporary
output <- opt$output
verbose <- opt$verbose
num_threads <- opt$cores
pattern_f <- opt$pattern
pattern_r <- gsub("R1", "R2", pattern_f)
trunc_f <- opt$truncf
trunc_r <- opt$truncr

#### Set number of threads ###############################################################
options(mc.cores = num_threads)
cat("Number of threads:", getOption("mc.cores"), "\n")

#### Read in the files #####################################################################
if (verbose == TRUE) print("Reading file and sample names...")

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=pattern_f, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=pattern_r, full.names = TRUE)) 

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#### Inspect read quality profiles ##########################################################
if (verbose == TRUE) print("Inspect read quality profiles...")

p1<-plotQualityProfile(fnFs[1:2], aggregate=TRUE) +
  geom_hline(yintercept = 27) +
  labs(title="R1")
p2<-plotQualityProfile(fnRs[1:2], aggregate=TRUE) +
  geom_hline(yintercept = 27) +
  labs(title="R2") 

pdf(paste0(tmp_path, "/quality_profiles.pdf"))
p1
p2
dev.off()

if (verbose == TRUE) print("Quality profiles saved to /intermediate/quality_profiles.pdf")


#### Filter and trim ########################################################################
if (verbose == TRUE) print("Filter and trim...")

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(tmp_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(tmp_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(trunc_f,trunc_r),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=num_threads, verbose=TRUE)
head(out)


file_name <- "trim_stats.tsv"
file_path <- file.path(tmp_path, file_name)
write.table(out, file_path, sep = "\t", row.names = FALSE)

saveRDS(out, paste0(tmp_path, "/trim_stats.rds"))

#### Learn the Error Rates ##################################################################
if (verbose == TRUE) print("Learn the Error Rates...")


errF <- learnErrors(filtFs, multithread=num_threads, verbose=TRUE) # Forward
errR <- learnErrors(filtRs, multithread=num_threads, verbose=TRUE) # Reverse

p3 <- plotErrors(errF, nominalQ=TRUE) +
  labs(title="R1")
p4 <- plotErrors(errR, nominalQ=TRUE) +
  labs(title="R2")

pdf(paste0(tmp_path, "/ErrorRates.pdf"))
p3
p4
dev.off()

if (verbose == TRUE) print("Error plots saved to /intermediate/ErrorRates.pdf")

#### Sample Inference #######################################################################
if (verbose == TRUE) print("Perform sample inference...")

dadaFs <- dada(filtFs, err=errF, multithread=num_threads, verbose=TRUE) # Forward
dadaRs <- dada(filtRs, err=errR, multithread=num_threads, verbose=TRUE) # Reverse

print("Foward reads")
dadaFs[[1]]

print("Reverse reads")
dadaRs[[1]]

#### Merge paired reads #####################################################################
if (verbose == TRUE) print("Merging forward and reverse reads...")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

newdir <- paste0(tmp_path, "/merge_stats")
dir.create(newdir, showWarnings = FALSE)
for (i in seq_along(mergers)) {
  file_name <- paste0(sample.names[i], "_merge_stats.tsv")
  file_path <- file.path(newdir, file_name)
  write.table(mergers[[i]], file_path, sep = "\t", row.names = FALSE)
}

saveRDS(mergers, paste0(tmp_path, "/merge_stats.rds"))

if (verbose == TRUE) print("Merging statistics saved to /intermediate/merge_stats/.")

#### Construct sequence table ###################################################################
if (verbose == TRUE) print("Making a sequence table...")

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Save Sequence table
file_name <- "sequence_table.tsv"
file_path <- file.path(tmp_path, file_name)
write.table(seqtab, file_path, sep = "\t", row.names = FALSE)

saveRDS(seqtab, paste0(tmp_path, "/sequence_table.rds"))

if (verbose == TRUE) print("Sequencing table saved to /intermediate/sequence_table.tsv")

#### Remove Chimeras #########################################################################
if (verbose == TRUE) print("Remove Chimeras...")

# Remove chimeras from table
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=num_threads, verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, paste0(tmp_path, "/sequence_table_nochim.rds"))

# Create fasta of all unique sequences
dna_sequences <- colnames(seqtab.nochim)
sequence_names <- paste0("ASV", 1:length(dna_sequences))
names(dna_sequences) <- sequence_names

dna_string_set <- DNAStringSet(dna_sequences)
writeXStringSet(dna_string_set, file.path(tmp_path, "dna-sequences.fasta"), format = "fasta")

#Save Sequence table
colnames(seqtab.nochim) <- sequence_names
feature.table <- t(seqtab.nochim)

# If there is only one sample, the sample names disappear, so we add them back
if (length(sample.names) == 1) {
  colnames(feature.table) <- sample.names
}

head(feature.table)

file_name <- "feature-table.tsv"
file_path <- file.path(tmp_path, file_name)
write.table(feature.table, file_path, sep = "\t", row.names = TRUE, quote = FALSE)

if (verbose == TRUE) print("Sequencing table without chimeras saved to /intermediate/feature-table.tsv")

#### Track reads through the pipeline ############################################################

getN <- function(x) sum(getUniques(x))

# Check if there is more than one sample
if (length(sample.names) > 1) {
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
} else {
  track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
}

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Save No-Chimera Sequence table
file_name <- "read_stats.tsv"
file_path <- file.path(tmp_path, file_name)
write.table(track, file_path, sep = "\t", row.names = FALSE)

saveRDS(track, paste0(tmp_path, "/read_stats.rds"))

if (verbose == TRUE) print("Read stats saved to /intermediate/read_stats.tsv")


#### Assign Taxonomy ##############################################################################
# Naive Bayesian classifier method
seqtab.nochim <- readRDS(paste0(tmp_path, "/sequence_table_nochim.rds"))
taxa <- assignTaxonomy(colnames(seqtab.nochim), paste0(tmp_path, "/silva_nr99_v138.1_train_set.fa.gz"), verbose = TRUE)

# replace sequences with sequence names
seqnames <- names(dna_sequences)[match(data.frame(dna_sequences)$dna_sequences, rownames(taxa))]
taxa <- data.frame(taxa)
taxa$ID <- seqnames
taxa <- select(taxa, ID, everything())

taxa2 <- taxa
rownames(taxa2) <- NULL

head(taxa2)

# Save Taxonomy
file_name <- "taxonomy.tsv"
file_path <- file.path(tmp_path, file_name)
write.table(taxa2, file_path, sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(taxa, paste0(tmp_path, "/taxonomy.rds"))

if (verbose == TRUE) print("Taxonomy saved to /intermediate/taxonomy.tsv")
