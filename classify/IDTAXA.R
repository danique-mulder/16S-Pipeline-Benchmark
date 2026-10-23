
# Runs the DECIPHER IDTAXA classification on a FASTA file and writes a taxonomy table.

# Load packages
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(c('BiocManager'), quiet = TRUE, update = FALSE, ask = FALSE)
}

if(!requireNamespace("DECIPHER", quietly = TRUE)) {
  BiocManager::install(c('DECIPHER'), quiet = TRUE, update = FALSE, ask = FALSE)
}

require(DECIPHER)
require(optparse)
require(tibble)
require(dplyr)

#### Parse command line arguments 
# Define the options
option_list <- list(
  make_option(c("-s", "--seq"), type = "character", default = NULL, help = "fasta file with sequences"),
  make_option(c("-t", "--trainset"), type = "character", default = NULL, help = "path to trainset"),
  make_option(c("-n", "--nametrain"), type = "character", default = NA, help = "abbreviation of the training set name"),
  make_option(c("-o", "--output"), type = "character", default = NULL, help = "output directory"),
  make_option(c("-c", "--cpu"), type = "integer", default = 24, help = "number of threads to use [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Access the parsed arguments
seq <- opt$seq
trainset <- opt$trainset
outputdir <- opt$output
cpu <- opt$cpu
tname <- opt$nametrain

# Set number of threadqs
options(mc.cores = cpu)
cat("Number of threads:", getOption("mc.cores"), "\n")

# load the sequences from the file
seqs <- readDNAStringSet(seq)

# remove any gaps (if needed)
seqs <- RemoveGaps(seqs)

# load a training set object (trainingSet)
# see http://DECIPHER.codes/Downloads.html
load(trainset)

# classify the sequences
ids <- IdTaxa(seqs,
   trainingSet,
   strand="both", # or "top" if same as trainingSet
   threshold=60, # 60% confidence threshold
   processors=NULL) # use all available processors

# saveRDS(ids, file = file.path(outputdir, "idtaxa_silva_ids.rds"))

ranks <- c("domain", "phylum", "class", "order", "family", "genus") # ranks of interest

# Get dataframe with taxonomic information
taxid <- do.call(rbind, lapply(names(ids), function(x) { 
  id <- x
  tax <- ids[[x]]$taxon[-1]
  length(tax) <- 7 # Ensure tax is length 7, filling with NA if shorter
  c(id, tax) 
})) %>% data.frame()

colnames(taxid) <- c("ID", "domain", "phylum", "class", "order", "family", "genus", "species")

# taxid <- t(sapply(ids, function(x) {
#         m <- match(ranks, x$rank)
#         taxa <- x$taxon[m]
#         # taxa[startsWith(taxa, "unclassified_")] <- NA
#         taxa
# }))

# taxid <- data.frame(taxid) %>% rownames_to_column() 
# colnames(taxid) <- c("ID", ranks)

# Get dataframe with confidence scores
confidence <- do.call(rbind, lapply(names(ids), function(x) { 
  id <- x
  conf <- ids[[x]]$confidence[-1]
  length(conf) <- 7 # Ensure tax is length 7, filling with NA if shorter
  c(id, conf) 
})) %>% data.frame()

colnames(confidence) <- c("ID", "domain", "phylum", "class", "order", "family", "genus", "species")

# confidence <- t(sapply(ids, function(x) {
#         m <- match(ranks, x$rank)
#         conf <- x$confidence[m]
#         # taxa[startsWith(taxa, "unclassified_")] <- NA
#         conf
# }))

# confidence <- data.frame(confidence) %>% rownames_to_column()
# colnames(confidence) <- c("ID", ranks)

head(taxid)
head(confidence)

# Save taxid to a tsv file
file_name <- paste0("idtaxa", ifelse(grep("_otu", seq), "_otu", ""), ifelse(grep("concat", seq), "_concat", ""), ifelse(is.na(tname), "", paste0("_", tname)), ".tsv")
output_file <- file.path(outputdir, file_name)
write.table(taxid, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

file_name <- paste0("idtaxa_confidence", ifelse(grep("_otu", seq), "_otu", ""), ifelse(grep("concat", seq), "_concat", ""), ifelse(is.na(tname), "", paste0("_", tname)), ".tsv")
output_file <- file.path(outputdir, file_name)
write.table(confidence, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

