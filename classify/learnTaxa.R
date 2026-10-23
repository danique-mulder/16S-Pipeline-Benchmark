
# Trains the IDTAXA classifiers (LearnTaxa) on the V4 reference sets for SILVA and Greengenes2.

.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))

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

db <- c("silva", "gg") 

for (d in db) {
  # Silva
  if (d %in% "silva") {
    fasta <- "/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/silva-138-99-seqs-515-806.fasta"
    taxonomy <- "/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/silva-138-99-tax-515-806.tsv"

    outrds <- "/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/idtaxa-trainset-silva-138-99-515-806.rds"
    outRData <- "/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/trainset.RData"

  } else {
    fasta <- "/gpfs/work2/0/gs/danique/pipelines/ref/greengenes2_sequences.fasta"
    taxonomy <- "/gpfs/work2/0/gs/danique/pipelines/ref/greengenes2_taxonomy_v4.tsv"

    outrds <- "/gpfs/work2/0/gs/danique/pipelines/ref/greengenes2_sequences.rds"
    outRData <- "/gpfs/work2/0/gs/danique/pipelines/ref/trainset_gg.RData"
  }

  # load the sequences from the file
  seqs <- readDNAStringSet(fasta)
  seqs <- RemoveGaps(seqs) # remove gaps 

  tax_id <- read.table(taxonomy, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  identical(names(seqs),tax_id$Feature.ID) # Check if sequence IDs and taxonomy IDs are the same/same order

  tax <- paste0("Root; ", tax_id$Taxon)

  comp <- sapply(strsplit(tax, "; "), function(x) { any(grepl("s__", x))})
  sum(!comp) # any entries without species level? 

  trainset <- LearnTaxa(seqs, tax)
  trainset$ranks <- c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  trainingSet <- trainset

  saveRDS(trainingSet, file = outrds)
  save(trainingSet, file = outRData)

}
