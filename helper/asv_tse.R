.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(dplyr)
library(tibble)
library(optparse)
library(readr)

option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Dataset name", metavar = "character"),
  make_option(c("-c", "--classifier"), type = "character", default = NULL,
              help = "Classifier name", metavar = "character"),
  make_option(c("-t", "--trim"), type = "numeric", default = NULL,
              help = "Trim length", metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dataset)) stop("--dataset / -d must be provided")
d <- opt$dataset
c <- opt$classifier
trim <- opt$trim

print(paste(d, c, trim))

# Output path and pipelines
output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
algorithms <- "/projects/0/gs/danique/pipelines/int_results/metadata/algorithms.tsv"
pipelines <- c("NGTax2", "eDNAFlow", "MetaWorks", "VSEARCH", "APSCALE", "AMPtk", "QIIME2Deblur", "LotuS2", "Natrix", "DADAist2", "dadasnake", 
        "nfcore_ampliseq", "QIIME2", "Tourmaline", "DADA2", "LotuS2_otu")
   
alloutput <- readRDS(file = paste0(output_path, "all_tables_combined.rds"))
algorithms_df <- read_tsv(algorithms) %>% select(where(~ !all(is.na(.)))) # algorithmic information

#  Load all pipeline output and dereplication table
if (!is.null(trim)) {
     uc <- read.table(file = paste0("/projects/0/gs/danique/pipelines/int_results/", d, "_trim", trim, "_mapping.uc"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
} else {
    uc <- read.table(file = paste0("/projects/0/gs/danique/pipelines/int_results/", d, "_mapping.uc"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
}
colnames(uc) <- c("type", "cluster", "length", "percentage", "strand", "alignment_start", "alignment_end","dk", "original", "new")

# Replace * (original == new) with original sequence name
uc$new[uc$new == "*"] <- uc$original[uc$new == "*"]

# Initialize combined feature and metadata tables
feat_tot <- data.frame(asv = unique(uc$new))
tax_tot <- data.frame(kingdom = character(), 
                      phylum = character(), 
                      class = character(), 
                      order = character(), 
                      family = character(), 
                      genus = character(), 
                      species = character())
meta_tot <- data.frame(sample = character(),
                       pipeline = character(),
                       ID = character(),
                       primer=character(),
                       qc=character(),
                       denoising=character(),
                       chimera=character(),
                       merging=character(),
                       post=character(),
                       tax=character(),
                       db=character())

final_asv <- unique(uc$new)
tot_asv <- length(final_asv)

for (p in pipelines) {

    if (p %in% "NGTax2" && !is.null(trim)) next
    if (p %in% c("LotuS2", "LotuS2_otu") && trim == 249) next

    print(paste(d,"|", p))

    if (is.null(alloutput[[d]][[p]][["feat"]])) next

    feat <- alloutput[[d]][[p]][["feat"]] %>%
        select(where(~ !all(is.na(.))))

    tax <- alloutput[[d]][[p]][[c]]
    tax <- tax[match(rownames(feat), rownames(tax)), , drop = FALSE]
    trn <- rownames(tax)

    tax <- sapply(colnames(tax), function(x) { gsub("^[a-z]__", "", tax[[x]]) }) %>% data.frame() # needed for idtaxa
    rownames(tax) <- trn

    identical(rownames(feat), rownames(tax)) %>% print()

    rn <- paste0(p, "_", rownames(feat))
    nrn <- sapply(rn, function(x) uc$new[uc$original %in% x][1]) # assign new rownames
    
    if (!is.null(trim) && length(unique(nrn)) < length(nrn)) {
        duplicates <- nrn[duplicated(nrn)]
        print(paste(length(duplicates), "duplicates for", p, "after trimming at", trim, "bp. Check taxonomy and agglomerate."))
        for (dup in duplicates) {
            dup_indices <- which(nrn == dup)
            gn <- tax[dup_indices, ]$genus %>% unique() 
            # print(gn)
            if (length(gn) > 1) message(paste("Warning: different tax for duplicated ASV", dup, "in pipeline", p))
        }

        # Agglomerate features with the same new name (same asv sequence). Sum the features
        feat$nrn <- nrn
        feat <- feat %>%
            group_by(nrn) %>%
            summarise(across(everything(), sum), .groups = "drop") %>%
            column_to_rownames("nrn")

        # Take tax from the first asv.
        tax$nrn <- nrn
        tax <- tax %>%
            group_by(nrn) %>%
            summarise(across(everything(), dplyr::first), .groups = "drop") %>%
            column_to_rownames("nrn")
    } else if (!is.null(trim)) {
        print(paste("No duplicates for", p, "after trimming at", trim))
        rownames(feat) <- nrn
        rownames(tax) <- nrn
    } else if (is.null(trim)) {
        rownames(feat) <- nrn
        rownames(tax) <- nrn
    }

    colnames(feat) <- paste0(p, "_", colnames(feat))

    meta <- data.frame(sample = sapply(strsplit(colnames(feat), "_"), tail, 1),
                       pipeline = rep(p, ncol(feat)))
    rownames(meta) <- colnames(feat)
    meta$ID <- colnames(feat)
    meta <- meta %>%
        dplyr::left_join(., algorithms_df, by = "pipeline") # add algorithm information

    feat$asv <- rownames(feat)

    # Combine tables
    feat_tot <- feat_tot %>%
        dplyr::left_join(., feat, by = "asv")
    meta_tot <- meta_tot %>% add_row(meta)

    tax_add <- tax %>% filter(rownames(tax) %in% final_asv)
    colnames(tax_add) <- colnames(tax_tot)
    
    if (nrow(tax_add) > 0) {
        final_asv <- final_asv[-which(final_asv %in% rownames(tax_add))]
        tax_tot <- rbind(tax_tot, tax_add)
    }

    # print(length(final_asv) + nrow(tax_tot))
    if ((length(final_asv) + nrow(tax_tot)) != tot_asv) break

}

# Combine into tse
feat_tot <- feat_tot %>% 
    column_to_rownames("asv") %>%
    as.matrix()

feat_tot[is.na(feat_tot)] <- 0

rownames(meta_tot) <- meta_tot$ID
identical(colnames(feat_tot), rownames(meta_tot)) # check if all samples are present in correct order

tax_tot <- tax_tot[match(rownames(feat_tot), rownames(tax_tot)), , drop = FALSE]
identical(rownames(feat_tot), rownames(tax_tot))

tse <- TreeSummarizedExperiment(
    assays =  SimpleList(counts = feat_tot),
    colData = DataFrame(meta_tot),
    rowData = DataFrame(tax_tot))

tse <- tse[rowSums(assay(tse), na.rm=TRUE) > 0, ] # remove all features that are absent from all samples (rowsums==0)

if (is.null(trim)) { 
    trim_char <- ""
} else {
    trim_char <- paste0("_trim", trim)
}

saveRDS(tse, file = paste0(output_path, "tse_", c, "_", d, trim_char, "_asv_allpipelines.rds"))
