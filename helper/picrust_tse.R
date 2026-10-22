.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(dplyr)
library(tibble)
library(optparse)

option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Dataset name", metavar = "character"),
  make_option(c("-c", "--classifier"), type = "character", default = NULL,
              help = "Classifier name", metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dataset)) stop("--dataset / -d must be provided")
d <- opt$dataset
lev <- opt$classifier

# Output path and pipelines
output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
pipelines <- c("NGTax2", "eDNAFlow", "MetaWorks", "VSEARCH", "APSCALE", "AMPtk", "QIIME2Deblur", "LotuS2", "Natrix", "DADAist2", "dadasnake", 
        "nfcore_ampliseq", "QIIME2", "Tourmaline", "DADA2", "LotuS2_otu")
   
alloutput <- readRDS(file = paste0(output_path, "picrust_combined.rds"))

datasets <- c("mock5", "mock12", "miseq_v4_sim", "novaseq_v4_sim", "gemsim", "mindset", "parkinson")

for (d in datasets) {
    tse_all <- list() # save all tse at asv level 
    feat_tot <- data.frame(func = character())
    meta_tot <- data.frame(sample = character(),
                        pipeline = character(),
                        ID = character())

    for (p in pipelines) { 

        print(paste0("Dataset: ", d, " | Level: ", lev, " | Pipeline: ", p))
        
        if (is.null(alloutput[[d]][[p]][[lev]])) next
        feat <- alloutput[[d]][[p]][[lev]] %>%
            column_to_rownames(var="func")

        tax <- data.frame(func = rownames(feat))
        rownames(tax) <- rownames(feat)

        meta <- data.frame(sample = colnames(feat),
                        pipeline = rep(p, ncol(feat))) 

        colnames(feat) <- paste0(p, "_", colnames(feat))
        rownames(meta) <- colnames(feat)
        meta$ID <- colnames(feat)

        print(identical(rownames(tax), rownames(feat)))

        tse <- TreeSummarizedExperiment(
            assays =  SimpleList(counts = as.matrix(feat)),
            colData = DataFrame(meta),
            rowData = DataFrame(tax)) # asv-level tse

        tse_all[[p]] <- tse # add to main list

        feat_tot <- feat_tot %>%
            dplyr::full_join(., rownames_to_column(feat, var="func"), by = "func")

        meta_tot <- meta_tot %>% add_row(meta)  
    }

    # Combine into one large tse 
    feat_tot <- feat_tot %>%
        column_to_rownames("func")

    feat_tot[is.na(feat_tot)] <- 0

    rownames(meta_tot) <- meta_tot$ID
    print(identical(colnames(feat_tot), rownames(meta_tot)))

    tax_tot <- data.frame(func = rownames(feat_tot))
    rownames(tax_tot) <- rownames(feat_tot)

    print(identical(rownames(feat_tot), rownames(tax_tot)))

    tse_full <- TreeSummarizedExperiment(
    assays =  SimpleList(counts = as.matrix(feat_tot)),
    colData = DataFrame(meta_tot),
    rowData = DataFrame(tax_tot)) 

    assays(tse_full)
    dim(tse_full)

    saveRDS(tse_full, file = paste0(output_path, "tse_", lev, "_", d, "picrust_allpipelines.rds"))
    saveRDS(tse_all, file = paste0(output_path, "tse_", lev, "_", d, "picrust_list_allpipelines.rds"))

}
