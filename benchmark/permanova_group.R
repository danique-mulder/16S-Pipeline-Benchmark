# Tests beta-diversity–phenotype associations (PERMANOVA, distance ~ phenotype) for each pipeline and for the ground truth.

.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(optparse)
library(tidyverse)
library(tibble)
library(ggpubr)
library(vegan)

option_list <- list(
  make_option(c("-c", "--classifier"), type = "character", default = NULL,
              help = "Classifier", metavar = "character"),
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Dataset name", metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))

# Check required arguments
if (is.null(opt$dataset)) {
  stop("--dataset (-d) must be provided")
}

if (is.null(opt$classifier)) {
  stop("--classifier (-c) must be provided (asv, q2, idtaxa or deftax)")
}

c <- opt$classifier
d <- opt$dataset

output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
pipelines <- c("NGTax2", "eDNAFlow", "MetaWorks", "VSEARCH", "APSCALE", "QIIME2Deblur", "AMPtk", "LotuS2", "Natrix", 
               "DADAist2", "dadasnake", "nfcore_ampliseq", "QIIME2", "Tourmaline", "DADA2", "LotuS2_otu", "Expected")

#### Run the permanova
set.seed(1234)
perm <- list()

if (c %in% "deftax") {
    tse <- readRDS(file = paste0(output_path, "tse_deftax_", d, "_meta_allpipelines.rds"))
} else if (c %in% c("q2", "idtaxa")) {
    tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_genus_meta_allpipelines.rds"))
} else if (c %in% "asv") {
    tse <- readRDS(file = paste0(output_path, "tse_q2_", d, "_asv_meta_alpha_allpipelines.rds"))
    if (d %in% c("miseq_v4_sim", "novaseq_v4_sim")) {
        tse_exp <- readRDS(file = paste0(output_path, "tse_exp_", d, ".rds"))
        meta <- colData(tse_exp) %>% as.data.frame()
        meta2 <- colData(tse) %>% as.data.frame() %>% select(sample, label) %>% unique()
        meta <- meta %>% left_join(meta2, by = "sample")
        colData(tse_exp) <- DataFrame(meta)
    }    
} 

for (p in pipelines) {
    if (c %in% "deftax") {
        if (is.null(tse[[p]])) next
        tse_sub <- tse[[p]]
        tse_sub <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance", pseudocount=0)
        tse_sub <- transformAssay(tse_sub, assay.type = "relabundance", method = "clr", pseudocount=TRUE)
    } else {
        if (p %in% "Expected" & d %in% c("miseq_v4_sim", "novaseq_v4_sim")) {
           tse_sub <- tse_exp
        } else {
           tse_sub <- tse[, tse$pipeline %in% p]
        }
        if (ncol(tse_sub) == 0) next 
        tse_sub <- tse_sub[rowSums(assay(tse_sub), na.rm = TRUE) > 0, ] # remove absent genera from tse
        if (!assayNames(tse_sub) %in% "relabundance") tse_sub <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance", pseudocount=0)
        tse_sub <- transformAssay(tse_sub, assay.type = "relabundance", method = "clr", pseudocount=TRUE)
    }

    print(paste(d, "|", c, "|", p))

    var <- "label"
    if (d %in% "mindset") var <- c(var, "RC1")
    
    for (mydist in c("aitchison", "jaccard")) {
        if (mydist %in% "aitchison") {
            atype = "clr"
            mymethod = "euclidean"
        } else {
            atype = "relabundance"
            mymethod = "jaccard"
        }

        for (v in var) {
            tse_sub <- tse_sub[, !is.na(tse_sub[[v]])]
            print(dim(tse_sub))
            perm[[mydist]][[v]][[p]] <- getPERMANOVA(
                tse_sub,
                assay.type = atype,
                formula = as.formula(paste("x ~", v)),
                method = mymethod,
                permutations = 999)
        }

        print(paste("Permanova for", p, mydist, "ran successfully"))
    }
}

print("saving permanova stats")
saveRDS(perm, file = paste0(output_path, "permanova_group_", d, "_", c, ".rds"))


    