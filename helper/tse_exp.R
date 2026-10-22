.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(dplyr)
library(tibble)
library(optparse)
library(readr)
library(tidyr) 

option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Path to expected tse", metavar = "character"))

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dataset)) stop("--dataset / -d must be provided")
d <- opt$dataset

print(paste(d))

# Output path and pipelines
output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"

t <- paste0("/projects/0/gs/danique/pipelines/int_results/routput/tse_idtaxa_", d, "_asv_allpipelines.rds")
tse <- readRDS(t)

e <- paste0("/projects/0/gs/danique/pipelines/int_results/routput/tse_exp_",d, ".rds")
tse_exp <- readRDS(e)  

r <-  paste0("/projects/0/gs/danique/pipelines/recall/recall_all_",d, ".tsv")
recall <-  read_tsv(r) %>% 
    select(1,2)

colnames(recall) <- c("ID", "match")
recall$ID2 <- gsub("^\\d+_", "", recall$ID)

if (d %in% "mock12") {
    recall <- recall %>% 
        mutate(ID2 = sapply(strsplit(ID2, "_"), function(x) paste0(x[1:2], collapse="_")))
}

recall2 <- recall %>%
  mutate(match = strsplit(match, ";")) %>%
  unnest(match) %>%
  mutate(match = trimws(match)) %>%
  distinct(ID2, match) 

df <- lapply(rownames(tse), function(rn) {
    newmatch <- recall2$ID2[recall2$match %in% rn]
    if (length(newmatch) == 0) {
        newmatch <- rn # keep old id for non matching 
        d <- data.frame(rn = rn, newmatch=newmatch, weight=1)
    } else if (length(newmatch) > 1) {
        # print(paste(rn, newmatch))
        d <- data.frame(rn = rep(rn, length(newmatch)), newmatch = newmatch, weight = 1/length(newmatch))
    } else {
        d <- data.frame(rn = rn, newmatch=newmatch, weight = 1)
    }
}) %>% bind_rows()

# Link expected IDs to observed
tax <- rowData(tse) %>% data.frame() %>% rownames_to_column(var = "rn")
feat <- assay(tse) %>% data.frame() %>% rownames_to_column(var = "rn")
meta <- colData(tse) %>% data.frame() 

tax2 <- df %>%
    dplyr::left_join(., tax, by = "rn") 
rownames(tax2) <- make.names(tax2$rn, unique=TRUE)  
tax2$rn <- NULL    
    
feat2 <- df %>% select(rn) %>% 
    dplyr::left_join(., feat, by = "rn")
rownames(feat2) <- rownames(tax2)
feat2$rn <- NULL

# Scale the abundance for ambiguous ASVs that occur multiple times (fraction)
feat3 <- feat2 %>% 
    mutate(across(where(is.numeric), ~ . * tax2$weight)) %>%
    mutate(newmatch = tax2$newmatch) %>% 
    group_by(newmatch) %>% 
    summarise(across(where(is.numeric), ~ sum(.))) %>%
    mutate(across(where(is.numeric), ~ . / sum(.)))

# Add expected 
feat_exp <- assay(tse_exp) %>% data.frame() %>%
    rownames_to_column(var="newmatch") %>%
    mutate(newmatch = gsub("^\\d+_", "", newmatch)) 
meta_exp <- colData(tse_exp) %>% data.frame()

feat_comb <- feat3 %>% 
    dplyr::full_join(., feat_exp, by = "newmatch")
feat_comb[is.na(feat_comb)] <- 0

# Agglomerate duplicates (for several IDs matching the same expected) 
feat_comb <- feat_comb %>% 
    group_by(newmatch) %>% 
    summarize(across(where(is.numeric), ~ sum(.))) 

feat_comb <- feat_comb %>% column_to_rownames(var = "newmatch")

meta_comb <- meta %>% 
    bind_rows(meta_exp)

tax_comb <- data.frame(id = rownames(feat_comb))
rownames(tax_comb) <- rownames(feat_comb)

tse_comb <- TreeSummarizedExperiment(
    assays =  SimpleList(relabundance = as.matrix(feat_comb)),
    colData = DataFrame(meta_comb),
    rowData = DataFrame(tax_comb))

sfile <- paste0(output_path, "tse_exp_", d, "_asv_allpipelines.rds")
saveRDS(tse_comb, file = sfile)
