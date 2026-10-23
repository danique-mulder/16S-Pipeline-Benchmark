# Applies the ASV-level filters (prevalence, total count, minimum and maximum length) 
# and recalculates Jaccard similarity and correlations after each filtering step.

.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(optparse)
library(tidyverse)
library(tibble)
library(ggpubr)
library(showtext)
library(vegan)
library(Matrix)

output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
sims <- c("miseq_v4_sim", "novaseq_v4_sim")
reals <- c("gemsim", "mindset", "parkinson")

pipelines <- c("Expected", "NGTax2", "eDNAFlow", "MetaWorks", "VSEARCH", "APSCALE", "QIIME2Deblur", "AMPtk", "LotuS2", "Natrix", "DADAist2", "dadasnake", 
        "nfcore_ampliseq", "QIIME2", "Tourmaline", "DADA2", "LotuS2_otu")

getabstat <- 0
filttab <- 0
calcmet <- 1

#############################################################
# init function

jaccard_sparse <- function(m) {
  # coerce to sparse column-compressed binary
  m <- as(m, "CsparseMatrix")

  A   <- crossprod(m) # interaction size (set_i ∩ set_j)
  s   <- Matrix::diag(A) # pipeline set size (set_i...n)

  # work on the non-zero entries only
  A_t  <- as(A, "TsparseMatrix") # triplet form (i=rows, j=columns, x=value)
  A_t@x <- A_t@x / (s[A_t@i + 1L] +
                      s[A_t@j + 1L] -
                      A_t@x)

  J <- as(A_t, "CsparseMatrix")
  diag(J) <- 1
  J
}

calcmetrics <- function(feat) {

    # Get pairwise Jaccard
    feat_sub <- feat[, !grepl("Expected|NGTax2|LotuS2|LotuS2_otu", colnames(feat))]
    feat_pa <- feat_sub
    feat_pa[feat_pa > 0] <- 1 # pa transform

    mat <- jaccard_sparse(feat_pa) %>% as.matrix()
    
    patterns <- paste0(pipelines, collapse = "|")  # "pipeline1|pipeline2|...|pipeline14"
    pipe_names <- gsub(paste0("^(.*(", patterns, "))(.*)$"), "\\1", rownames(mat))
    unique_pipelines <- unique(pipe_names)
    pipe_factor <- factor(pipe_names, levels = unique_pipelines)
    n_pipe <- length(unique_pipelines)

    patterns2 <- paste0(paste0(pipelines, "_"), collapse = "|")
    sample_names <- gsub(patterns2, "", rownames(mat))
    unique_samples <- unique(sample_names)

    result <- matrix(0, n_pipe, n_pipe, dimnames = list(unique_pipelines, unique_pipelines))
    for (k in seq_len(n_pipe)) {
        for (l in k:n_pipe) { # upper triangle + diagonal only
            pair_vals <- sapply(unique_samples, function(samp) {
            mat[paste0(unique_pipelines[k], "_", samp),
                paste0(unique_pipelines[l], "_", samp)]
            })
            result[k, l] <- result[l, k] <- mean(pair_vals)
        }
    } 
    
    diag(result) <- NA
    jacc <- rowMeans(result, na.rm = TRUE)

    # Get correlation
    sample_list <- split(colnames(feat_sub), sample_names)
    feat_clr <- t(compositions::clr(t(feat_sub + 1e-6)))

    # Fix pseudo-count first
    corr <- lapply(names(sample_list), function(samp) {
        s <- sample_list[[samp]]
        clr_tmp <- feat_clr[, s, drop = FALSE]   # features × pipelines
        abs_tmp <- feat_sub[, s, drop = FALSE]

        # Precompute all masks as a logical matrix (features × pipelines)
        nonzero <- abs_tmp != 0                  # logical, features × pipelines

        n <- ncol(clr_tmp)
        cor_mat <- matrix(1, n, n, dimnames = list(colnames(clr_tmp), colnames(clr_tmp)))

        for (i in seq_len(n)) {
            for (j in seq(i + 1L, n, length.out = max(0, n - i))) {
                # union mask: at least one pipeline non-zero
                mask <- nonzero[, i] | nonzero[, j]
                m <- sum(mask, na.rm = TRUE)

                if (m < 3L) {             # need at least 3 points for correlation
                    cor_mat[i, j] <- cor_mat[j, i] <- NA
                    next
                }

                is <- clr_tmp[mask,i]
                js <- clr_tmp[mask,j]

                cor_mat[i, j] <- cor_mat[j, i] <- cor(is, js,
                    method = "spearman", use = "complete.obs")
            }
        }

        rownames(cor_mat) <- colnames(cor_mat) <- gsub(paste0("_", samp), "", colnames(cor_mat))
        cor_mat
    })

    corr_array <- simplify2array(corr)
    avg_cor_mat <- apply(corr_array, 1:2, mean, na.rm = TRUE)
    diag(avg_cor_mat) <- NA

    df_cor <- data.frame(pipeline = rownames(avg_cor_mat),
                         cor = rowMeans(avg_cor_mat, na.rm=TRUE))

    df <- data.frame(pipeline = names(jacc),
                     jaccard=jacc) %>% 
          dplyr::left_join(df_cor, by = "pipeline")
    return(df)
}

#############################################################
####### Get abundance statistics

if (getabstat) {
    ab_stats <- list()
    seq_length <- read_tsv("/projects/0/gs/danique/pipelines/int_results/seq_length.tsv")
    seq_length$ID <- paste0(seq_length$Pipeline, "_", seq_length$ID)

    # Collect information about each feature (seq length, sum of relative abundance, sum of counts, prevalence etc.)
    for (d in c(sims, reals)) { 
    
    tse_asv <- readRDS(file = paste0(output_path, "tse_q2_", d, "_trim249_asv_allpipelines.rds"))

    lev <- "asv"
    tmp <- get(paste0("tse_", lev))
    ave_ab <- data.frame(pipeline = character(),
                        id = character(),
                        seq_length = numeric(),
                        sum = numeric(),
                        sum_count = numeric(),
                        max_relab = numeric(),
                        max_count = numeric(),
                        count_u5 = numeric(), # <= 5 is TRUE, >5 is FALSE 
                        prev_n = numeric(), # number of samples with count > 0
                        prev_1 = numeric(), # present in 1=TRUE, present in >1 sample=FALSE
                        ave = numeric(),
                        prev = numeric())
    
    seq_sub <- seq_length %>% filter(Dataset %in% d)

    for (p in pipelines) {
        
        if (p %in% "Expected") next
        print(paste("Dataset:", d, "Pipeline:", p, "Level:", lev))

        tmp_sub <- tmp[,tmp$pipeline %in% p]
        tmp_sub <- tmp_sub[rowSums(assay(tmp_sub), na.rm = TRUE) > 0, ]
        
        if (is.null(tmp_sub)) next

        tax <- rowData(tmp_sub) %>% data.frame() %>% 
        rownames_to_column("ID") %>%
        dplyr::left_join(., seq_sub, by = "ID") %>%
        column_to_rownames("ID")
        
        rowData(tmp_sub) <- DataFrame(tax)      
        
        if (ncol(tmp_sub) == 0) next
        tmp_sub <- transformAssay(tmp_sub, assay.type = "counts", method = "relabundance", pseudocount = FALSE)

        new_row <- data.frame(pipeline = p,
                            id = rownames(tmp_sub), 
                            seq_length = rowData(tmp_sub)$Length,
                            sum = rowSums(assay(tmp_sub, "relabundance")),
                            sum_count = rowSums(assay(tmp_sub, "counts")),
                            max_relab = apply(assay(tmp_sub, "relabundance"), 1, max),
                            max_count = apply(assay(tmp_sub, "counts"), 1, max),
                            count_u5 = rowSums(assay(tmp_sub, "counts")) <= 5,
                            ave = rowMeans(assay(tmp_sub, "relabundance")),
                            prev = rowSums(assay(tmp_sub) > 0) / ncol(tmp_sub) * 100,
                            prev_n = rowSums(assay(tmp_sub) > 0),
                            prev_1 = rowSums(assay(tmp_sub) > 0) == 1)
        
        ave_ab <- ave_ab %>%
        add_row(new_row)
        
    }
    ab_stats[[d]][[lev]] <- ave_ab
    }

    # Save 
    saveRDS(ab_stats, file = paste0(output_path, "ab_stats_real.rds"))
}

###### Filter the feature tables 
if (filttab) {

    print("Filtering feature tables...")
    ab_stats <- readRDS(file = paste0(output_path, "ab_stats_real.rds"))

    data_df_filt <- list()
    tse_filt <- list()

    for (d in c(sims, reals)) { #   

        tse_asv <- readRDS(file = paste0(output_path, "tse_q2_", d, "_trim249_asv_allpipelines.rds"))
        if (!any(assayNames(tse_asv) %in% "relabundance")) {
            tse_asv <- transformAssay(tse_asv, assay.type = "counts", method = "relabundance", pseudocount = FALSE)
        }

        assay_nofilt <- list()
        assay_nofilt[["1"]] <- assay(tse_asv, "relabundance") 

        assay_prev <- list()
        assay_count <- list()
        assay_short <- list()
        assay_long <- list()

        for (p in unique(ab_stats[[d]][["asv"]]$pipeline)) {

            print(paste("Dataset:", d, "| Pipeline:", p))
            data_df <- ab_stats[[d]][["asv"]] %>% filter(pipeline %in% p) 
            data_df_filt[[d]][[p]][["nofilt"]][["1"]] <- data_df

            # Prevalence filter
            print("Prevalence filter")
            for (i in 2:6) { # present in >-1-6 samples (1=no filtering, 6= present in 6 or more samples)
                i2 <- as.character(i)
                print(i2)
                data_df_filt[[d]][[p]][["prev"]][[i2]] <- data_df %>% filter(prev_n >= i)

                if (is.null(assay_prev[[i2]])) { 
                    assay_prev[[i2]]  <- assay_nofilt[["1"]]
                } 

                row_idx <- !rownames(assay_prev[[i2]]) %in% data_df_filt[[d]][[p]][["prev"]][[i2]]$id
                col_idx <- grepl(p, colnames(assay_prev[[i2]]))
                assay_prev[[i2]][row_idx, col_idx] <- 0
            }

            # Count filter 
            print("Count filter")
            for (i in 2:6) { # total count between 1-10 (1==1 or more = no filter)
                i2 <- as.character(i)
                print(i2)
                data_df_filt[[d]][[p]][["count"]][[i2]] <- data_df %>% filter(sum_count >= i)

                if (is.null(assay_count[[i2]])) { 
                    assay_count[[i2]]  <- assay_nofilt[["1"]]
                } 

                row_idx <- !rownames(assay_count[[i2]]) %in% data_df_filt[[d]][[p]][["count"]][[i2]]$id
                col_idx <- grepl(p, colnames(assay_count[[i2]]))
                assay_count[[i2]][row_idx, col_idx] <- 0

            }

            # Short seq
            print("Filter short sequences")
            for (i in 1:1) { # no filter (0) vs filter sequnces < 240 (1) 
                i2 <- as.character(i)
                print(i2)

                if (p %in% c("NGTax2", "AMPtk", "LotuS2", "LotuS2_otu", "QIIME2_Deblur")) { # These are pipelines that trim sequences to single length (no variation)
                    data_df_filt[[d]][[p]][["short"]][[i2]] <- data_df # no filtering, so continue 
                    next
                }

                data_df_filt[[d]][[p]][["short"]][[i2]] <- data_df %>% filter(seq_length > 240)

                if (is.null(assay_short[[i2]])) { 
                    assay_short[[i2]]  <- assay_nofilt[["1"]]
                } 

                row_idx <- !rownames(assay_short[[i2]]) %in% data_df_filt[[d]][[p]][["short"]][[i2]]$id
                col_idx <- grepl(p, colnames(assay_short[[i2]]))
                assay_short[[i2]][row_idx, col_idx] <- 0

            }

            # Long seq 
            print("Filter long sequences")
            for (i in 1:1) { # no filter (0) vs filter sequnces < 240
            
                i2 <- as.character(i)
                print(i2)

                if (p %in% c("NGTax2", "AMPtk", "LotuS2", "LotuS2_otu", "QIIME2_Deblur")) { # These are pipelines that trim sequences to single length (no variation)
                    data_df_filt[[d]][[p]][["long"]][[i2]] <- data_df # no filtering, so continue 
                    next
                }

                data_df_filt[[d]][[p]][["long"]][[i2]] <- data_df %>% filter(seq_length < 260)

                if (is.null(assay_long[[i2]])) { 
                    assay_long[[i2]]  <- assay_nofilt[["1"]]
                } 

                row_idx <- !rownames(assay_long[[i2]]) %in% data_df_filt[[d]][[p]][["long"]][[i2]]$id
                col_idx <- grepl(p, colnames(assay_long[[i2]]))
                assay_long[[i2]][row_idx, col_idx] <- 0
            }
        }

        saveRDS(assay_nofilt, file = paste0(output_path, "assay_nofilt_", d, ".rds"))
        saveRDS(assay_prev, file = paste0(output_path, "assay_prev_", d, ".rds"))
        saveRDS(assay_count, file = paste0(output_path, "assay_count_", d, ".rds"))
        saveRDS(assay_short, file = paste0(output_path, "assay_short_", d, ".rds"))
        saveRDS(assay_long, file = paste0(output_path, "assay_long_", d, ".rds"))

    }

    saveRDS(data_df_filt, file = paste0(output_path, "filtered_df.rds"))
}


##### Calculate metrics 
if (calcmet) {
    
    data_df_filt <- readRDS(file = paste0(output_path, "filtered_df.rds"))
    ab_stats <- readRDS(file = paste0(output_path, "ab_stats_real.rds"))

    # Calculate metrics 
    for (d in c(reals)) {

        # tse_filt <- readRDS(file = paste0(output_path, "filtered_tse.rds"))
        assay_nofilt <- readRDS(file = paste0(output_path, "assay_nofilt_", d, ".rds"))
        assay_prev <- readRDS(file = paste0(output_path, "assay_prev_", d, ".rds"))
        assay_count <- readRDS(file = paste0(output_path, "assay_count_", d, ".rds"))
        assay_short <- readRDS(file = paste0(output_path, "assay_short_", d, ".rds"))
        assay_long <- readRDS(file = paste0(output_path, "assay_long_", d, ".rds"))

        for (m in c("nofilt", "prev", "count", "short", "long")) {

            prev_df <- data.frame(dataset = character(),
                pipeline = character(),
                filttype = character(),
                threshold = numeric(),
                jaccard = numeric(),
                cor = numeric())

            feat_all <- get(paste0("assay_", m))
            for (i2 in names(feat_all)) {

                threshold <- as.numeric(i2)
                print(paste("Dataset:", d, " | Metric:", m, " | Index:", i2))
                
                feat <- feat_all[[i2]] 

                # calculate correlation and jaccard
                metrics <- calcmetrics(feat) ### Revise for matrix 
                metrics <- metrics %>%
                    mutate(dataset = d,
                        filttype = m,
                        threshold = threshold)

                prev_df <- prev_df %>%
                    bind_rows(metrics)
            }

            saveRDS(prev_df, file = paste0(output_path, "filtering_settings_real_", d, "_", m, ".rds"))
        }
    }
}
