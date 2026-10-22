.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia, quietly=TRUE)
library(optparse, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(tibble, quietly=TRUE)


option_list <- list(
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Dataset name", metavar = "character"),
  make_option(c("-e", "--expected"), type = "character", default = NULL,
              help = "Including expected values?", metavar = "character"),
  make_option(c("-l", "--level", type = "character", default = NULL,
              help = "Level: asv (asv, asv199, asv249), classifier (q2, idtaxa, silva-gg, deftax), picrust (KO/cyc)", metavar = "character")))

opt <- parse_args(OptionParser(option_list = option_list))

# Check required arguments
if (is.null(opt$dataset) || is.null(opt$level) || is.null(opt$expected)) {
  stop("Both --dataset (-d), --expected (-e), and --level (-l) must be provided")
}

if (!opt$level %in% c("asv", "asv249", "asv199", "q2", "idtaxa", "q2_gg", "idtaxa_gg", "deftax", "KO", "cyc")) {
  stop("--level (-l) must be one of: asv, asv249, asv199, q2, idtaxa, q2_gg, idtaxa_gg, deftax, KO, cyc")
}

d <- opt$dataset
c <- opt$level
e <- opt$expected

output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
avg_corr <- list()
avg_union_corr <- list()
avg_part_corr <- list()

if (e == 0) {
    if (c %in% "asv") {
        tse <- readRDS(file = paste0(output_path, "tse_q2_", d, "_asv_allpipelines.rds"))
    } else if (c %in% c("q2", "idtaxa", "q2_gg", "idtaxa_gg", "deftax")) {
        tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_genus_allpipelines.rds"))
    } else if (c %in% c("KO", "cyc")) {
        tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "picrust_allpipelines.rds"))
    } else if (c %in% c("asv199", "asv249")) {
        trim <- gsub("asv", "", c)
        tse <- readRDS(file = paste0(output_path, "tse_q2_", d, "_trim", trim, "_asv_allpipelines.rds")) 
    }
} else {
    if (c %in% "asv") {
        tse <- readRDS(file = paste0(output_path, "tse_exp_", d, "_asv_allpipelines.rds"))
    } else if (c %in% c("q2", "idtaxa")) {
        tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_genus_exp_allpipelines.rds"))  
    } else if (c %in% c("KO", "cyc")) {
        tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_picrust_exp_allpipelines.rds"))  
    } 
}

sample_list <- split(colnames(tse), tse$sample)
transformations <- c("relabundance", "clr", "log")

p <- ""
objname <- paste0("tse", p)
tse_tmp <- get(objname)

if (!any(assayNames(tse) %in% "relabundance")) {
    tse_tmp <- transformAssay(tse_tmp, assay.type = "counts", method = "relabundance", pseudocount = FALSE)
}

tse_tmp <- transformAssay(tse_tmp, assay.type = "relabundance", method = "clr", pseudocount = TRUE)
tse_tmp <- transformAssay(tse_tmp, assay.type = "relabundance", method = "log", pseudocount = TRUE)

assign(objname, tse_tmp)

p2 <- ifelse(p %in% "", "all", gsub("_", "", p)) 

for (t in transformations) {

    print(paste(c, d, p2, t))

    corr_list <- lapply(sample_list, function(s) {
        mat_tmp <- assay(tse_tmp, t)
        mat_tmp <- mat_tmp[, s, drop=FALSE]
        
        mat_tmp_abs <- assay(tse_tmp, "relabundance")
        mat_tmp_abs <- mat_tmp_abs[, s, drop=FALSE]

        # Compute correlation between pipelines 
        n_pipelines <- ncol(mat_tmp)
        
        # Initialize both correlation matrices
        cor_nonzero <- matrix(NA, nrow = n_pipelines, ncol = n_pipelines) # either or both of the pipeline are >0
        cor_union <- matrix(NA, nrow = n_pipelines, ncol = n_pipelines) # on only shared features
        cor_partial <- matrix(NA, nrow = n_pipelines, ncol = n_pipelines) # adjust for presence/absence patterns
        rownames(cor_nonzero) <- colnames(cor_nonzero) <- colnames(mat_tmp)
        rownames(cor_union) <- colnames(cor_union) <- colnames(mat_tmp)
        rownames(cor_partial) <- colnames(cor_partial) <- colnames(mat_tmp)
        
        for (i in 1:n_pipelines) {
            for (j in i:n_pipelines) {
                # Get features where at least one pipeline has non-zero
                mask <- (mat_tmp_abs[, i] != 0) | (mat_tmp_abs[, j] != 0)
                cor_val <- cor(mat_tmp[mask, i], mat_tmp[mask, j], 
                                method = "pearson", use = "complete.obs")
                cor_nonzero[i, j] <- cor_val
                cor_nonzero[j, i] <- cor_val

                #  Union correlation
                # Get features where both pipelines have non-zero
                mask <- (mat_tmp_abs[, i] != 0) & (mat_tmp_abs[, j] != 0)
                if (sum(mask) == 0) {
                    print(paste("No overlapping features between", colnames(mat_tmp)[i],"and",  colnames(mat_tmp)[j]))
                    cor_union[i,j] <- NA
                    cor_union[j,i] <- NA
                } else {
                    cor_val <- cor(mat_tmp[mask, i], mat_tmp[mask, j], 
                                    method = "pearson", use = "complete.obs")
                    cor_union[i, j] <- cor_val
                    cor_union[j, i] <- cor_val
                }
                
                # Partial correlation controlling for sparsity
                if (!t %in% "pa") {
                    # Create sparsity indicators (using absolute counts)
                    sparse_i <- as.numeric(mat_tmp_abs[, i] == 0)
                    sparse_j <- as.numeric(mat_tmp_abs[, j] == 0)
                    
                    # Need variation in both variables and sparsity indicators
                    if (length(unique(mat_tmp[, i])) > 1 && 
                        length(unique(mat_tmp[, j])) > 1 &&
                        length(unique(sparse_i)) > 1 && 
                        length(unique(sparse_j)) > 1) {
                    
                    # Residualize against sparsity patterns
                    resid_i <- residuals(lm(mat_tmp[, i] ~ sparse_i + sparse_j))
                    resid_j <- residuals(lm(mat_tmp[, j] ~ sparse_i + sparse_j))
                    
                    # Correlate residuals
                    cor_val <- cor(resid_i, resid_j, method = "pearson", use = "complete.obs")
                    cor_partial[i, j] <- cor_val
                    cor_partial[j, i] <- cor_val
                    }
                } else { 
                    cor_partial <-  NA
                }
            }
        }
        
        list(nonzero_pairs = cor_nonzero, 
             union = cor_union,
             partial_sparsity = cor_partial)
    })

    pipeline_names <- tse_tmp$pipeline %>% unique()
    n <- length(pipeline_names)
    avg_corr[[p2]][[d]][[c]][[t]] <- matrix(NA, n, n, dimnames=list(pipeline_names, pipeline_names))    
    avg_union_corr[[p2]][[d]][[c]][[t]] <- matrix(NA, n, n, dimnames=list(pipeline_names, pipeline_names))
    avg_part_corr[[p2]][[d]][[c]][[t]] <- matrix(NA, n, n, dimnames=list(pipeline_names, pipeline_names))

    for(i in 1:n) {
        for(j in i:n) {

            cor_values <- sapply(names(corr_list), function(m) ifelse(!nrow(corr_list[[m]][["nonzero_pairs"]]) < j, corr_list[[m]][["nonzero_pairs"]][i,j], NA))
            avg_corr[[p2]][[d]][[c]][[t]][i,j] <- mean(cor_values, na.rm=TRUE)
            avg_corr[[p2]][[d]][[c]][[t]][j,i] <- avg_corr[[p2]][[d]][[c]][[t]][i,j]

            cor_union_values <- sapply(names(corr_list), function(m) ifelse(!nrow(corr_list[[m]][["union"]]) < j, corr_list[[m]][["union"]][i,j], NA))
            if (any(!is.na(cor_union_values))) {
                avg_union_corr[[p2]][[d]][[c]][[t]][i,j] <- mean(cor_union_values, na.rm=TRUE)
                avg_union_corr[[p2]][[d]][[c]][[t]][j,i] <- avg_union_corr[[p2]][[d]][[c]][[t]][i,j]
            } else {
                avg_union_corr[[p2]][[d]][[c]][[t]][i,j] <- NA
                avg_union_corr[[p2]][[d]][[c]][[t]][j,i] <- NA
            }

            cor_part_values <- sapply(names(corr_list), function(m) ifelse(!nrow(corr_list[[m]][["partial_sparsity"]]) < j, corr_list[[m]][["partial_sparsity"]][i,j], NA))
            avg_part_corr[[p2]][[d]][[c]][[t]][i,j] <- mean(cor_part_values, na.rm=TRUE)
            avg_part_corr[[p2]][[d]][[c]][[t]][j,i] <- avg_part_corr[[p2]][[d]][[c]][[t]][i,j]
        }
    }

} # transformation

cor_list <- list(avg_corr = avg_corr,
                 avg_union_corr = avg_union_corr,
                 avg_part_corr = avg_part_corr)

exp <- ifelse(e == 1 & c %in% "asv", "_exp", "")
saveRDS(cor_list, file = paste0(output_path, "cor_list", c,exp, "_", d, ".rds"))