# Calculates abundance correlations between pipelines and against the expected community, 
# for both the union and the intersection sets.

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
avg_intersection_corr <- list()
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
        cor_nonzero <- matrix(NA, nrow = n_pipelines, ncol = n_pipelines) # either or both of the pipeline are >0 --> cor (union)
        cor_intersection <- matrix(NA, nrow = n_pipelines, ncol = n_pipelines) # on only shared features --> cor (intersection)
        rownames(cor_nonzero) <- colnames(cor_nonzero) <- colnames(mat_tmp)
        rownames(cor_intersection) <- colnames(cor_intersection) <- colnames(mat_tmp)
        
        for (i in 1:n_pipelines) {
            for (j in i:n_pipelines) {
                # Get features where at least one pipeline has non-zero
                mask <- (mat_tmp_abs[, i] != 0) | (mat_tmp_abs[, j] != 0)
                cor_val <- cor(mat_tmp[mask, i], mat_tmp[mask, j], 
                                method = "pearson", use = "complete.obs")
                cor_nonzero[i, j] <- cor_val
                cor_nonzero[j, i] <- cor_val

                #  Intersection correlation
                # Get features where both pipelines have non-zero
                mask <- (mat_tmp_abs[, i] != 0) & (mat_tmp_abs[, j] != 0)
                if (sum(mask) == 0) {
                    print(paste("No overlapping features between", colnames(mat_tmp)[i],"and",  colnames(mat_tmp)[j]))
                    cor_intersection[i,j] <- NA
                    cor_intersection[j,i] <- NA
                } else {
                    cor_val <- cor(mat_tmp[mask, i], mat_tmp[mask, j], 
                                    method = "pearson", use = "complete.obs")
                    cor_intersection[i, j] <- cor_val
                    cor_intersection[j, i] <- cor_val
                }
                
                          }
        }
        
        list(nonzero_pairs = cor_nonzero, 
             intersection = cor_intersection)
    })

    pipeline_names <- tse_tmp$pipeline %>% unique()
    n <- length(pipeline_names)
    avg_corr[[p2]][[d]][[c]][[t]] <- matrix(NA, n, n, dimnames=list(pipeline_names, pipeline_names))    
    avg_intersection_corr[[p2]][[d]][[c]][[t]] <- matrix(NA, n, n, dimnames=list(pipeline_names, pipeline_names))
  
    for(i in 1:n) {
        for(j in i:n) {

            cor_values <- sapply(names(corr_list), function(m) ifelse(!nrow(corr_list[[m]][["nonzero_pairs"]]) < j, corr_list[[m]][["nonzero_pairs"]][i,j], NA))
            avg_corr[[p2]][[d]][[c]][[t]][i,j] <- mean(cor_values, na.rm=TRUE)
            avg_corr[[p2]][[d]][[c]][[t]][j,i] <- avg_corr[[p2]][[d]][[c]][[t]][i,j]

            cor_intersection_values <- sapply(names(corr_list), function(m) ifelse(!nrow(corr_list[[m]][["intersection"]]) < j, corr_list[[m]][["intersection"]][i,j], NA))
            if (any(!is.na(cor_intersection_values))) {
                avg_intersection_corr[[p2]][[d]][[c]][[t]][i,j] <- mean(cor_intersection_values, na.rm=TRUE)
                avg_intersection_corr[[p2]][[d]][[c]][[t]][j,i] <- avg_intersection_corr[[p2]][[d]][[c]][[t]][i,j]
            } else {
                avg_intersection_corr[[p2]][[d]][[c]][[t]][i,j] <- NA
                avg_intersection_corr[[p2]][[d]][[c]][[t]][j,i] <- NA
            }

        }
    }

} # transformation

cor_list <- list(avg_corr = avg_corr,
                 avg_intersection_corr = avg_intersection_corr)

exp <- ifelse(e == 1 & c %in% "asv", "_exp", "")
saveRDS(cor_list, file = paste0(output_path, "cor_list", c,exp, "_", d, ".rds"))


# for (d in datasets$datasets) {
#     for (c in c("asv", "asv249", "asv199", "q2", "idtaxa", "q2_gg", "idtaxa_gg", "deftax", "KO", "cyc")) {
#         for (e in c(0,1)) {
    
#             exp <- ifelse(e == 1 & c %in% "asv", "_exp", "")
#             fname <- paste0(output_path, "cor_list", c,exp, "_", d, ".rds")

#             if (file.exists(fname)) {
#                 cors <- readRDS(fname) 
#                 cors[["avg_part_corr"]] <- NULL
#                 names(cors)[names(cors) %in% "avg_union_corr"] <- "avg_intersection_corr"

#                 saveRDS(cors, fname)
                    
#             }

#         }
    
#     }
# }
