.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(mia)
library(optparse)
library(tidyverse)
library(tibble)
library(ggpubr)
library(showtext)
library(vegan)

option_list <- list(
  make_option(c("-c", "--classifier"), type = "character", default = NULL,
              help = "Classifier", metavar = "character"),
  make_option(c("-d", "--dataset"), type = "character", default = NULL,
              help = "Dataset name", metavar = "character"),
  make_option(c("-l", "--level", type = "character", default = NULL,
              help = "Taxonomic level: asv or genus", metavar = "character")))

opt <- parse_args(OptionParser(option_list = option_list))

# Check required arguments
if (is.null(opt$dataset) || is.null(opt$level)) {
  stop("Both --dataset (-d) and --level (-l) must be provided")
}

if (opt$level %in% "genus" && is.null(opt$classifier)) {
  stop("--classifier (-c) must be provided when --level is set to genus")
}

if (!opt$level %in% c("asv", "asv249", "asv199", "asv_exp", "genus", "exp", "KO", "cyc", "exp_KO","exp_cyc")) {
  stop("--level should be set to asv, genus, or exp")
}

c <- opt$classifier
d <- opt$dataset
lev <- opt$level

# Start MDS + Permanova
output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"
plot_path <- "/projects/0/gs/danique/pipelines/int_results/figures/"

# For plot order and colors 
pipelines <- c("NGTax2", "eDNAFlow", "MetaWorks", "VSEARCH", "APSCALE", "QIIME2Deblur", "AMPtk", "LotuS2", "Natrix", "DADAist2", "dadasnake", 
        "nfcore_ampliseq", "QIIME2", "Tourmaline", "DADA2", "LotuS2_otu")

#### 
set.seed(1234)

# ASV-level 
if (lev %in% c("asv", "asv249", "asv199")) {

  print("ASV-level")

  if (lev %in% "asv") {
    tse <- readRDS(file = paste0(output_path, "tse_idtaxa_", d, "_asv_allpipelines.rds"))
  } else {
    trim <- gsub("asv", "", lev)
    tse <- readRDS(file = paste0(output_path, "tse_q2_", d, "_trim", trim, "_asv_allpipelines.rds")) 
  }
  
  # assay(tse, "counts")[is.na(assay(tse, "counts"))] <- 0
  tse <- transformAssay(tse, assay.type = "counts", method = "clr", pseudocount = TRUE)
  tse <- transformAssay(tse, assay.type = "counts", method = "pa", pseudocount = TRUE)

  if (lev %in% "asv") {
    tse <- tse[, !tse$pipeline %in% c("NGTax2", "AMPtk", "LotuS2", "LotuS2_otu", "QIIME2Deblur")]
  } else if (lev %in% "asv249") {
    tse <- tse[, !tse$pipeline %in% c("NGTax2", "LotuS2", "LotuS2_otu")]
  } else {
    tse <- tse[, !tse$pipeline %in% c("NGTax2")]
  }
  
  tse$sample <- str_remove(tse$sample, tse$pipeline)
  tse$sample <- str_remove(tse$sample, "_")

  meta <- data.frame(colData(tse))

  # Calculate distances
  print("Get distances")

  for (mydist in c("aitchison", "jaccard")) {

    if (mydist %in% "aitchison") {
      atype <- "clr"
      method = "euclidean"
    } else {
      atype = "counts"
      method = "jaccard"
    }

    mat <- getDissimilarity(tse, method = method, assay.type = atype) %>% as.matrix()

    m <- lapply(unique(meta$sample), function(s) {
      mat_tmp <- mat[grepl(s, rownames(mat)), grepl(s, colnames(mat))]
      diag(mat_tmp) <- NA

      m <- data.frame(t(rowMeans(mat_tmp, na.rm = TRUE))) 
      colnames(m) <- str_remove_all(colnames(m), paste0("_", s))
      return(m)
    }) %>% bind_rows()

    m$sample <- unique(meta$sample)

    # Run MDS on distance matrix
    mat <- as.dist(mat)
    mds_res <- cmdscale(mat, eig = TRUE) 
    
    eig_vec <- mds_res$points %>% data.frame()
    colnames(eig_vec) <- c("PC1", "PC2")

    eig_vec <- eig_vec %>%
      rownames_to_column(var = "ID") %>%
      dplyr::left_join(., meta, by = "ID") %>% 
      column_to_rownames(var = "ID")

    # Get eigen values
    eig_vals <- mds_res$eig
    positive_eig <- eig_vals[eig_vals > 0]

    # Variance explained
    varexp_1 <- positive_eig[1] / sum(positive_eig)
    var_exp2 <- positive_eig[2] / sum(positive_eig)

    eig_vec$pipeline <- factor(eig_vec$pipeline, levels = pipelines)

    m_save <- list(m = m,
                   mat = mat,
                   mds = mds_res,
                   eig_vec = eig_vec)

    saveRDS(m_save, file = paste0(output_path, "dist_mds_", mydist, "_", d, "_", lev, ".rds"))

    w <- abs(diff(c(min(eig_vec$PC1), max(eig_vec$PC1)))) * 0.02
    h <- abs(diff(c(min(eig_vec$PC2), max(eig_vec$PC2)))) * 0.02


  }
}

if (lev %in% c("asv_exp")) {

  # load tse
  tse <- readRDS(file = paste0(output_path, "tse_exp_", d, "_asv_allpipelines.rds"))
  
  # transform
  tse <- transformAssay(tse, assay.type = "relabundance", method = "clr", pseudocount = TRUE)

  tse$sample <- str_remove(tse$sample, tse$pipeline)
  tse$sample <- str_remove(tse$sample, "_")

  meta <- data.frame(colData(tse))

  # Calculate distances
  print("Get distances")

  for (mydist in c("aitchison", "jaccard")) {

    if (mydist %in% "aitchison") {
      atype <- "clr"
      method = "euclidean"
    } else {
      atype = "relabundance"
      method = "jaccard"
    }

    mat <- getDissimilarity(tse, method = method, assay.type = atype) %>% as.matrix()

    mat_exp <- data.frame(mat) %>% 
      select(all_of(contains("Expected_"))) %>%
      filter(!grepl("Expected_", rownames(.)))

    m_save <- list(mat = mat,
                   mat_exp = mat_exp)

    saveRDS(m_save, file = paste0(output_path, "dist_mds_", mydist, "_", d, "_", lev, ".rds"))

  }
}


# Genus-level
if (lev %in% c("genus", "KO", "cyc")) {
  
  print(paste0(lev, " level"))
  
  if (lev %in% "genus") {
    tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_genus_allpipelines.rds"))
  } else {
     tse <- readRDS(file = paste0(output_path, "tse_", lev, "_", d, "picrust_allpipelines.rds"))
  }
  
  # tse <- tse[, !tse$pipeline %in% c("NGTax2")] # "Expected"
  tse <- transformAssay(tse, assay.type = "counts", method = "clr", pseudocount = TRUE)
  tse <- transformAssay(tse, assay.type = "counts", method = "pa", pseudocount = TRUE)
  meta <- data.frame(colData(tse))

  # Get distances
  for (mydist in c("jaccard", "aitchison")) {

    print("Get distances")

    if (mydist %in% "aitchison") {
      atype="clr"
      method="euclidean"
    } else {
      atype="counts"
      method="jaccard"
    }
        
    mat <- getDissimilarity(tse, method = method, assay.type = atype) %>% as.matrix()

    m <- lapply(unique(meta$sample), function(s) {
      mat_tmp <- mat[grepl(s, rownames(mat)), grepl(s, colnames(mat))]
      diag(mat_tmp) <- NA

      m <- data.frame(t(rowMeans(mat_tmp, na.rm = TRUE))) 
      colnames(m) <- str_remove_all(colnames(m), paste0("_", s))
      return(m)
    }) %>% bind_rows()

    m$sample <- unique(meta$sample)

    # Run MDS on distance matrix
    mat <- as.dist(mat)
    mds_res <- cmdscale(mat, eig = TRUE) 

    eig_vec <- mds_res$points %>% data.frame()
    colnames(eig_vec) <- c("PC1", "PC2")

    eig_vec <- eig_vec %>%
      rownames_to_column(var = "ID") %>%
      dplyr::left_join(., meta, by = "ID") %>% 
      column_to_rownames(var = "ID")

    # Get eigen values
    eig_vals <- mds_res$eig
    positive_eig <- eig_vals[eig_vals > 0]

    # Variance explained
    varexp_1 <- positive_eig[1] / sum(positive_eig)
    var_exp2 <- positive_eig[2] / sum(positive_eig)

    eig_vec$pipeline <- factor(eig_vec$pipeline, levels = pipelines)

    # Contributors to mds axes
    if (lev %in% "genus") {
      sample_order <- rownames(mds_res$points)
      if (mydist %in% "jaccard") atype <- "pa"
      X <- assay(tse, assay = atype)[, sample_order, drop = FALSE]
      fit <- envfit(mds_res$points, t(X))
      contrib <- fit$vectors$arrows

      m_save <- list(m = m,
                    mat = mat,
                    mds = mds_res,
                    eig_vec = eig_vec,
                    contrib = contrib)
    } else {
      m_save <- list(m = m,
                    mat = mat,
                    mds = mds_res,
                    eig_vec = eig_vec)
    }

    saveRDS(m_save, file = paste0(output_path, "dist_mds_", c, "_", mydist, "_", d,  "_", lev, ".rds"))

    w <- abs(diff(c(min(eig_vec$PC1), max(eig_vec$PC1)))) * 0.02
    h <- abs(diff(c(min(eig_vec$PC2), max(eig_vec$PC2)))) * 0.02

    # Permanova + leave-one-out
    res <- list()
    print("Permanova with all pipelines")
    
    b <- ifelse(mydist %in% "jaccard", TRUE, FALSE)

    res[[mydist]][["all"]] <- getPERMANOVA(
            tse,
            assay.type = atype,
            formula = x ~ pipeline,
            method = method,
            binary = b,
            permutations = 999,
            strata = factor(tse$sample))

    for (p in unique(tse$pipeline)) {
        tse_tmp <- tse[, !tse$pipeline %in% p]
        tse_tmp$pipeline <- factor(tse_tmp$pipeline)

        print(paste("Excluding pipeline: ", p))	

        res[[mydist]][[p]] <- getPERMANOVA(
            tse_tmp,
            assay.type = atype,
            formula = x ~ pipeline,
            method = method,
            binary = b,
            permutations = 999,
            strata = factor(tse_tmp$sample))
    }

    saveRDS(res, file = paste0(output_path, "permanova_", lev, c, "_",mydist, "_", d, "_leaveoneout.rds"))
  }
}


# Expected level
if (lev %in% c("exp", "exp_KO", "exp_cyc")) {
  print("MDS with expected composition level")

  if (lev %in% c("exp_KO", "exp_cyc")) {
    c <- gsub("exp_", "", lev)
    lev <- "exp"
    tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_picrust_exp_allpipelines.rds")) 
    tse <- transformAssay(tse, assay.type = "counts", method = "relab", pseudocount = FALSE)
  } else if (c %in% c("q2", "idtaxa")) {
    tse <- readRDS(file = paste0(output_path, "tse_", c, "_", d, "_genus_exp_allpipelines.rds"))  
    tse <- transformAssay(tse, assay.type = "relabundance", method = "relab", pseudocount = FALSE)
  }

  meta <- data.frame(colData(tse))

  pipelines <- c("Expected", pipelines)
  cols_p <- c("Expected" = "black", cols_p)

  # Get distances
  for (mydist in c("aitchison", "jaccard")) {
    print("Get distances")

    if (mydist %in% "aitchison") {
      tse <- transformAssay(tse, assay.type = "relab", method = "clr", pseudocount = TRUE)
      atype="clr"
      method="euclidean"
    } else {
      tse <- transformAssay(tse, assay.type = "relab", method = "pa", pseudocount = FALSE)
      atype="pa"
      method="jaccard"
    }

    mat <- getDissimilarity(tse, method = method, assay.type = atype) %>% as.matrix()
    mat2 <- mat
    mat <- as.data.frame(mat) %>% 
      dplyr::select(all_of(contains("Expected"))) %>%
      filter(!grepl("Expected", rownames(.)))

    m <- data.frame(pipeline = pipelines)
    
    for (s in unique(meta$sample)) {
      mat_tmp <- mat %>% 
        filter(grepl(s, rownames(.))) %>%
        dplyr::select(all_of(contains(s))) %>%
        mutate(pipeline = str_remove(rownames(.), paste0("_", s)))

      m <- m %>% left_join(., mat_tmp, by = "pipeline")
    }

    colnames(m) <- str_remove(colnames(m), paste0("Expected", "_"))
    
    # Run MDS on distance matrix
    mat <- mat2
    mat <- as.dist(mat)
    mds_res <- cmdscale(mat, eig = TRUE)

    eig_vec <- mds_res$points %>% data.frame()
    colnames(eig_vec) <- c("PC1", "PC2")

    eig_vec <- eig_vec %>%
      rownames_to_column(var = "ID") %>%
      dplyr::left_join(., meta, by = "ID") %>% 
      column_to_rownames(var = "ID")

    # Get eigen values
    eig_vals <- mds_res$eig
    positive_eig <- eig_vals[eig_vals > 0]

    # Variance explained
    varexp_1 <- positive_eig[1] / sum(positive_eig)
    var_exp2 <- positive_eig[2] / sum(positive_eig)

    eig_vec$pipeline <- factor(eig_vec$pipeline, levels = pipelines)

    w <- abs(diff(c(min(eig_vec$PC1), max(eig_vec$PC1)))) * 0.02
    h <- abs(diff(c(min(eig_vec$PC2), max(eig_vec$PC2)))) * 0.02

    meta_eig <- eig_vec %>%
      rownames_to_column(var = "ID") %>%
      mutate(alpha = factor(ifelse(pipeline %in% "Expected", "Expected", "Observed")))

    m_save <- list(m = m,
                   mat = mat,
                   mds = mds_res,
                   eig_vec = eig_vec,
                   meta_eig = meta_eig)

    saveRDS(m_save, file = paste0(output_path, "dist_mds_", c, "_", mydist, "_", d,  "_exp.rds"))

        
  }
}
