
# Anonymize sample names in `alloutput`
# Output used for public release, so that no real sample IDs are present in the data 

# Every sample gets a name "<dataset><n>" 
# One mapping is built per dataset (common mapping for simulations) and applied to every pipeline, so a sample has the same anonymized name in reads_in, all feature tables, raw files, etc.
.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(dplyr)
library(tibble)
library(readr)

output_path <- "/projects/0/gs/danique/pipelines/int_results/routput/"

# Collect real sample ids 
# Primary source: reads_in$sample_names
# Extra IDs: column names of any pipeline's feature table not already covered.
collect_sample_ids <- function(x) {
  ids <- character(0)
  if (!is.null(x[["reads_in"]])) ids <- as.character(x[["reads_in"]][["sample_names"]])

  for (el in x) {
    if (is.list(el) && !is.data.frame(el) && !is.null(el[["feat"]])) {
      extra <- colnames(el[["feat"]])
      # skip names that are just make.names()-mangled versions of known IDs
      extra <- extra[!extra %in% ids & !extra %in% make.names(ids)]
      ids <- c(ids, extra)
    }
  }
  ids <- unique(ids)
  ids[!is.na(ids) & ids != ""]
}

# Build mapping for one dataset
make_mapping <- function(ids, prefix) {
  anon <- paste0(prefix, seq_along(ids))
  map  <- setNames(anon, ids)

  mangled <- make.names(ids)
  add <- mangled != ids & !mangled %in% ids
  map <- c(map, setNames(anon[add], mangled[add]))

  list(map = map,
       key = data.frame(original = ids, anonymized = anon, stringsAsFactors = FALSE))
}

# Replacement helpers 
# Exact matches (column names, row names, cell values)
recode <- function(v, map) {
  if (is.null(v)) return(v)
  hit <- v %in% names(map)
  v[hit] <- map[v[hit]]
  v
}

# IDs embedded in longer strings (raw fastq file names).
# Longest IDs first (so "S12" is handled before "S1"), via placeholder tokens so
# an anonymized name can never be matched again by a later ID.
recode_substr <- function(v, map) {
  keys <- names(map)[order(nchar(names(map)), decreasing = TRUE)]
  tok  <- paste0("\u0001", seq_along(keys), "\u0002")
  for (i in seq_along(keys)) v <- gsub(keys[i], tok[i], v, fixed = TRUE)
  for (i in seq_along(keys)) v <- gsub(tok[i], map[[keys[i]]], v, fixed = TRUE)
  v
}

# Walk the nested list. Data frames: names, row names, character/factor columns.

anonymize <- function(x, map) {
  if (is.data.frame(x)) {
    names(x)    <- recode(names(x), map)
    rownames(x) <- recode(rownames(x), map)
    for (col in seq_along(x)) {
      if (is.character(x[[col]])) x[[col]] <- recode(x[[col]], map)
      if (is.factor(x[[col]]))    levels(x[[col]]) <- recode(levels(x[[col]]), map)
    }
  } else if (is.matrix(x)) {
    if (!is.null(dimnames(x))) dimnames(x) <- lapply(dimnames(x), recode, map = map)
    if (is.character(x)) x[] <- recode(x, map)
  } else if (is.list(x) && !isS4(x)) {
    if (!is.null(names(x))) names(x) <- recode(names(x), map)
    for (i in seq_along(x)) {
      if (!is.null(x[[i]])) x[[i]] <- anonymize(x[[i]], map)
    }
  } else if (is.character(x)) {
    x <- recode(x, map)
  } else if (is.factor(x)) {
    levels(x) <- recode(levels(x), map)
  }
  x
}

# Apply per dataset 
alloutput <- readRDS(paste0(output_path, "all_tables_combined.rds"))

alloutput_anon <- alloutput
anon <- list()
anon_keys <- list()

for (d in c("gemsim", "mindset")) {
  ids <- collect_sample_ids(alloutput[[d]])
  m <- make_mapping(ids, d)
  anon[[d]] <- m
}

ids <- collect_sample_ids(alloutput[["miseq_v4_sim"]])
m <- make_mapping(ids, "sim")
anon[["miseq_v4_sim"]] <- anon[["novaseq_v4_sim"]]  <- m

for (d in names(alloutput)) {
    
  if (grepl("mock", d)) {
    alloutput_anon[[d]][["raw"]] <- NULL      
    message("Skipping ", d)
    next
  }

  m <- anon[[d]]  
  anon_keys[[d]] <- m$key
  x <- anonymize(alloutput[[d]], m$map)
  x[["raw"]] <- NULL # remove fastq paths

  alloutput_anon[[d]] <- x
  message("Anonymized ", d)

}

# Check for leftover real IDs
all_strings <- function(x) {
  out <- character(0)
  if (isS4(x)) return(out)
  if (is.data.frame(x) || is.matrix(x)) {
    out <- c(out, rownames(x), colnames(x))
    if (is.data.frame(x)) {
      for (col in seq_along(x)) {
        if (is.character(x[[col]])) out <- c(out, x[[col]])
        if (is.factor(x[[col]]))    out <- c(out, levels(x[[col]]))
      }
    } else if (is.character(x)) out <- c(out, as.vector(x))
  } else if (is.list(x)) {
    out <- c(out, names(x))
    for (el in x) out <- c(out, all_strings(el))
  } else if (is.character(x)) {
    out <- c(out, x)
  } else if (is.factor(x)) {
    out <- c(out, levels(x))
  }
  unique(out[!is.na(out)])
}

for (d in names(alloutput_anon)) {
  print(d)
  s <- all_strings(alloutput_anon[[d]])
  orig <- anon_keys[[d]]$original
  
  exact <- orig[orig %in% s | make.names(orig) %in% s]
  
  # substring check only for IDs long enough to be meaningful
  long <- orig[nchar(orig) >= 4]
  sub  <- long[vapply(long, function(id) any(grepl(id, s, fixed = TRUE)), logical(1))]
  left <- unique(c(exact, sub))
  if (length(left) > 0) {
    warning("Dataset ", d, ": ", length(left), " real ID(s) still present, e.g. ",
            paste(head(left, 5), collapse = ", "))
  }
}

# Save
# saveRDS(alloutput_anon, file = paste0(output_path, "all_tables_combined_anonymized.rds"))
pipeline_output <- alloutput_anon
save(pipeline_output, file = paste0(output_path, "pipeline_output.RData"))

# Keep the key PRIVATE: store it outside anything you publish.
key_dir <- paste0(output_path, "anonymization_keys/")
dir.create(key_dir, showWarnings = FALSE, recursive = TRUE)

for (d in names(anon_keys)) {
  write_tsv(anon_keys[[d]], paste0(key_dir, d, "_sample_key.tsv"))
}

saveRDS(anon_keys, paste0(key_dir, "all_sample_keys.rds"))

##############################################################################################
# Anonymize sample names in simulation files
# Output used for public release, so that no real sample IDs are present in the data 

# Count table 
sim_dir <- "/projects/0/gs/danique/pipelines/simulations/input"
counts_file=file.path(sim_dir, "implanted_counts.tsv") 
counts <- read_tsv(counts_file)

# Mapping
key_dir <- paste0(output_path, "anonymization_keys/")
anon_keys <- readRDS(paste0(key_dir, "all_sample_keys.rds"))
mapping <- anon_keys[["miseq_v4_sim"]]

df <- data.frame(old = colnames(counts)) %>%
  left_join(., mapping, by = c("old" = "original"))

df$anonymized[1] <- "featureID"
colnames(counts) <- df$anonymized

write_tsv(counts, file.path(sim_dir, "implanted_counts_anon.tsv"))


# Fastq files
dir.create(file.path(sim_dir, "fastq_anon"))

for (d in c("miseq", "novaseq")) {
  
  anon_files <- file.path(sim_dir, "fastq_anon", d)
  dir.create(anon_files)

  fastq_dir <- file.path(sim_dir, "fastq", d)
  fastq_files <- list.files(fastq_dir, full.names = TRUE)

  fn  <- basename(fastq_files)
  ids <- sub("_L\\d{3}_R[12]_\\d{3}\\.fastq\\.gz$", "", fn)   # "s23042600364"
  idx <- match(ids, mapping$original)

  if (anyNA(idx)) {
    warning(d, ": no mapping for ", paste(unique(ids[is.na(idx)]), collapse = ", "))
  }
  
  keep <- !is.na(idx)

  new_fn <- paste0(mapping$anonymized[idx[keep]],
                   substring(fn[keep], nchar(ids[keep]) + 1))   # keep "_L001_R2_001.fastq.gz"
  
  new_filenames <- file.path(anon_files, new_fn)

  file.copy(fastq_files[keep], new_filenames, overwrite = TRUE)

}