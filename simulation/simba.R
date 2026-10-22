# Install simbaR (only when changes were made)
# devtools::install("/projects/0/gs/danique/pipelines/simulations/SIMBA")

# Load packages
library("tidyverse") 
library("simbaR")
library("compositions")
library("optparse")
library("rhdf5")

option_list <- list(
  make_option(c("--count"), type = "character", help = "Path to count table file", metavar = "FILE"),
  make_option(c("--metadata"), type = "character", help = "Path to metadata file", metavar = "FILE"),
  make_option(c("--taxonomy"), type = "character", help = "Path to taxonomy file", metavar = "FILE"),
  make_option(c("--abundance"), type = "character", help = "Abundance scaling, comma-separated string, e.g. '1.5,2'"),
  make_option(c("--prevalence"), type = "double", help = "Prevalence scaling, e.g. 0.1"),
  make_option(c("--feature_type"), type = "character", help = "Feature type, e.g. 'middle_high'"),
  make_option(c("--output_count"), type = "character", help = "Output implanted counts file", metavar = "FILE"),
  make_option(c("--list"), type = "character", help = "Output overview list (.rds)", metavar = "FILE"),
  make_option(c("--tmp"), type = "character", help = "Temporary directory path", metavar = "DIR"),
  make_option(c("--keep"), type = "logical", help = "Whether to keep temporary files (TRUE/FALSE)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#########################################################################################
#                                     Importing data                                    #
#########################################################################################

feat <- opt$count %>% read.table(header = TRUE, row.names = 1, sep = "\t", comment.char = "")
meta <- opt$metadata %>% read.table(header = TRUE, row.names = 1, sep = "\t", comment.char = "")
tax <- opt$taxonomy %>% read.table(header = TRUE, row.names = 1, sep = "\t", comment.char = "")

tax <- tax[match(rownames(feat), rownames(tax)), , drop = FALSE]
identical(rownames(feat), rownames(tax))

t <- data.frame(genus = sapply(tax$Taxon, function(x) str_split(x, "; ")[[1]][6]) %>% unlist())
t$genus <- str_remove_all(t$genus, "g__")
t$genus <- str_replace_all(t$genus, "_", ".")
t$fulltax <- tax$Taxon
rownames(t) <- rownames(tax)
t$featnames <- paste0(t$genus, "_", rownames(t))
rownames(feat) <- t$featnames

meta$Individual_ID <- rownames(meta)
meta$Timepoint <- 1

feat[1:3, 1:3]
meta[1:2, 1:2]
tax[1:3, 1:1]
t[1:3, 1:3]

# matching sample names
identical(colnames(feat),rownames(meta))
meta <- meta[match(colnames(feat), rownames(meta)),]
identical(colnames(feat),rownames(meta)) 

feat_b <- feat

#########################################################################################
#                                      data checks                                      #
#########################################################################################

feat_t <- as.data.frame(t(feat_b))
genera <- sapply(strsplit(colnames(feat_t), "_"), `[`, 1)
gen_t <- t(rowsum(t(feat_t), genera, reorder = FALSE))

prev <- data.frame(p = colSums(gen_t > 0) / nrow(gen_t))

prev_b <- prev

#########################################################################################
#                           Setting parameters for simulation                           #
#########################################################################################

# filtering parameters
ab.cutoff <- 0
prev.cutoff <- 0
log.n0 <- 1e-05 # for reality check

# other parameters for the simulations
prop.markers <- c(0.1) # 5% and 10% of genera
class.balance <- 0.5
subset_size <- nrow(meta)

# effect sizes
ab.scale <- opt$abundance %>% strsplit(",") %>% unlist() %>% as.numeric()
prev.scale <- opt$prevalence %>% strsplit(",") %>% unlist() %>% as.numeric()
feature.type <- opt$feature_type 
feature.level <- "genus" # asv or genus

#########################################################################################
#                                         Implantation                                  #
#########################################################################################

# Remove the simulation file if it exists
remove_sim <- function(sim.name) {
  if (file.exists(sim.implant)) {
    file.remove(sim.implant)
    message("Removed existing simulation file: ", sim.implant)
  }
}

sim.implant <- paste0(opt$tmp, '/sim.h5')
remove_sim(sim.implant)

set.seed(1234) # for reproducibility

create.data.simulation(feat = feat,
                meta = meta,
                sim.location = sim.implant,
                sim.type ='cross-section',
                sim.method ='resampling',
                filt.params = list(ab.cutoff=as.numeric(ab.cutoff),
                                   prev.cutoff=prev.cutoff,
                                   log.n0=as.numeric(log.n0)),
                sim.params = list(ab.scale=ab.scale,
                                  prev.scale=prev.scale,
                                  prop.markers=prop.markers,
                                  conf=NULL,
                                  class.balance=class.balance,
                                  feature.type=feature.type,
                                  feature.level=feature.level,
                                  repeats=10,
                                  balanced=TRUE))


#######################################################################
#######################################################################
#######################################################################

# Load the simulation results
temp <- rhdf5::h5ls(sim.implant, recursive=FALSE)
simulation_names <- setdiff(temp$name, c('original_data', 'simulation_parameters'))

all_simulations <- lapply(simulation_names, function(x) {
  h5read(sim.implant, name = x)
})

names(all_simulations) <- simulation_names

# Assess whether we can detect the implanted effects
df <- data.frame()
results_gen <- list()

for (x in simulation_names) {
  sim <- all_simulations[[x]]
  meta <- data.frame(labels = sim$labels, row.names = sim$sample_names)
  
  feat <- t(sim$features) 
  colnames(feat) <- sim$feature_names
  feat_clr <- clr(feat + 0.5) # clr transformation
  
  # Genus-level analysis
  genera <- sapply(strsplit(colnames(feat), "_"), `[`, 1)
  gen <- t(rowsum(t(feat), genera, reorder = FALSE))
  gen_clr <- clr(gen + 0.5) # clr transformation for genus-level

  prev <- apply(gen, 2, function(x) sum(x > 0) / nrow(gen_clr)) > 0.1 
  gen_clr_prev <- gen_clr[, prev]  

  mods_gen <- lapply(colnames(gen_clr_prev), function(g) lm(gen_clr_prev[,g] ~ meta$labels))
  results_gen[[x]] <- do.call(rbind, lapply(mods_gen, function(m) {
    coefs <- summary(m)$coefficients
    if ("meta$labels" %in% rownames(coefs)) coefs["meta$labels",] else rep(NA,4)
  }))
  results_gen[[x]] <- data.frame(results_gen[[x]], padj = p.adjust(results_gen[[x]][,4], "fdr"))
  rownames(results_gen[[x]]) <- colnames(gen_clr_prev)
  
  # Calculate metrics
  calc_metrics <- function(results, true_markers) {
    sig <- rownames(results)[results$padj < 0.05]
    TP <- sum(sig %in% true_markers)
    FP <- sum(!sig %in% true_markers)
    FN <- sum(!true_markers %in% sig)
    TN <- nrow(results) - TP - FP - FN
    data.frame(TPR = TP/(TP + FN), FPR = FP/(FP + TN))
  }

  # Genus-level metrics (using genus-matched markers)
  true_markers_gen <- genera[sapply(genera, function(g) any(grepl(g, paste0(sim$marker_idx, "_*"), fixed = TRUE)))] %>% unique()
  df <- rbind(df, data.frame(
    simulation = x,
    level = "genus",
    calc_metrics(results_gen[[x]], true_markers_gen)
  ))
}

df %>% filter(level %in% "genus") %>% arrange(desc(TPR))
saveRDS(sim, file = paste0(opt$tmp, "/implantation_lmstats.rds"))

# Select the simulation with the best performance 
bestsim <- df %>% filter(TPR == 1, FPR == 0) %>% slice(1) %>% pull(simulation)
sim <- all_simulations[[bestsim]]

markers <- sim$marker_idx
markers_gen <- sapply(strsplit(sim$marker_idx, "_"), `[`, 1) %>% unique()

sim_feat <- t(sim$features)

id <- sapply(strsplit(sim$feature_names, "_"), `[`, 2)
genera <- sapply(strsplit(sim$feature_names, "_"), `[`, 1)
sim_gen <- t(rowsum(t(sim_feat), genera, reorder = FALSE))
sim_gen_clr <- clr(sim_gen + 0.5)
sim_gen_rel <- sim_gen / rowSums(sim_gen) * 100

ave_abund <- data.frame(genus = colnames(sim_gen_rel), 
                        abund_rel = colMeans(sim_gen_rel, na.rm = TRUE),
                        abund_clr = colMeans(sim_gen_clr, na.rm = TRUE))

prev <- as.data.frame(colMeans(sim_gen > 0, na.rm = TRUE)) %>%
  rownames_to_column("genus")
colnames(prev) <- c("genus", "prev")

best_results <- ave_abund %>%
  left_join(prev, by = "genus") %>%
  left_join(rownames_to_column(results_gen[[bestsim]], "genus"), by = "genus") %>%
  mutate(implanted = ifelse(genus %in% markers_gen, "yes", "no"),
         detected = ifelse(implanted == "yes" & padj < 0.05, "yes", "no"),
         tp = case_when(
           prev < 0.1 ~ "not tested",
           implanted == "yes" & detected == "yes" ~ "implanted + detected",
           implanted == "yes" & detected == "no" ~ "implanted + not detected",
           implanted == "no" & detected == "yes" ~ "not implanted + detected",
           implanted == "no" & detected == "no" ~ "not implanted + not detected"))

write.table(best_results, file = paste0(opt$tmp, "/lm_stats_implantation.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)           
           
# Save the implanted counts
rownames(sim_feat) <- sim$sample_names
colnames(sim_feat) <- id

sim_feat_save <- as.data.frame(t(sim_feat)) %>%
  rownames_to_column("featureID") 

write.table(sim_feat_save, file = paste0(opt$tmp, "/implanted_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)    

# Write metafile
meta <- data.frame(sample_names = sim$sample_names, labels = ifelse(sim$labels == 1, 1, 0))
exp <- data.frame(marker_idx=sim$marker_idx, 
                  genus=sapply(strsplit(sim$marker_idx, "_"), `[`, 1),
                  ab_shift=1.5, 
                  prev_shift=0.2)

write.table(meta, file = paste0(opt$tmp, "/implanted_meta.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)              
write.table(exp, file = paste0(opt$tmp, "/expected_associations.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)  

# Plot the detection of implanted features
plt <- ggplot() +
  geom_point(data = subset(best_results, implanted %in% "no"),
             aes(x = prev, y = abund_clr, color = tp), size = 2, alpha = 0.9) +
  geom_point(data = subset(best_results, implanted %in% "yes"),
             aes(x = prev, y = abund_clr, color = tp), size = 2, alpha = 1) +
  scale_color_manual(values = c(
    "not tested" = "#eeecec", 
    "implanted + detected" = "#7BC17E",
    "implanted + not detected" = "#e59295",
    "not implanted + detected" = "#d3494e",
    "not implanted + not detected" = "#d4d4d4")) +
  theme_minimal() +
  labs(color = "Implanted", x = "Prevalence", y = "CLR Abundance")

ggsave(paste0(opt$tmp, "/abun_prev_implanted.pdf"), plot = plt, width = 8, height = 4)
