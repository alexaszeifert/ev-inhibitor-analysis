library(tidyverse)
library(openxlsx)
library(pheatmap)
library(grid)

# load miRNA RPMs and standardize sample names
df <- read.xlsx("PATH/TO/mirna_summary_RPMs.xlsx")
df$miRNA <- sub("^hsa-", "", df$miRNA)
colnames(df) <- c('miRNA', 'spEV + Calpeptin', 'spEV + CytochalasinD', 'spEV (dH2O)',
                    'spEV (DMSO)','spEV + GW4869',  'spEV + Nexinhib20', 
                    'spEV + Pantethine', 'spEV + R5421', 'spEV + Y27632')
rpms <- df # keep original variable use downstream

# define set order and sample groups for annotations
abundance_cols <- c('spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421', 'spEV + GW4869',
                    'spEV + Nexinhib20', 'spEV + CytochalasinD', 
                    'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632')

control_columns <- c('spEV (DMSO)', 'spEV (dH2O)')
case_columns <- setdiff(abundance_cols, control_columns)
sample_groups <- data.frame(
  Group = ifelse(
    abundance_cols %in% control_columns, 
    "Control", 
    "Case"
  ), row.names = abundance_cols
)

# ------------------------------------------------------------------------------
# Figure 11B: top 100 by total RPM (log2+1) ------------------------------------
rpms$total_abundance <- rowSums(rpms[abundance_cols], na.rm = TRUE)
top100_abundant <- head(rpms$miRNA[order(rpms$total_abundance, decreasing = TRUE)], 100)
mat <- as.matrix(rpms[match(top100_abundant, rpms$miRNA), abundance_cols])
rownames(mat) <- top100_abundant
mat_log2 <- log2(mat + 1)

gt <- pheatmap(
  mat_log2,
  scale = "row",            
  cluster_rows = TRUE,
  cluster_cols = FALSE,      
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 100 most abundant miRNAs across treatments (log2 RPM)",
  annotation_col = sample_groups,
  annotation_colors = list(Group = c(Case = "#20A387FF", Control = "#333333"))
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)
ggsave("../figs/mirna_heatmap_top100_abundance_log2rpm.pdf", plot=gt, width = 10, height = 10)

# ------------------------------------------------------------------------------
# vehicle-adjusted deltas: build log2(treat) − log2(vehicle) matrix once
mat <- as.matrix(df[abundance_cols])
rownames(mat) <- make.unique(df$miRNA)
mat_log2  <- log2(mat + 1)

control_of <- c(                             # map each treatment to its vehicle
  'spEV + Calpeptin'     = 'spEV (DMSO)',
  'spEV + R5421'         = 'spEV (DMSO)',
  'spEV + GW4869'        = 'spEV (DMSO)',
  'spEV + Nexinhib20'    = 'spEV (DMSO)',
  'spEV + CytochalasinD' = 'spEV (DMSO)', 
  'spEV + Pantethine'    = 'spEV (dH2O)',
  'spEV + Y27632'        = 'spEV (dH2O)'
)
treat_cols <- names(control_of)

mat_log2_adj <- sapply(treat_cols, function(col) {
  mat_log2[, col] - mat_log2[, control_of[[col]]]
})
mat_log2_adj <- as.matrix(mat_log2_adj)
rownames(mat_log2_adj) <- rownames(mat_log2)

# column annotations/colors reused for all delta heatmaps ----------------------
annot_all <- data.frame(
  Vehicle = c(
    'spEV (DMSO)'          = "DMSO",
    'spEV + Calpeptin'     = "DMSO",
    'spEV + R5421'         = "DMSO",
    'spEV + GW4869'        = "DMSO",
    'spEV + Nexinhib20'    = "DMSO",
    'spEV + CytochalasinD' = "DMSO",
    'spEV (dH2O)'          = "dH2O",
    'spEV + Pantethine'    = "dH2O",
    'spEV + Y27632'        = "dH2O"
  )[colnames(mat_log2)],
  row.names = colnames(mat_log2)
)
annot_all$Vehicle <- factor(annot_all$Vehicle, levels = c("DMSO","dH2O"))
ann_colors <- list(
  Vehicle = c(DMSO = "#F4A6B5", 
              dH2O = "#B5A6F4") 
)

# delta heatmap 1: top 100 by mean absolute delta ------------------------------
mean_adj <- rowMeans(abs(mat_log2_adj), na.rm = TRUE)
top_mean_ids <- names(sort(mean_adj, decreasing = TRUE))[1:100]
m_delta_mean <- mat_log2_adj[top_mean_ids, treat_cols, drop = FALSE]
lim <- max(abs(range(m_delta_mean, finite = TRUE)))
col_delta <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(255)
breaks_delta <- seq(-lim, lim, length.out = 256)

gt <- pheatmap(
  m_delta_mean,
  scale = "none",               
  color = col_delta,
  breaks = breaks_delta,
  cluster_rows = TRUE,
  cluster_cols = FALSE,     
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 100 miRNAs with highest \nmean treatment–vehicle difference (log2 scale)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mirna_heatmap_top100_mean_diff_log2rpm.pdf", plot=gt, width = 10, height = 10)

# Figure 11C: top 100 by SD of delta
sds_adj  <- apply(mat_log2_adj, 1, sd,  na.rm = TRUE)
top_ids  <- names(sort(sds_adj, decreasing = TRUE))[1:100]
m_delta  <- mat_log2_adj[top_ids, treat_cols, drop = FALSE]
lim      <- max(abs(range(m_delta, finite = TRUE)))
col_delta <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(255)
breaks_delta <- seq(-lim, lim, length.out = 256)

gt <- pheatmap(
  m_delta,
  scale = "none",               
  color = col_delta,
  breaks = breaks_delta,
  cluster_rows = TRUE,
  cluster_cols = FALSE,     
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 100 miRNAs with highest variability in \ntreatment–vehicle differences (log2 scale, SD)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mirna_heatmap_top100_sd_diff_log2rpm.pdf", plot=gt, width = 10, height = 10)

# delta heatmap 3: top 100 by variance of delta

var_adj <- apply(mat_log2_adj, 1, var, na.rm = TRUE)
top_ids <- names(sort(var_adj, decreasing = TRUE))[1:100]
m_delta <- mat_log2_adj[top_ids, treat_cols, drop = FALSE]
lim <- max(abs(range(m_delta, finite = TRUE)))
col_delta <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(255)
breaks_delta <- seq(-lim, lim, length.out = 256)

gt <- pheatmap(
  m_delta,
  scale = "none",               
  color = col_delta,
  breaks = breaks_delta,
  cluster_rows = TRUE,
  cluster_cols = FALSE,     
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 100 miRNAs with highest variability in \ntreatment–vehicle differences (log2 scale, variance)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mirna_heatmap_top100_var_diff_log2rpm.pdf", plot=gt, width = 10, height = 10)

