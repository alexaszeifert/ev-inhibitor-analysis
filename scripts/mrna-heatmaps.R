library(tidyverse)
library(pheatmap)
library(grid)
library(org.Hs.eg.db)
library(AnnotationDbi)

df <- read_csv('PATH/TO/rna_expression.csv', show_col_types = FALSE) %>%
  dplyr::select('rna_id', 'gene_id',
                'tpm_spEV_DMSO', 'tpm_spEV_50ug_ml_Calpeptin', 'tpm_spEV_50uM_R5421',
                'tpm_spEV_10uM_GW4869', 'tpm_spEV_10uM_Nexinhib20', 'tpm_spEV_1uM_CytochalasinD',
                'tpm_spEV_dH2O', 'tpm_spEV_10uM_Pantethine', 'tpm_spEV_10uM_Y27632') %>%
  filter(!grepl("^XM_|^rna-XM_", rna_id))

# define columns, labels, and sample group annotation (control vs case)
tpm_cols <- c('tpm_spEV_DMSO', 'tpm_spEV_50ug_ml_Calpeptin', 
              'tpm_spEV_50uM_R5421', 'tpm_spEV_10uM_GW4869',
              'tpm_spEV_10uM_Nexinhib20', 'tpm_spEV_1uM_CytochalasinD',
              'tpm_spEV_dH2O', 'tpm_spEV_10uM_Pantethine',
              'tpm_spEV_10uM_Y27632')
control_columns <- c('tpm_spEV_DMSO', 'tpm_spEV_dH2O')
sample_groups <- data.frame(
  Group = ifelse(tpm_cols %in% control_columns, "Control", "Case"),
  row.names = tpm_cols
)
label_map <- c(
  'tpm_spEV_DMSO'              = 'spEV (DMSO)',
  'tpm_spEV_50ug_ml_Calpeptin' = 'spEV + Calpeptin',
  'tpm_spEV_50uM_R5421'        = 'spEV + R5421',
  'tpm_spEV_10uM_GW4869'       = 'spEV + GW4869',
  'tpm_spEV_10uM_Nexinhib20'   = 'spEV + Nexinhib20',
  'tpm_spEV_1uM_CytochalasinD' = 'spEV + CytochalasinD',
  'tpm_spEV_dH2O'              = 'spEV (dH2O)',
  'tpm_spEV_10uM_Pantethine'   = 'spEV + Pantethine',
  'tpm_spEV_10uM_Y27632'       = 'spEV + Y27632'
)

# map ENTREZ: gene names and drop hemoglobin entries
entrez_ids <- as.character(df$gene_id)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
gene_names <- mapIds(org.Hs.eg.db,
                     keys = entrez_ids,
                     column = "GENENAME",
                     keytype = "ENTREZID",
                     multiVals = "first")
df$gene_symbol <- gene_symbols
df$gene_name <- gene_names
df$rna_id_clean <- sub("^rna-(NM_)", "\\1", df$rna_id)
df <- df[!grepl("hemoglobin", df$gene_name, ignore.case = TRUE), ]

# ------------------------------------------------------------------------------
# Figure 10B: top 100 by total TPM (row-z-scored)
df$total_tpm <- as.numeric(rowSums(df[tpm_cols], na.rm = TRUE))
top100 <- df %>% 
  dplyr::arrange(dplyr::desc(total_tpm)) %>% 
  dplyr::slice_head(n = 100)
mat <- as.matrix(top100[tpm_cols])
rownames(mat) <- make.unique(top100$rna_id_clean)
mat_log2 <- log2(mat + 1)

labs       <- unname(label_map[colnames(mat)])

gt <- pheatmap(
  mat_log2,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 100 most abundant mRNAs across treatments (log2 TPM)",
  annotation_col = sample_groups,
  annotation_colors = list(Group = c(Case = "#20A387FF", Control = "#333333"))
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)
ggsave("../figs/mrna_heatmap_top100_abundance_log2tpm.pdf", plot = gt, width = 10, height = 10)

# ------------------------------------------------------------------------------
# vehicle-adjusted deltas: build log2(treat) − log2(vehicle) matrix once
mat <- as.matrix(df[tpm_cols])
rownames(mat) <- make.unique(df$rna_id_clean)
mat_log2  <- log2(mat + 1)
colnames(mat_log2) <- unname(label_map[colnames(mat_log2)])

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
  Group   = c(Case = "#20A387", Control = "#333333"),
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
  main = "Top 100 mRNAs with highest \nmean treatment–vehicle difference (log2 scale)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mrna_heatmap_top100_mean_diff_log2tpm.pdf", plot = gt, width = 10, height = 10)

# Figure 10C: top 100 by SD of delta
sds_adj <- apply(mat_log2_adj, 1, sd, na.rm = TRUE)
top_ids <- names(sort(sds_adj, decreasing = TRUE))[1:100]
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
  main = "Top 100 mRNAs with highest variability in \ntreatment–vehicle differences (log2 scale, SD)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mrna_heatmap_top100_sd_diff_log2tpm.pdf", plot=gt, width = 10, height = 10)

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
  main = "Top 100 mRNAs with highest variability in \ntreatment–vehicle differences (log2 scale, variance)",
  annotation_col    = annot_all[treat_cols, , drop = FALSE],
  annotation_colors = ann_colors,
  # visually separate DMSO block from dH2O block
  gaps_col = 5, 
  width = 10, height = 10
)$gtable

grid::grid.newpage()
grid::grid.draw(gt)

ggsave("../figs/mrna_heatmap_top100_var_diff_log2tpm.pdf", plot=gt, width = 10, height = 10)
