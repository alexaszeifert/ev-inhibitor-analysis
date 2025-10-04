library(tidyverse)
library(readxl)
library(pheatmap)
library(grid)

# load & clean protein data
df <- read_excel('PATH/TO/protein_peptide-workfile.xlsx') %>% 
  dplyr::select(c('Description', 'Contaminant', 'Checked', 'Biological Process', 
                             'Cellular Component', 'Molecular Function',
                             'Reactome Pathways', 'Neutrophil Granule protein',
                             'Abundances (Grouped): spEV_Calpeptin',
                             'Abundances (Grouped): spEV_CytochalasinD',
                             'Abundances (Grouped): spEV_dH2O',
                             'Abundances (Grouped): spEV_DMSO',
                             'Abundances (Grouped): spEV_GW4869',
                             'Abundances (Grouped): spEV_Nexinhib20',
                             'Abundances (Grouped): spEV_Pantethine',
                             'Abundances (Grouped): spEV_R5421',
                             'Abundances (Grouped): spEV_Y27632')) %>% 
                               filter(Contaminant == FALSE)
df <- df[!is.na(df$Checked), ]
colnames(df) <- c('protein_name', 'contaminant', 'checked', 'biological_process',
                  'cellular_component', 'molecular_function', 'reactome_pathways', 
                  'neutrophil_granule_protein',
                  'spEV + Calpeptin', 'spEV + CytochalasinD', 
                  'spEV (dH2O)', 'spEV (DMSO)', 'spEV + GW4869', 
                  'spEV + Nexinhib20', 'spEV + Pantethine', 
                  'spEV + R5421', 'spEV + Y27632')
df$protein_name <- gsub("\\s*\\[.*\\]$", "", df$protein_name)

# define sample sets
abundance_cols <- c('spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421', 'spEV + GW4869',
                    'spEV + Nexinhib20', 'spEV + CytochalasinD', 
                    'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632')

control_columns <- c('spEV (DMSO)', 'spEV (dH2O)')
case_columns <- c('spEV + Calpeptin', 'spEV + R5421', 'spEV + GW4869',
                  'spEV + Nexinhib20', 'spEV + CytochalasinD',
                  'spEV + Pantethine', 'spEV + Y27632')

sample_groups <- data.frame(
  Group = ifelse(
    abundance_cols %in% control_columns, 
    "Control", 
    "Case"), row.names = abundance_cols)

# ------------------------------------------------------------------------------
# Figure 9B: top 100 most abundant proteins ------------------------------------
df$total_abundance <- rowSums(df[abundance_cols], na.rm = TRUE)
top_proteins <- df %>%
  dplyr::arrange(desc(total_abundance)) %>%
  dplyr::slice_head(n = 100)

mat <- as.matrix(top_proteins[abundance_cols])
rownames(mat) <- top_proteins$protein_name

mat_log2 <- log2(mat + 1)

gt <- 
  pheatmap(mat_log2,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           fontsize_row = 6,
           fontsize_col = 10,
           angle_col = 45,
           main = "Top 100 most abundant proteins across treatments (log2 abundance)",
           annotation_col = sample_groups,
           annotation_colors = list(Group = c(Case = "#20A387FF", Control = "#333333")))$gtable

ggsave("../figs/protein_heatmap_top100_abundance_log2abundance.pdf", plot=gt, width = 10, height = 10)

# Figure 9C neutrophil granule proteins ----------------------------------------
neutro_proteins <- df %>% 
  filter(neutrophil_granule_protein == TRUE)

mat <- as.matrix(neutro_proteins[abundance_cols])
rownames(mat) <- neutro_proteins$protein_name
mat_log2 <- log2(mat + 1)

gt <- pheatmap(
  mat_log2,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  main = "Neutrophil granule proteins across treatments \n(log2 abundance)",
  annotation_col = sample_groups,
  annotation_colors = list(Group = c(Case = "#20A387FF", Control = "#333333"))
)$gtable

ggsave("../figs/protein_heatmap_neutrophil_granule_log2abundance.pdf", plot = gt, width = 10, height = 10)
