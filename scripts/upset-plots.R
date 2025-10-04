library(UpSetR)
library(openxlsx)

# miRNA ------------------------------------------------------------------------
# read data and standardize columns
df <- read.xlsx("PATH/TO/mirna_summary_RPMs.xlsx") # input placeholder
df$miRNA <- sub("^hsa-", "", df$miRNA)
colnames(df) <- c('miRNA', 'spEV + Calpeptin', 'spEV + CytochalasinD', 'spEV (dH2O)',
                    'spEV (DMSO)','spEV + GW4869',  'spEV + Nexinhib20', 
                    'spEV + Pantethine', 'spEV + R5421', 'spEV + Y27632')
# choose set order and build presence/absence matrix (RPM != 0)
abundance_cols <- c('spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421', 'spEV + GW4869',
                    'spEV + Nexinhib20', 'spEV + CytochalasinD', 
                    'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632')
binary_mat <- as.data.frame(df[abundance_cols] != 0)
binary_mat$miRNA <- df$miRNA
binary_mat[abundance_cols] <- lapply(binary_mat[abundance_cols], as.integer)

# plot UpSet with highlights
abundance_cols_upset <- rev(c(
  'spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421',
  'spEV + GW4869', 'spEV + Nexinhib20', 'spEV + CytochalasinD',
  'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632'
))
pdf("../figs/mirna_upset_highlighted.pdf", width = 20, height = 10)
UpSetR::upset(binary_mat,
      sets = abundance_cols_upset,
      text.scale = 2.3,
      keep.order = TRUE,
      order.by = "freq",
      queries = list(
        list(query = intersects, params = list("spEV + Calpeptin"), color = "#34FEAD", active = TRUE),
        list(query = intersects, params = list("spEV + Pantethine"), color = "#BC9BFE", active = TRUE),
        list(query = intersects, params = list("spEV + R5421"), color = "#F8B939", active = TRUE),
        list(query = intersects, params = list("spEV + GW4869"), color = "#688DFF", active = TRUE),
        list(query = intersects, params = list("spEV + CytochalasinD"), color = "#B97B3D", active = TRUE),
        list(query = intersects, params = list("spEV (dH2O)"), color = "#404040", active = TRUE),
        list(query = intersects, params = list("spEV + Nexinhib20"), color = "#FE6161", active = TRUE),
        list(query = intersects, params = list("spEV + Y27632"), color = "#E78AC3", active = TRUE)
      )
)
grid.text("miRNA Overlap Across EV Samples", 
          x = 0.5, y = 0.99, gp = gpar(fontsize = 19, fontface = "bold"), just = "center")
dev.off()

# mRNA -------------------------------------------------------------------------
# load expression, keep relevant columns, and drop predicted accessions
df <- read_csv('PATH/TO/rna_expression.csv') # input placeholder
df <- df %>% dplyr::select(c('rna_id', 'gene_id', 'tpm_spEV_DMSO', 'tpm_spEV_50ug_ml_Calpeptin', 
                             'tpm_spEV_50uM_R5421', 'tpm_spEV_10uM_GW4869',
                             'tpm_spEV_10uM_Nexinhib20', 'tpm_spEV_1uM_CytochalasinD',
                             'tpm_spEV_dH2O', 'tpm_spEV_10uM_Pantethine',
                             'tpm_spEV_10uM_Y27632'))
tpm_cols <- c('tpm_spEV_DMSO', 'tpm_spEV_50ug_ml_Calpeptin', 
              'tpm_spEV_50uM_R5421', 'tpm_spEV_10uM_GW4869',
              'tpm_spEV_10uM_Nexinhib20', 'tpm_spEV_1uM_CytochalasinD',
              'tpm_spEV_dH2O', 'tpm_spEV_10uM_Pantethine',
              'tpm_spEV_10uM_Y27632')
filtered_df <- df %>% filter(!grepl("^XM_|^rna-XM_", rna_id))
filtered_df$rna_id_clean <- sub("^rna-(NM_)", "\\1", filtered_df$rna_id)

# map ENTREZ IDs to gene symbol/name and remove hemoglobin entries
entrez_ids <- as.character(filtered_df$gene_id)
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

filtered_df$gene_symbol <- gene_symbols
filtered_df$gene_name <- gene_names
filtered_df <- filtered_df[!grepl("hemoglobin", filtered_df$gene_name, ignore.case = TRUE), ]

# relabel columns for plotting, build presence/absence (TPM > 0), and order sets
set_labels <- c(
  'tpm_spEV_DMSO' = 'spEV (DMSO)',
  'tpm_spEV_50ug_ml_Calpeptin' = 'spEV + Calpeptin', 
  'tpm_spEV_50uM_R5421' = 'spEV + R5421', 
  'tpm_spEV_10uM_GW4869' = 'spEV + GW4869',
  'tpm_spEV_10uM_Nexinhib20' = 'spEV + Nexinhib20', 
  'tpm_spEV_1uM_CytochalasinD' = 'spEV + CytochalasinD',
  'tpm_spEV_dH2O' = 'spEV (dH2O)', 
  'tpm_spEV_10uM_Pantethine' = 'spEV + Pantethine',
  'tpm_spEV_10uM_Y27632' = 'spEV + Y27632'
)
binary_mat <- as.data.frame(filtered_df[tpm_cols] > 0)
binary_mat[tpm_cols] <- lapply(binary_mat[tpm_cols], as.integer)
colnames(binary_mat)[colnames(binary_mat) %in% names(set_labels)] <- set_labels[colnames(binary_mat)[colnames(binary_mat) %in% names(set_labels)]]
tpm_cols_upset <- unname(set_labels[tpm_cols])
binary_mat <- binary_mat[, tpm_cols_upset]

# plot UpSet with highlights
pdf("../figs/mrna_upset_highlighted.pdf", width = 20, height = 10)
UpSetR::upset(binary_mat,
              sets = rev(tpm_cols_upset),
              text.scale = 1.75,
              keep.order = TRUE,
              order.by = "freq",
              queries = list(
                list(query = intersects, params = list("spEV + Calpeptin"), color = "#37fea0", active = TRUE),
                list(query = intersects, params = list("spEV + Pantethine"), color = "#e3539e", active = TRUE),
                list(query = intersects, params = list("spEV + R5421"), color = "#f6b000", active = TRUE),
                list(query = intersects, params = list("spEV + GW4869"), color = "#2880ff", active = TRUE),
                list(query = intersects, params = list("spEV + CytochalasinD"), color = "#b26f3b", active = TRUE),
                list(query = intersects, params = list("spEV (dH2O)"), color = "#404040", active = TRUE),
                list(query = intersects, params = list("spEV + Nexinhib20"), color = "#ff4f4e", active = TRUE),
                list(query = intersects, params = list("spEV + Y27632"), color = "#b500fe", active = TRUE)
              )
)
grid.text("mRNA Overlap Across EV Samples", 
          x = 0.5, y = 0.99, gp = gpar(fontsize = 19, fontface = "bold"), just = "center")
dev.off()

# protein ----------------------------------------------------------------------
# load and tidy protein sheet, keep non-contaminants and standardize columns
df <- read_excel('PATH/TO/protein_peptide-workfile.xlsx') # input placeholder
df <- df[!is.na(df$Checked), ]
df <- df %>% dplyr::select(c('Description', 'Contaminant', 'Biological Process', 
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
                             'Abundances (Grouped): spEV_Y27632'))
colnames(df) <- c('protein_name', 'contaminant', 'biological_process',
                  'cellular_component', 'molecular_function', 'reactome_pathways', 'neutrophil_granule_protein',
                  'spEV + Calpeptin', 'spEV + CytochalasinD', 
                  'spEV (dH2O)', 'spEV (DMSO)', 'spEV + GW4869', 
                  'spEV + Nexinhib20', 'spEV + Pantethine', 
                  'spEV + R5421', 'spEV + Y27632')
df$protein_name <- gsub("\\s*\\[.*\\]$", "", df$protein_name)
df <- df %>% filter(contaminant == FALSE)

# build presence/absence (non-NA)
abundance_cols <- c('spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421', 'spEV + GW4869',
                    'spEV + Nexinhib20', 'spEV + CytochalasinD', 
                    'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632')

binary_mat <- as.data.frame(!is.na(df[abundance_cols]))
binary_mat$protein <- df$protein_name
binary_mat[abundance_cols] <- lapply(binary_mat[abundance_cols], as.integer)

abundance_cols_upset <- rev(c(
  'spEV (DMSO)', 'spEV + Calpeptin', 'spEV + R5421',
  'spEV + GW4869', 'spEV + Nexinhib20', 'spEV + CytochalasinD',
  'spEV (dH2O)', 'spEV + Pantethine', 'spEV + Y27632'
))

# plot UpSet with highlights
pdf("../figs/protein_upset_highlighted.pdf", width = 20, height = 10)
UpSetR::upset(binary_mat,
      sets = abundance_cols_upset,
      text.scale = 2.3,
      keep.order = TRUE,
      order.by = "freq",
      queries = list(
        list(query = intersects, params = list("spEV + Calpeptin"), color = "#34FEAD", active = TRUE),
        list(query = intersects, params = list("spEV + Pantethine"), color = "#BC9BFE", active = TRUE),
        list(query = intersects, params = list("spEV + R5421"), color = "#F8B939", active = TRUE),
        list(query = intersects, params = list("spEV + GW4869"), color = "#688DFF", active = TRUE),
        list(query = intersects, params = list("spEV + CytochalasinD"), color = "#B97B3D", active = TRUE),
        list(query = intersects, params = list("spEV (dH2O)"), color = "#404040", active = TRUE),
        list(query = intersects, params = list("spEV + Nexinhib20"), color = "#FE6161", active = TRUE)
      )
)
grid.text("Protein Overlap Across EV Samples", 
          x = 0.5, y = 0.99, gp = gpar(fontsize = 19, fontface = "bold"), just = "center")
dev.off()

