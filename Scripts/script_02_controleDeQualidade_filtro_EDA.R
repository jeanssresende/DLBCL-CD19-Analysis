################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 02 – Controle de Qualidade, Filtragem Clínica e EDA Global
# Datasets: GSE197977 (ZUMA-1) e GSE248835 (ZUMA-7)
################################################################################

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

options(stringsAsFactors = FALSE)

# Criar diretórios para resultados
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/qc")) dir.create("results/qc", recursive = TRUE)
if (!dir.exists("data/filtered")) dir.create("data/filtered", recursive = TRUE)

## =============================================================================
## 1. Carregamento dos Dados Processados
## =============================================================================
message("Carregando dados processados...")
zuma1      <- readRDS("./data/ZUMA1_GSE197977_processed.rds")
zuma7_mrna <- readRDS("./data/ZUMA7_GSE248835_mRNA_processed.rds")
zuma7_sigs <- readRDS("./data/ZUMA7_GSE248835_Signatures_processed.rds")

## =============================================================================
## 2. Harmonização Clínica e Filtros Estritos (Baseline vs Post & Braço CAR-T)
## =============================================================================
message("Aplicando filtros clínicos e binarizando desfechos...")

# --- ZUMA-1 (BASELINE ESTRITO) ---
idx_z1_base <- grep("^BASELINE$", zuma1$clinical$visit_ch1, ignore.case = TRUE)
z1_clin_base <- zuma1$clinical[idx_z1_base, ] %>%
  mutate(
    response_group = case_when(
      grepl("CR", bestresponse_ch1, ignore.case = TRUE) ~ "CR",
      grepl("PR|SD|PD", bestresponse_ch1, ignore.case = TRUE) ~ "Non-CR",
      TRUE ~ "Unknown"
    ),
    response_original = bestresponse_ch1 
  )
z1_expr_base <- zuma1$expression_log2[, idx_z1_base]

# --- ZUMA-1 (POST-INFUSION / DAY 7-14 & WK4) ---
idx_z1_post <- grep("DAY7-14|FCBWK4", zuma1$clinical$visit_ch1, ignore.case = TRUE)
z1_clin_post <- zuma1$clinical[idx_z1_post, ] %>%
  mutate(
    response_group = case_when(
      grepl("CR", bestresponse_ch1, ignore.case = TRUE) ~ "CR",
      grepl("PR|SD|PD", bestresponse_ch1, ignore.case = TRUE) ~ "Non-CR",
      TRUE ~ "Unknown"
    ),
    response_original = bestresponse_ch1 
  )
z1_expr_post <- zuma1$expression_log2[, idx_z1_post]


# --- ZUMA-7 (mRNA e Signatures) ---
idx_z7_axi <- grep("Axicabtagene Ciloleucel", zuma7_mrna$clinical$treatment_arm_ch1, ignore.case = TRUE)
z7_clin_axi <- zuma7_mrna$clinical[idx_z7_axi, ] %>%
  mutate(
    response_group = case_when(
      ongoing_response_ch1 == "Ongoing Response" ~ "Ongoing",
      ongoing_response_ch1 %in% c("Nonresponders", "Relapsed") ~ "Others",
      TRUE ~ "Unknown"
    ),
    response_original = ongoing_response_ch1
  )
z7_expr_axi <- zuma7_mrna$expression_log2[, idx_z7_axi]
z7_expr_sigs <- zuma7_sigs$expression_log2[, idx_z7_axi]

## =============================================================================
## 3. Funções de QC de Distribuição (Boxplot & Density)
## =============================================================================
message("Gerando gráficos de distribuição global...")

plot_qc_distribution <- function(expr_matrix, dataset_name) {
  expr_long <- as.data.frame(expr_matrix) %>% 
    pivot_longer(cols = everything(), names_to = "Sample", values_to = "Log2_Expression")
  
  p_box <- ggplot(expr_long, aes(x = Sample, y = Log2_Expression)) +
    geom_boxplot(outlier.size = 0.5, color = "#2c3e50", fill = "#bdc3c7", alpha = 0.7) +
    theme_classic() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = paste(dataset_name, "- Boxplot Global"), x = "Amostras", y = "Expressão log2")
  
  p_dens <- ggplot(expr_long, aes(x = Log2_Expression, group = Sample)) +
    geom_density(alpha = 0.15, color = "#2980b9") +
    theme_classic() +
    labs(title = paste(dataset_name, "- Densidade Global"), x = "Expressão log2", y = "Densidade")
  
  ggsave(paste0("results/qc/", dataset_name, "_boxplot.pdf"), p_box, width = 10, height = 5)
  ggsave(paste0("results/qc/", dataset_name, "_density.pdf"), p_dens, width = 8, height = 5)
}

plot_qc_distribution(z1_expr_base, "ZUMA1_Baseline")
plot_qc_distribution(z1_expr_post, "ZUMA1_Post")
plot_qc_distribution(z7_expr_axi, "ZUMA7_AxiCel_mRNA")
plot_qc_distribution(z7_expr_sigs, "ZUMA7_AxiCel_Signatures")

## =============================================================================
## 4. Correlação entre amostras (Spearman com Escala Absoluta de -1 a 1)
## =============================================================================
message("Calculando matrizes de correlação globais...")

calculate_sample_correlation_pdf <- function(expr_matrix, dataset_name) {
  sample_cor <- cor(expr_matrix, method = "spearman", use = "pairwise.complete.obs")
  my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
  my_breaks <- seq(-1, 1, length.out = 51)
  
  pdf(file = paste0("results/qc/", dataset_name, "_correlation_spearman.pdf"), width = 10, height = 8)
  pheatmap(sample_cor, main = paste(dataset_name, "- Correlação Global (Spearman)"), 
           show_rownames = FALSE, show_colnames = FALSE, 
           color = my_colors, breaks = my_breaks) 
  dev.off()
}

calculate_sample_correlation_pdf(z1_expr_base, "ZUMA1_Baseline")
calculate_sample_correlation_pdf(z1_expr_post, "ZUMA1_Post")
calculate_sample_correlation_pdf(z7_expr_axi, "ZUMA7_AxiCel_mRNA")
calculate_sample_correlation_pdf(z7_expr_sigs, "ZUMA7_AxiCel_Signatures")

## =============================================================================
## 5. Análise de Componentes Principais (PCA por Variáveis Clínicas)
## =============================================================================
message("Rodando PCAs por agrupamento clínico...")

plot_pca <- function(expr_matrix, clinical_data, color_var, dataset_name) {
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  pca_df <- data.frame(Sample = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2])
  clin_sub <- clinical_data %>% rownames_to_column("Sample")
  pca_df <- left_join(pca_df, clin_sub, by = "Sample")
  
  p <- ggplot(pca_df, aes_string(x = "PC1", y = "PC2", color = color_var)) +
    geom_point(size = 3.5, alpha = 0.85) +
    theme_classic() +
    theme(legend.position = "right", plot.title = element_text(face="bold", hjust=0.5)) +
    labs(title = paste(dataset_name, "- PCA por", color_var),
         x = paste0("PC1 (", var_explained[1], "%)"), y = paste0("PC2 (", var_explained[2], "%)"))
  
  ggsave(paste0("results/qc/", dataset_name, "_PCA_", color_var, ".pdf"), p, width = 7, height = 5)
}

plot_pca(z1_expr_base, z1_clin_base, "response_group", "ZUMA1_Baseline")
#plot_pca(z1_expr_post, z1_clin_post, "response_group", "ZUMA1_Post")
plot_pca(z7_expr_axi, z7_clin_axi, "response_group", "ZUMA7_AxiCel_mRNA")
plot_pca(z7_expr_sigs, z7_clin_axi, "response_group", "ZUMA7_AxiCel_Signatures")

# --- PCA LONGITUDINAL (Avaliação do Impacto Clínico Temporal do CAR-T) ---
message("Gerando PCA longitudinal clínica (ZUMA-1)...")
common_genes <- intersect(rownames(z1_expr_base), rownames(z1_expr_post))
z1_combined_expr <- cbind(z1_expr_base[common_genes,], z1_expr_post[common_genes,])

z1_combined_clin <- rbind(
  z1_clin_base %>% select(geo_accession, response_group, visit_ch1),
  z1_clin_post %>% select(geo_accession, response_group, visit_ch1)
)

pca_long <- prcomp(t(z1_combined_expr), scale. = TRUE)
var_long <- round(100 * (pca_long$sdev^2 / sum(pca_long$sdev^2)), 1)
df_long <- data.frame(Sample = rownames(pca_long$x), PC1 = pca_long$x[, 1], PC2 = pca_long$x[, 2]) %>%
  left_join(z1_combined_clin, by = c("Sample" = "geo_accession"))

p_long <- ggplot(df_long, aes(x = PC1, y = PC2, color = visit_ch1, shape = response_group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = values <- c("BASELINE" = "#16a085", "DAY7-14" = "#d35400", "FCBWK4" = "#8e44ad")) +
  theme_classic() +
  labs(title = "ZUMA-1: Impacto Temporal do CAR-T no Transcriptoma Global",
       x = paste0("PC1 (", var_long[1], "%)"), y = paste0("PC2 (", var_long[2], "%)"))

ggsave("results/qc/ZUMA1_Longitudinal_PCA.pdf", p_long, width = 8, height = 6)

## =============================================================================
## 6. Agrupamento Não-Supervisionado de Amostras (Heatmap de Variância)
## =============================================================================
message("Gerando Heatmaps de Agrupamento de Amostras (Versão com Proteção de Variância)...")

plot_top_variable_genes <- function(expr_matrix, clinical_data, dataset_name, annot_cols, top_n = 500) {
  
  # 1. Calcular a variância por gene e remover estritamente variância zero ou NAs
  gene_var <- apply(expr_matrix, 1, var, na.rm = TRUE)
  gene_var <- gene_var[!is.na(gene_var) & gene_var > 0]
  
  # Se a coorte tiver menos genes válidos que o top_n (comum em subgrupos pequenos), ajustamos o teto
  n_genes_validos <- length(gene_var)
  top_n_actual <- min(top_n, n_genes_validos)
  
  if (top_n_actual == 0) {
    warning(paste("Aviso: Nenhum gene com variância maior que zero encontrado em:", dataset_name))
    return(NULL)
  }
  
  # Selecionar os top genes baseados na variância real existente
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n_actual]
  expr_top <- expr_matrix[top_genes, , drop = FALSE]
  
  # Remoção estrita de qualquer linha com NA remanescente para blindar o hclust
  expr_top <- na.omit(expr_top)
  
  # 2. Garantir que os dados clínicos usem o geo_accession como Rownames (Alinhamento estrito)
  clinical_df <- as.data.frame(clinical_data)
  if ("geo_accession" %in% colnames(clinical_df)) {
    rownames(clinical_df) <- clinical_df$geo_accession
  }
  
  # 3. Selecionar apenas as colunas de metadados clínicos desejadas
  annotation_df <- clinical_df[, annot_cols, drop = FALSE]
  
  # 4. Tratamento de strings vazias ou nulas nos metadados para evitar falhas de legenda
  annotation_df <- annotation_df %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ case_when(
      is.na(.) ~ "Unknown",
      . == "" ~ "Unknown",
      toupper(.) == "NA" ~ "Unknown",
      tolower(.) == "null" ~ "Unknown",
      toupper(.) == "UNCLASSIFIED" ~ "Unknown",
      TRUE ~ .
    ))) %>%
    as.data.frame()
  
  rownames(annotation_df) <- rownames(clinical_df)
  
  # 5. Forçar o alinhamento exato com as colunas da matriz de expressão fatiada
  annotation_df <- annotation_df[colnames(expr_top), , drop = FALSE]
  
  # Configuração de cores e quebras simétricas de Z-score
  my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
  my_breaks <- seq(-3, 3, length.out = 51)
  
  # 6. Geração do arquivo PDF protegida por um capturador de erros (tryCatch)
  pdf(file = paste0("results/qc/", dataset_name, "_SampleClustering_Top", top_n, ".pdf"), width = 10, height = 8)
  
  tryCatch({
    pheatmap(expr_top, 
             scale = "row", 
             show_rownames = FALSE, 
             show_colnames = FALSE,
             annotation_col = annotation_df, 
             main = paste0(dataset_name, " - Agrupamento por Variância (Top ", top_n_actual, " Genes)"),
             color = my_colors, 
             breaks = my_breaks)
  }, error = function(e) {
    message(paste("Erro crítico detectado ao plotar o heatmap de", dataset_name, ":", e$message))
  })
  
  dev.off()
}

# Execução das chamadas com as anotações mapeadas
zuma1_annot <- c("response_group", "response_original", "molecular_subgroup_ch1")
zuma7_annot <- c("response_group", "response_original", "cell_of_origin_ch1")

plot_top_variable_genes(z1_expr_base, z1_clin_base, "ZUMA1_Baseline", zuma1_annot)
plot_top_variable_genes(z1_expr_post, z1_clin_post, "ZUMA1_Post", zuma1_annot)
plot_top_variable_genes(z7_expr_axi, z7_clin_axi, "ZUMA7_AxiCel_mRNA", zuma7_annot)

# --- HEATMAP DAS ASSINATURAS DO TME (ZUMA-7) ---
# Justificativa: Assinaturas não são genes isolados, são módulos celulares/clínicos macro.
message("Gerando Heatmap das Módulos Clínicos do TME (ZUMA-7)...")
z7_expr_sigs <- zuma7_sigs$expression_log2[, idx_z7_axi]

pdf(file = "results/qc/ZUMA7_AxiCel_All_IO360_Signatures.pdf", width = 11, height = 9)
pheatmap(z7_expr_sigs, scale = "row", show_rownames = TRUE, show_colnames = FALSE,
         annotation_col = z7_clin_axi[, zuma7_annot, drop=FALSE],
         main = "ZUMA-7: Perfis Arquiteturais do TME (IO360 Signatures)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         breaks = seq(-3, 3, length.out = 51), fontsize_row = 8)
dev.off()

## =============================================================================
## 7. Exportação das Matrizes Purificadas
## =============================================================================
message("Salvando dados limpos...")
z1_base_clean <- list(expr = z1_expr_base, clin = z1_clin_base)
z1_post_clean <- list(expr = z1_expr_post, clin = z1_clin_post)
z7_mrna_clean <- list(expr = z7_expr_axi, clin = z7_clin_axi)
z7_sigs_clean <- list(expr = z7_expr_sigs, clin = z7_clin_axi)

saveRDS(z1_base_clean, file = "data/filtered/ZUMA1_Baseline_Clean.rds")
saveRDS(z1_post_clean, file = "data/filtered/ZUMA1_Post_Clean.rds")
saveRDS(z7_mrna_clean, file = "data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
saveRDS(z7_sigs_clean, file = "data/filtered/ZUMA7_AxiCel_Sigs_Clean.rds")

cat("\n==========================================================\n")
cat("Controle de Qualidade e EDA Clínica Concluídos!\n")
cat("Nenhum gene individual foi avaliado nesta etapa.\n")
cat("==========================================================\n")

