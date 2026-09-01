################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 06 - Escudos Tumorais (Exaustão, Checkpoints Imunes e HLA)
# Objetivo:
#   1. Estatística e Violinos Individuais de TODOS os genes de Exaustão e HLA
#   2. Heatmaps Complexos com Metadados Clínicos (Tumor Burden, CRS, NE, etc)
#   3. Estatística FDR (Asteriscos) embutida nas linhas dos Heatmaps
################################################################################

## =============================================================================
## 1. Pacotes e Diretórios
## =============================================================================
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
dir.create("results/exhaustion_hla/single_plots", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados Purificados (Fase 02)
## =============================================================================
message("Carregando matrizes purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
z7_sig  <- readRDS("data/filtered/ZUMA7_AxiCel_Sigs_Clean.rds")

pal_binaria <- c("CR" = "#4DBBD5FF", "Non-CR" = "#E64B35FF", "Ongoing" = "#00A087FF", "Others" = "#3C5488FF")
pal_granular <- c("CR" = "#4DBBD5FF", "PR" = "#F39B7FFF", "SD" = "#7E6148FF", "PD" = "#E64B35FF",
                  "Ongoing Response" = "#00A087FF", "Relapsed" = "#8491B4FF", "Nonresponders" = "#3C5488FF", "Missing" = "gray80")

## =============================================================================
## 3. Definição dos Grupos de Genes
## =============================================================================
exhaustion_signature_groups <- list(
  "Inhibitory receptors" = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD244", "BTLA", "CD160"),
  "Exhaustion TF program" = c("TOX", "TOX2", "EOMES", "TBX21", "BATF", "PRDM1", "TCF7", "CXCR5"),
  "Metabolic / Terminal" = c("ENTPD1", "NT5E", "CXCL13", "LAYN"),
  "Cytotoxic activation" = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GNLY", "PRF1", "NKG7", "MKI67")
)

hla_signature_groups <- list(
  "Classical HLA-I" = c("HLA-A", "HLA-B", "HLA-C"),
  "HLA-I processing" = c("B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5"),
  "Classical HLA-II" = c("HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
  "Non-classical HLA-I" = c("HLA-E", "HLA-F", "HLA-G")
)

# Compilar todos os genes para o loop de Violinos Individuais
all_target_genes <- unique(c(unlist(exhaustion_signature_groups), unlist(hla_signature_groups)))

## =============================================================================
## 4. Construção dos DataFrames para Gráficos Individuais
## =============================================================================
build_plot_df <- function(dataset_list) {
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  
  genes_present <- intersect(all_target_genes, rownames(expr))
  df_genes <- t(expr[genes_present, , drop = FALSE]) %>% as.data.frame() %>% rownames_to_column("sample")
  return(clin %>% left_join(df_genes, by = "sample"))
}

df_z1 <- build_plot_df(z1_base)
df_z7 <- build_plot_df(z7_mrna)

# Incorporar Checkpoints IO360 no ZUMA-7
io360_checkpoints <- c("CTLA4.IO360", "PD-L1.IO360", "PD-L2.IO360", "TIGIT.IO360", "Treg.IO360", "Exhausted CD8.IO360")
valid_io360 <- intersect(io360_checkpoints, rownames(z7_sig$expr))
df_z7_io360 <- t(z7_sig$expr[valid_io360, , drop = FALSE]) %>% as.data.frame() %>% rownames_to_column("sample")
df_z7 <- df_z7 %>% left_join(df_z7_io360, by = "sample")

## =============================================================================
## 5. Estatística e Violinos Individuais (Gene a Gene)
## =============================================================================
message("Gerando gráficos de Violino para todos os genes (isso pode demorar um pouco)...")

run_stats_and_plot <- function(df, features, cohort_name) {
  valid_features <- intersect(features, colnames(df))
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>% select(group = response_group, value = all_of(feature)) %>% filter(group != "Unknown") %>% drop_na()
    if (nlevels(as.factor(tmp$group)) != 2) return(NULL)
    test <- wilcox.test(value ~ as.factor(group), data = tmp)
    tibble(feature = feature, group_1 = levels(as.factor(tmp$group))[1], group_2 = levels(as.factor(tmp$group))[2], p_value = test$p.value)
  }) %>% mutate(padj_BH = p.adjust(p_value, method = "BH"), label_padj = paste0("FDR = ", signif(padj_BH, 2)))
  
  write_csv(stats, paste0("results/exhaustion_hla/", cohort_name, "_All_Genes_Stats.csv"))
  
  plot_single <- function(feature_name) {
    plot_df <- df %>% select(sample, group = response_group, response_original, score = all_of(feature_name)) %>%
      filter(group != "Unknown") %>% drop_na(score) %>%
      mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original))
    
    stat_row <- stats %>% filter(feature == feature_name) %>% mutate(group1 = group_1, group2 = group_2)
    y_pos <- max(plot_df$score, na.rm = TRUE) + 0.10 * diff(range(plot_df$score, na.rm = TRUE))
    
    p <- ggplot(plot_df, aes(x = group, y = score)) +
      geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", scale = "width") +
      geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
      geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85) +
      scale_fill_manual(values = pal_binaria) + scale_color_manual(values = pal_granular) +
      theme_classic(base_size = 14) + labs(title = feature_name, x = NULL, y = "Expression (Log2)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none")
    
    if(nrow(stat_row) > 0) p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, tip.length = 0.015)
    
    safe_filename <- gsub("[^A-Za-z0-9]", "_", feature_name)
    ggsave(paste0("results/exhaustion_hla/single_plots/", cohort_name, "_", safe_filename, "_Violin.pdf"), p, width = 4.5, height = 5)
  }
  walk(valid_features, plot_single)
}

run_stats_and_plot(df_z1, all_target_genes, "ZUMA1_Baseline")
run_stats_and_plot(df_z7, c(all_target_genes, valid_io360), "ZUMA7_mRNA")

## =============================================================================
## 6. Preparação dos Metadados Clínicos para o Heatmap (CORRIGIDO)
## =============================================================================
# ZUMA-1: Metadados
anno_z1 <- z1_base$clin %>%
  rownames_to_column("sample") %>% # 1. Salva os IDs das amostras numa coluna
  filter(response_group != "Unknown") %>%
  mutate(Tumor_Burden_Log10 = log10(as.numeric(baseline_tumour_burden_mm2_ch1) + 1)) %>%
  select(sample, # Garante que a coluna 'sample' continue no dataframe
         Response = response_group, 
         Clinical_Class = response_original, 
         Molecular_Subgroup = molecular_subgroup_ch1, 
         Neurotox_Grade = worst_grade_of_ne_ch1, 
         Tumor_Burden_Log10) %>%
  column_to_rownames("sample") # 2. Transforma a coluna de volta no nome das linhas

# ZUMA-7: Metadados (Corrigindo nomes sujos do GEO)
anno_z7 <- z7_mrna$clin %>%
  rownames_to_column("sample") %>% # 1. Salva os IDs
  filter(response_group != "Unknown") %>%
  mutate(
    Tumor_Burden_Log10 = log10(as.numeric(baseline_tumor_burden_spd_ch1) + 1),
    EFS_Event = as.factor(event_free_survival_event_ch1)
  ) %>%
  select(sample,
         Response = response_group, 
         Clinical_Class = response_original, 
         Cell_of_Origin = cell_of_origin_ch1, 
         Grade3_NE = grade3_ne_ch1, 
         Grade3_CRS = grade3_crs_ch1, 
         EFS_Event, 
         Tumor_Burden_Log10) %>%
  column_to_rownames("sample") # 2. Transforma de volta
## =============================================================================
## 7. Função do Heatmap (SIMPLIFICADA - SEM ESTATÍSTICA LATERAL)
## =============================================================================
zscore_and_cap <- function(mat, cap = 2) {
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled[mat_scaled > cap] <- cap
  mat_scaled[mat_scaled < -cap] <- -cap
  return(mat_scaled)
}

build_complex_heatmap_v2 <- function(dataset, signature_list, anno_df, title, prefix) {
  
  expr <- dataset$expr
  genes_present <- unlist(signature_list, use.names = FALSE) %>% intersect(rownames(expr))
  
  if(length(genes_present) == 0) return(NULL)
  
  # 1. Alinhamento de dados
  valid_samples <- intersect(rownames(anno_df), colnames(expr))
  if(length(valid_samples) == 0) return(NULL)
  
  anno_df <- anno_df[valid_samples, , drop = FALSE]
  anno_df <- anno_df[!is.na(anno_df$Response), , drop = FALSE]
  valid_samples <- rownames(anno_df)
  
  mat_scaled <- zscore_and_cap(expr[genes_present, valid_samples, drop = FALSE])
  message("✅ Desenhando ", prefix, " | Matriz: ", nrow(mat_scaled), " genes x ", ncol(mat_scaled), " amostras.")
  
  # 2. Anotação Topo (Metadados Variáveis)
  min_tb <- min(anno_df$Tumor_Burden_Log10, na.rm = TRUE)
  max_tb <- max(anno_df$Tumor_Burden_Log10, na.rm = TRUE)
  if(is.infinite(min_tb) | is.na(min_tb)) min_tb <- 0 
  if(is.infinite(max_tb) | is.na(max_tb)) max_tb <- 1
  col_burden <- colorRamp2(c(min_tb, max_tb), c("white", "purple"))
  
  resp_levels <- unique(as.character(anno_df$Response))
  class_levels <- unique(as.character(anno_df$Clinical_Class))
  pal_resp <- pal_binaria[names(pal_binaria) %in% resp_levels]
  pal_class <- pal_granular[names(pal_granular) %in% class_levels]
  
  ha_list <- list(Response = pal_resp, Clinical_Class = pal_class, Tumor_Burden_Log10 = col_burden)
  top_ha <- HeatmapAnnotation(df = anno_df, col = ha_list, na_col = "gray90", annotation_name_side = "left")
  
  # 3. Organização Biológica (Divisão de Linhas e Colunas)
  gene_group <- rep(names(signature_list), lengths(signature_list))
  names(gene_group) <- unlist(signature_list, use.names = FALSE)
  valid_groups <- intersect(names(signature_list), unique(gene_group[genes_present]))
  row_split <- factor(gene_group[genes_present], levels = valid_groups)
  
  if(length(resp_levels) < 2) {
    column_split <- NULL
  } else {
    column_split <- factor(anno_df$Response, levels = resp_levels)
  }
  
  col_fun_heatmap <- colorRamp2(c(-2, 0, 2), c("#313695", "#FFFFBF", "#A50026"))
  
  # 4. Renderização Limpa e Direta
  pdf(paste0("results/exhaustion_hla/", prefix, "_Heatmap.pdf"), width = 10, height = 8)
  
  ht <- Heatmap(
    mat_scaled, name = "Z-score", col = col_fun_heatmap,
    top_annotation = top_ha, 
    row_split = row_split, column_split = column_split,
    row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
    cluster_rows = FALSE, cluster_column_slices = FALSE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9), row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = title
  )
  
  invisible(draw(ht, merge_legends = TRUE))
  dev.off()
}

# Gerar os Heatmaps
build_complex_heatmap_v2(z1_base, exhaustion_signature_groups, anno_z1, "ZUMA-1 T-cell Exhaustion", "ZUMA1_Exhaustion")
build_complex_heatmap_v2(z7_mrna, exhaustion_signature_groups, anno_z7, "ZUMA-7 T-cell Exhaustion", "ZUMA7_Exhaustion")
build_complex_heatmap_v2(z1_base, hla_signature_groups, anno_z1, "ZUMA-1 HLA & Antigen Presentation", "ZUMA1_HLA")
build_complex_heatmap_v2(z7_mrna, hla_signature_groups, anno_z7, "ZUMA-7 HLA & Antigen Presentation", "ZUMA7_HLA")
