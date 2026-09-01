################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 03 - Microambiente e Composição Celular (Deconvolução)
# Objetivo:
#   1. Estimar composição celular por MCP-counter
#   2. Calcular scores manuais baseados em assinaturas
#   3. Criar proxies de conteúdo B/tumoral e razão tumor-imune
#   4. Comparar grupos clínicos (Wilcoxon + FDR)
#   5. Gerar heatmaps e boxplots
################################################################################

## =============================================================================
## 1. Pacotes
## =============================================================================
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(patchwork)
library(MCPcounter)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

# Criar diretório de saída
dir.create("results/cellular_composition", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar objetos curados (Do Script 02)
## =============================================================================
message("Carregando matrizes purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z1_post <- readRDS("data/filtered/ZUMA1_Post_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")

## =============================================================================
## 3. Paletas de cores unificadas
## =============================================================================
pal_clinica <- c(
  "CR" = brewer.pal(3, "Set1")[2],       # Azul
  "Non-CR" = brewer.pal(3, "Set1")[1],   # Vermelho
  "Ongoing" = brewer.pal(3, "Set2")[1],  # Verde
  "Others" = brewer.pal(3, "Set2")[2]    # Laranja
)

heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

## =============================================================================
## 4. Assinaturas manuais de interesse
## =============================================================================
signature_list <- list(
  B_lineage_proxy = c("CD19", "MS4A1", "CD22", "CD79A", "CD79B", "PAX5", "POU2AF1", "BTK", "BLK", "FCRL2"),
  T_cell = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "LCK", "ZAP70"),
  Cytotoxic = c("CD8A", "PRF1", "GNLY", "GZMA", "GZMB", "GZMH", "NKG7"),
  T_cell_exhaustion = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "ENTPD1", "EOMES", "TBX21"),
  Myeloid = c("CD68", "CD163", "CSF1R", "ITGAM", "FCGR2A", "FCGR3A", "MRC1", "C1QA", "C1QB", "TREM2"),
  HLA_class_I = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5"),
  HLA_class_II = c("HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "HLA-DRB1"),
  IFN_gamma = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "GBP2")
)

## =============================================================================
## 5. Deconvolução por MCP-counter
## =============================================================================
run_mcp_counter <- function(dataset_list, cohort_name) {
  message("Rodando MCP-counter para: ", cohort_name)
  
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  
  mcp_scores <- MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
  
  mcp_df <- as.data.frame(t(mcp_scores)) %>%
    rownames_to_column("sample") %>%
    left_join(clin, by = "sample") %>%
    mutate(cohort = cohort_name)
  
  write_csv(mcp_df, paste0("results/cellular_composition/", cohort_name, "_MCPcounter_scores.csv"))
  return(mcp_df)
}

mcp_z1_base <- run_mcp_counter(z1_base, "ZUMA1_Baseline")
mcp_z1_post <- run_mcp_counter(z1_post, "ZUMA1_Post")
mcp_z7_mrna <- run_mcp_counter(z7_mrna, "ZUMA7_AxiCel_mRNA")

## =============================================================================
## 6. Cálculo de Scores Manuais
## =============================================================================
calculate_signature_scores <- function(dataset_list, cohort_name, signatures) {
  message("Calculando scores manuais para: ", cohort_name)
  
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  score_df <- data.frame(sample = colnames(expr))
  
  for (sig_name in names(signatures)) {
    genes_available <- intersect(signatures[[sig_name]], rownames(expr))
    if (length(genes_available) >= 2) {
      score_df[[sig_name]] <- colMeans(expr[genes_available, , drop = FALSE], na.rm = TRUE)
    } else {
      score_df[[sig_name]] <- NA_real_
    }
  }
  
  score_df <- score_df %>%
    mutate(
      immune_proxy = rowMeans(select(., any_of(c("T_cell", "Cytotoxic", "T_cell_exhaustion", "Myeloid", "HLA_class_I", "HLA_class_II", "IFN_gamma"))), na.rm = TRUE),
      tumor_immune_ratio_proxy = B_lineage_proxy - immune_proxy
    ) %>%
    left_join(clin, by = "sample") %>%
    mutate(cohort = cohort_name)
  
  write_csv(score_df, paste0("results/cellular_composition/", cohort_name, "_Manual_Scores.csv"))
  return(score_df)
}

sig_z1_base <- calculate_signature_scores(z1_base, "ZUMA1_Baseline", signature_list)
sig_z1_post <- calculate_signature_scores(z1_post, "ZUMA1_Post", signature_list)
sig_z7_mrna <- calculate_signature_scores(z7_mrna, "ZUMA7_AxiCel_mRNA", signature_list)

## =============================================================================
## 7. Estatística de Comparação de Grupos (Wilcoxon + FDR)
## =============================================================================
run_group_stats <- function(df, feature_cols, cohort_name, analysis_name) {
  
  # INTERSEÇÃO DE SEGURANÇA: Só tenta testar as células que realmente existem no painel
  valid_features <- intersect(feature_cols, colnames(df))
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>%
      select(group = response_group, value = all_of(feature)) %>%
      filter(!is.na(group), !is.na(value), group != "Unknown") %>%
      mutate(group = as.factor(group))
    
    if (nlevels(tmp$group) != 2) return(NULL)
    
    groups <- levels(tmp$group)
    test <- wilcox.test(value ~ group, data = tmp)
    
    summary_group <- tmp %>%
      group_by(group) %>%
      summarise(median = median(value, na.rm = TRUE), .groups = "drop")
    
    tibble(
      cohort = cohort_name,
      feature = feature,
      group_1 = groups[1], group_2 = groups[2],
      median_group_1 = summary_group$median[1],
      median_group_2 = summary_group$median[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(
      padj_BH = p.adjust(p_value, method = "BH"),
      label_padj = paste0("FDR = ", signif(padj_BH, 2))
    ) %>%
    arrange(padj_BH)
  
  write_csv(stats, paste0("results/cellular_composition/", cohort_name, "_", analysis_name, "_Stats.csv"))
  return(stats)
}

# Definir colunas numéricas de interesse
mcp_cols <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "B lineage", "Monocytic lineage", "Myeloid dendritic cells", "Endothelial cells", "Fibroblasts")
man_cols <- c(names(signature_list), "immune_proxy", "tumor_immune_ratio_proxy")

# Rodar Estatísticas
stat_mcp_z1b <- run_group_stats(mcp_z1_base, mcp_cols, "ZUMA1_Baseline", "MCP")
stat_mcp_z1p <- run_group_stats(mcp_z1_post, mcp_cols, "ZUMA1_Post", "MCP")
stat_mcp_z7  <- run_group_stats(mcp_z7_mrna, mcp_cols, "ZUMA7_mRNA", "MCP")

stat_man_z1b <- run_group_stats(sig_z1_base, man_cols, "ZUMA1_Baseline", "Manual")
stat_man_z1p <- run_group_stats(sig_z1_post, man_cols, "ZUMA1_Post", "Manual")
stat_man_z7  <- run_group_stats(sig_z7_mrna, man_cols, "ZUMA7_mRNA", "Manual")

## =============================================================================
## 8. Figuras Prontas para Publicação (Violin + Boxplot + Granular Jitter)
## =============================================================================
message("Gerando Figuras de Alta Qualidade (Violin Plots)...")

# 8.1 Paletas de Cores Sofisticadas (Estilo Nature/Cell)
# Cores para preencher o Violino/Boxplot (Grupos Binários)
pal_binaria <- c(
  "CR" = "#4DBBD5FF",       # Azul claro
  "Non-CR" = "#E64B35FF",   # Vermelho
  "Ongoing" = "#00A087FF",  # Verde esmeralda
  "Others" = "#3C5488FF"    # Azul escuro
)

# Cores para os pontos dos pacientes (Resposta Original Granular)
pal_granular <- c(
  "CR" = "#4DBBD5FF", 
  "PR" = "#F39B7FFF",       # Pêssego
  "SD" = "#7E6148FF",       # Marrom/Dourado
  "PD" = "#E64B35FF",       
  "Ongoing Response" = "#00A087FF", 
  "Relapsed" = "#8491B4FF", # Azul acinzentado
  "Nonresponders" = "#3C5488FF", 
  "Missing" = "gray80"      # Cinza para ausentes
)

# 8.2 A Função de Plotagem Avançada
plot_composition_violin <- function(df, stats_df, feature_cols, cohort_name, analysis_name) {
  
  valid_features <- intersect(feature_cols, colnames(df))
  
  # Preparação dos dados
  plot_df <- df %>%
    select(sample, group = response_group, response_original, all_of(valid_features)) %>%
    filter(group != "Unknown") %>%
    pivot_longer(cols = all_of(valid_features), names_to = "feature", values_to = "score") %>%
    drop_na(score)
  
  plot_df <- plot_df %>%
    mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original))
  
  # Cálculo da altura Global para a barra de P-value (Para funcionar bem com a escala fixa)
  global_max <- max(plot_df$score, na.rm = TRUE)
  global_range <- diff(range(plot_df$score, na.rm = TRUE))
  
  label_df <- plot_df %>%
    group_by(feature) %>%
    # A barra de p-value vai flutuar logo acima do maior ponto DENTRO daquele painel específico
    summarise(y.position = max(score) + (0.05 * global_range), .groups = "drop") %>%
    left_join(stats_df, by = "feature") %>%
    mutate(group1 = group_1, group2 = group_2) %>% 
    drop_na(padj_BH)
  
  # Construção da Figura
  p <- ggplot(plot_df, aes(x = group, y = score)) +
    
    geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", trim = FALSE, scale = "width") +
    geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
    geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85, stroke = 0.5) +
    
    # MUDANÇA AQUI: Removemos o scales="free_y". Agora todos dividem o mesmo eixo Y!
    facet_wrap(~feature, ncol = 3) +
    
    scale_fill_manual(values = pal_binaria, name = "Resposta Global:") +
    scale_color_manual(values = pal_granular, name = "Classificação Clínica:") +
    
    stat_pvalue_manual(label_df, label = "label_padj", y.position = "y.position", 
                       tip.length = 0.015, size = 3.5, fontface = "bold") +
    
    theme_classic(base_size = 14) +
    labs(title = paste(cohort_name, "-", analysis_name), x = NULL, y = "Estimativa de Abundância Celular") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  ggsave(paste0("results/cellular_composition/", cohort_name, "_", analysis_name, "_Violin.pdf"), p, width = 13, height = 11)
  ggsave(paste0("results/cellular_composition/", cohort_name, "_", analysis_name, "_Violin.png"), p, width = 13, height = 11, dpi = 300)
}

# 8.3 Execução das Plotagens
plot_composition_violin(mcp_z1_base, stat_mcp_z1b, mcp_cols, "ZUMA1_Baseline", "MCPcounter")
plot_composition_violin(sig_z1_base, stat_man_z1b, man_cols, "ZUMA1_Baseline", "Manual_Scores")
plot_composition_violin(mcp_z7_mrna, stat_mcp_z7, mcp_cols, "ZUMA7_mRNA", "MCPcounter")
plot_composition_violin(sig_z7_mrna, stat_man_z7, man_cols, "ZUMA7_mRNA", "Manual_Scores")

cat("\n==========================================================\n")
cat("Violin Plots de Alta Qualidade Gerados com Sucesso!\n")
cat("==========================================================\n")
