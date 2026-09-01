################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 03 - Microambiente e Composição Celular (Apenas MCP-counter)
# Objetivo:
#   1. Estimar composição celular por MCP-counter
#   2. Comparar Baseline ZUMA-1 (CR vs Non-CR)
#   3. Comparar ZUMA-1 Pareado (Post-CART vs Baseline)
#   4. Comparar ZUMA-7 Axi-cel (Ongoing vs Others)
#   5. Comparar ZUMA-7 SOC (Ongoing vs Others)
#   6. Gerar estatísticas (Wilcoxon + FDR) e gráficos de alta qualidade
################################################################################

## =============================================================================
## 1. Pacotes e Configurações
## =============================================================================
library(tidyverse)
library(ggpubr)
library(MCPcounter)
library(RColorBrewer)

options(stringsAsFactors = FALSE)
dir.create("results/cellular_composition", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Matrizes Purificadas
## =============================================================================
message("Carregando matrizes purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z1_post <- readRDS("data/filtered/ZUMA1_Post_Clean.rds")
z7_axi  <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
z7_soc  <- readRDS("data/filtered/ZUMA7_SOC_mRNA_Clean.rds")

# Padronização: Garantir que o SOC tenha a coluna 'response_group'
if("ongoing_2grps_ch1" %in% colnames(z7_soc$clin) && !"response_group" %in% colnames(z7_soc$clin)) {
  z7_soc$clin$response_group <- z7_soc$clin$ongoing_2grps_ch1
}

## =============================================================================
## 3. Paletas de Cores
## =============================================================================
pal_binaria <- c(
  "CR" = "#4DBBD5FF", "Non-CR" = "#E64B35FF",
  "Ongoing" = "#00A087FF", "Others" = "#3C5488FF"
)

pal_granular <- c(
  "CR" = "#4DBBD5FF", "PR" = "#F39B7FFF", "SD" = "#7E6148FF", "PD" = "#E64B35FF",
  "Ongoing Response" = "#00A087FF", "Relapsed" = "#8491B4FF", "Nonresponders" = "#3C5488FF", 
  "Missing" = "gray80"
)

# Cores para a análise pareada (Tempo)
pal_tempo <- c("Baseline" = "#8491B4FF", "Post-CART" = "#DC0000FF")

## =============================================================================
## 4. Deconvolução por MCP-counter
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
mcp_z7_axi  <- run_mcp_counter(z7_axi,  "ZUMA7_AxiCel")
mcp_z7_soc  <- run_mcp_counter(z7_soc,  "ZUMA7_SOC")

mcp_cols <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", 
              "B lineage", "Monocytic lineage", "Myeloid dendritic cells", 
              "Endothelial cells", "Fibroblasts", "Neutrophils")

## =============================================================================
## 5. Estatística de Comparação de Grupos (Não-Pareada)
## =============================================================================
run_unpaired_stats <- function(df, feature_cols, cohort_name) {
  
  valid_features <- intersect(feature_cols, colnames(df))
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>%
      select(group = response_group, value = all_of(feature)) %>%
      filter(!is.na(group), group != "Unknown", group != "Missing") %>%
      drop_na(value) %>%
      mutate(group = as.factor(group))
    
    if (nlevels(tmp$group) != 2) return(NULL)
    
    groups <- levels(tmp$group)
    test <- wilcox.test(value ~ group, data = tmp)
    
    summary_group <- tmp %>%
      group_by(group) %>%
      summarise(median = median(value, na.rm = TRUE), .groups = "drop")
    
    tibble(
      cohort = cohort_name, feature = feature,
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
    arrange(p_value)
  
  write_csv(stats, paste0("results/cellular_composition/", cohort_name, "_Unpaired_Stats.csv"))
  return(stats)
}

message("\nCalculando estatísticas não-pareadas...")
stat_z1_base <- run_unpaired_stats(mcp_z1_base, mcp_cols, "ZUMA1_Baseline")
stat_z7_axi  <- run_unpaired_stats(mcp_z7_axi,  mcp_cols, "ZUMA7_AxiCel")
stat_z7_soc  <- run_unpaired_stats(mcp_z7_soc,  mcp_cols, "ZUMA7_SOC")

## =============================================================================
## 6. Estatística Pareada (ZUMA-1 Post vs Baseline)
## =============================================================================

# 1. Extrair o ID do doente da coluna 'title' (Ex: "Sample 54 Patient 01" -> "Patient 01")
mcp_z1_base <- mcp_z1_base %>%
  mutate(patient_id = str_extract(title, "Patient \\d+"))

mcp_z1_post <- mcp_z1_post %>%
  mutate(patient_id = str_extract(title, "Patient \\d+"))

# 2. Verificar quantos pares perfeitos conseguimos formar
message("Calculando estatísticas pareadas (ZUMA-1)...")

# 3. Rodar a função pareada usando a nova coluna 'patient_id'
stat_z1_paired <- run_paired_stats(mcp_z1_base, mcp_z1_post, mcp_cols, id_col = "patient_id")

# ATENÇÃO: Ajuste 'id_col' para o nome exato da coluna que identifica o paciente
run_paired_stats <- function(df_base, df_post, feature_cols, id_col = "patient_id") {
  
  valid_features <- intersect(feature_cols, intersect(colnames(df_base), colnames(df_post)))
  
  stats <- map_dfr(valid_features, function(feature) {
    
    base_tmp <- df_base %>% select(all_of(id_col), base_val = all_of(feature))
    post_tmp <- df_post %>% select(all_of(id_col), post_val = all_of(feature))
    
    # Manter apenas pacientes que têm amostra nos dois tempos
    paired_df <- inner_join(base_tmp, post_tmp, by = id_col) %>% drop_na()
    
    if (nrow(paired_df) < 3) return(NULL) # Mínimo de n para testar
    
    # Teste de Wilcoxon pareado
    test <- wilcox.test(paired_df$post_val, paired_df$base_val, paired = TRUE)
    
    tibble(
      cohort = "ZUMA1_Paired", feature = feature,
      N_pairs = nrow(paired_df),
      median_Baseline = median(paired_df$base_val),
      median_Post = median(paired_df$post_val),
      p_value = test$p.value
    )
  }) %>%
    mutate(
      padj_BH = p.adjust(p_value, method = "BH"),
      label_padj = paste0("FDR = ", signif(padj_BH, 2))
    ) %>%
    arrange(p_value)
  
  write_csv(stats, "results/cellular_composition/ZUMA1_Paired_Stats.csv")
  return(stats)
}

message("Calculando estatísticas pareadas (ZUMA-1)...")
# NOTA: Troque "patient_id" pelo nome correto da coluna identificadora dos pacientes no ZUMA-1
stat_z1_paired <- run_paired_stats(mcp_z1_base, mcp_z1_post, mcp_cols, id_col = "patient_id")

## =============================================================================
## 7. Funções de Plotagem
## =============================================================================

# 7.1 Violinos para Análises Não-Pareadas
plot_unpaired_violin <- function(df, stats_df, feature_cols, cohort_name) {
  valid_features <- intersect(feature_cols, colnames(df))
  
  if (!"response_original" %in% colnames(df)) df$response_original <- "Missing"
  
  plot_df <- df %>%
    select(sample, group = response_group, response_original, all_of(valid_features)) %>%
    filter(group != "Unknown", group != "Missing", !is.na(group)) %>%
    pivot_longer(cols = all_of(valid_features), names_to = "feature", values_to = "score") %>%
    drop_na(score) %>%
    mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original))
  
  global_range <- diff(range(plot_df$score, na.rm = TRUE))
  
  label_df <- plot_df %>%
    group_by(feature) %>%
    summarise(y.position = max(score) + (0.05 * global_range), .groups = "drop") %>%
    left_join(stats_df, by = "feature") %>%
    mutate(group1 = group_1, group2 = group_2) %>% 
    drop_na(padj_BH)
  
  p <- ggplot(plot_df, aes(x = group, y = score)) +
    geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", trim = FALSE, scale = "width") +
    geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
    geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85, stroke = 0.5) +
    facet_wrap(~feature, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = pal_binaria, name = "Resposta Global:") +
    scale_color_manual(values = pal_granular, name = "Classificação Clínica:") +
    theme_classic(base_size = 14) +
    labs(title = paste(cohort_name, "(MCP-counter)"), x = NULL, y = "Estimativa de Abundância Celular") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      legend.position = "bottom", legend.box = "vertical",
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  if (nrow(label_df) > 0) {
    p <- p + stat_pvalue_manual(label_df, label = "label_padj", y.position = "y.position", 
                                tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  ggsave(paste0("results/cellular_composition/", cohort_name, "_Violin.pdf"), p, width = 12, height = 10)
}

# 7.2 Gráfico de Linhas Pareadas (Slope Graph) para ZUMA-1 Post vs Base
plot_paired_lines <- function(df_base, df_post, stats_df, feature_cols, id_col = "patient_id") {
  valid_features <- intersect(feature_cols, intersect(colnames(df_base), colnames(df_post)))
  
  # Preparar dados em formato longo vinculados pelo paciente
  df_b <- df_base %>% select(all_of(id_col), group = response_group, all_of(valid_features)) %>% mutate(Timepoint = "Baseline")
  df_p <- df_post %>% select(all_of(id_col), group = response_group, all_of(valid_features)) %>% mutate(Timepoint = "Post-CART")
  
  plot_df <- bind_rows(df_b, df_p) %>%
    filter(group != "Unknown", group != "Missing", !is.na(group)) %>%
    pivot_longer(cols = all_of(valid_features), names_to = "feature", values_to = "score") %>%
    drop_na(score) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("Baseline", "Post-CART")))
  
  # Garantir que só plota pacientes com os dois pontos de tempo
  complete_patients <- plot_df %>% count(!!sym(id_col), feature) %>% filter(n == 2) %>% pull(!!sym(id_col))
  plot_df <- plot_df %>% filter(!!sym(id_col) %in% complete_patients)
  
  global_range <- diff(range(plot_df$score, na.rm = TRUE))
  
  label_df <- plot_df %>%
    group_by(feature) %>%
    summarise(y.position = max(score) + (0.05 * global_range), .groups = "drop") %>%
    left_join(stats_df, by = "feature") %>%
    mutate(group1 = "Baseline", group2 = "Post-CART") %>% 
    drop_na(padj_BH)
  
  p <- ggplot(plot_df, aes(x = Timepoint, y = score)) +
    geom_line(aes(group = !!sym(id_col), color = group), alpha = 0.4, linewidth = 0.8) +
    geom_boxplot(aes(fill = Timepoint), width = 0.2, alpha = 0.6, outlier.shape = NA, color = "black") +
    geom_point(aes(fill = Timepoint), shape = 21, size = 2, color = "white", alpha = 0.9) +
    facet_wrap(~feature, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = pal_tempo, name = "Tempo:") +
    scale_color_manual(values = pal_binaria, name = "Resposta Clínica:") +
    theme_classic(base_size = 14) +
    labs(title = "ZUMA-1: Pareado (Post-CART vs Baseline)", x = NULL, y = "Estimativa de Abundância Celular (MCP-counter)") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = "gray80"),
      legend.position = "bottom", legend.box = "vertical",
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  if (nrow(label_df) > 0) {
    p <- p + stat_pvalue_manual(label_df, label = "label_padj", y.position = "y.position", 
                                tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  ggsave("results/cellular_composition/ZUMA1_Paired_SlopeGraph.pdf", p, width = 12, height = 10)
}

## =============================================================================
## 8. Executar Plotagens
## =============================================================================
message("\nGerando figuras...")
plot_unpaired_violin(mcp_z1_base, stat_z1_base, mcp_cols, "ZUMA1_Baseline")
plot_unpaired_violin(mcp_z7_axi,  stat_z7_axi,  mcp_cols, "ZUMA7_AxiCel")
plot_unpaired_violin(mcp_z7_soc,  stat_z7_soc,  mcp_cols, "ZUMA7_SOC")

# NOTA: Troque "patient_id" caso sua matriz use outro nome para o ID do paciente
plot_paired_lines(mcp_z1_base, mcp_z1_post, stat_z1_paired, mcp_cols, id_col = "patient_id")

cat("\n==========================================================\n")
cat("MCP-counter e Análises Concluídas com Sucesso!\n")
cat("==========================================================\n")