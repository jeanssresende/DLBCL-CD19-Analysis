################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 05 - Disponibilidade do Alvo e Genes da Linhagem B (Individualizados)
# Objetivo:
#   1. Avaliar genes originais validadas por coorte (ZUMA-1 Base/Post, ZUMA-7 Axi-cel/SOC)
#   2. Calcular estatística (Wilcoxon Não-Pareado e Pareado) com correção FDR restrita
#   3. Gerar Figuras de Alta Qualidade (Violin + Boxplot + Slope Graphs Pareados)
################################################################################

## =============================================================================
## 1. Pacotes e Diretórios
## =============================================================================
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

setwd("/media/jean/OneDrive/PesquisaCientifica/BigData/Projetos/cart_cd19_dlbcl_062026")
options(stringsAsFactors = FALSE)
dir.create("results/target_bcell_individual", showWarnings = FALSE, recursive = TRUE)
dir.create("results/target_bcell_individual/plots", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados Purificados
## =============================================================================
message("Carregando matrizes de expressão...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z1_post <- readRDS("data/filtered/ZUMA1_Post_Clean.rds") # Adicionado ZUMA-1 Post
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
z7_soc  <- readRDS("data/filtered/ZUMA7_SOC_mRNA_Clean.rds") 

# CORREÇÃO: Criar a coluna 'response_group' no SOC a partir da coluna existente
if("ongoing_2grps_ch1" %in% colnames(z7_soc$clin) && !"response_group" %in% colnames(z7_soc$clin)) {
  z7_soc$clin$response_group <- z7_soc$clin$ongoing_2grps_ch1
}

## =============================================================================
## 3. Paletas de Cores Estilo Publicação
## =============================================================================
pal_binaria <- c(
  "CR" = "#4DBBD5", "Non-CR" = "#E64B35",
  "Ongoing" = "#00A087", "Others" = "#3C5488"
)

pal_granular <- c(
  "CR" = "#4DBBD5", "PR" = "#F39B7F", "SD" = "#7E6148", "PD" = "#E64B35",
  "Ongoing Response" = "#00A087", "Relapsed" = "#8491B4", "Nonresponders" = "#3C5488", 
  "Missing" = "gray80"
)

# Paleta exata para a Evolução do CR Pareado (Azul Claro -> Azul Safira)
pal_tempo_cr <- c("Baseline" = "#4DBBD5", "Post-CART" = "#0072B5")

## =============================================================================
## 4. Definição das Listas de Genes Específicas por Estudo
## =============================================================================
# Lista validada pelos autores do ZUMA-1
z1_genes <- c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "CD22", "ST6GAL1", "POU2AF1", "BTK")

# Lista validada pelos autores do ZUMA-7 (B cell GES)
z7_genes <- c("BLK", "CD19", "MS4A1", "TNFRSF17", "FCRL2", "FAM30A", "PNOC", "SPIB", "TCL1A")

## =============================================================================
## 5. Função de Extração de Dados
## =============================================================================
build_target_df <- function(dataset_list, cohort_name, target_genes) {
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  
  # Garantir que só procuramos genes que realmente existem nesta matriz
  genes_present <- intersect(target_genes, rownames(expr))
  message(paste(cohort_name, "- Genes encontrados:", length(genes_present), "de", length(target_genes)))
  
  df <- t(expr[genes_present, , drop = FALSE]) %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(clin, by = "sample") %>%
    mutate(cohort = cohort_name)
  
  # Extração do ID do paciente (Crucial para a análise pareada do ZUMA-1)
  if("title" %in% colnames(df)) {
    df <- df %>% mutate(patient_id = str_extract(title, "Patient \\d+"))
  }
  
  return(list(data = df, genes = genes_present))
}

data_z1_base <- build_target_df(z1_base, "ZUMA1_Baseline", z1_genes)
data_z1_post <- build_target_df(z1_post, "ZUMA1_Post", z1_genes)
data_z7_axi  <- build_target_df(z7_mrna, "ZUMA7_AxiCel", z7_genes)
data_soc     <- build_target_df(z7_soc, "ZUMA7_SOC", z7_genes)

## =============================================================================
## 6. Estatística Gene a Gene (Wilcoxon + FDR focado)
## =============================================================================

# 6.1 Função para Análise Não-Pareada
run_focused_stats <- function(df_obj, cohort_name) {
  df <- df_obj$data
  valid_features <- df_obj$genes
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>% 
      select(group = response_group, value = all_of(feature)) %>% 
      filter(group != "Unknown" & group != "Missing" & !is.na(group)) %>% 
      drop_na(value)
    
    tmp$group <- as.factor(tmp$group)
    if (nlevels(tmp$group) != 2) return(NULL)
    
    test <- wilcox.test(value ~ group, data = tmp)
    tibble(
      cohort = cohort_name, feature = feature,
      group_1 = levels(tmp$group)[1], group_2 = levels(tmp$group)[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(
      padj_BH = p.adjust(p_value, method = "BH"), 
      label_padj = paste0("FDR = ", signif(padj_BH, 2))
    ) %>%
    arrange(p_value)
  
  write_csv(stats, paste0("results/target_bcell_individual/Stats_", cohort_name, "_Unpaired.csv"))
  return(stats)
}

# 6.2 Função para Análise Pareada (ZUMA-1 Post vs Baseline)
run_paired_focused_stats <- function(df_base_obj, df_post_obj, cohort_name, id_col = "patient_id") {
  df_base <- df_base_obj$data
  df_post <- df_post_obj$data
  valid_features <- intersect(df_base_obj$genes, df_post_obj$genes)
  
  stats <- map_dfr(valid_features, function(feature) {
    b_tmp <- df_base %>% select(all_of(id_col), base_val = all_of(feature))
    p_tmp <- df_post %>% select(all_of(id_col), post_val = all_of(feature))
    
    # Junção restrita a pacientes com amostras nos dois tempos
    paired_df <- inner_join(b_tmp, p_tmp, by = id_col) %>% drop_na()
    
    if (nrow(paired_df) < 3) return(NULL)
    
    # Teste de Wilcoxon Pareado
    test <- wilcox.test(paired_df$post_val, paired_df$base_val, paired = TRUE)
    
    tibble(
      cohort = cohort_name, feature = feature,
      N_pairs = nrow(paired_df),
      p_value = test$p.value
    )
  }) %>%
    mutate(
      padj_BH = p.adjust(p_value, method = "BH"), 
      label_padj = paste0("FDR = ", signif(padj_BH, 2))
    ) %>%
    arrange(p_value)
  
  write_csv(stats, paste0("results/target_bcell_individual/Stats_", cohort_name, "_Paired.csv"))
  return(stats)
}

message("\nCalculando estatísticas Não-Pareadas...")
stat_z1_base <- run_focused_stats(data_z1_base, "ZUMA1_Baseline")
stat_z7_axi  <- run_focused_stats(data_z7_axi, "ZUMA7_AxiCel")
stat_soc     <- run_focused_stats(data_soc, "ZUMA7_SOC")

message("Calculando estatísticas Pareadas...")
stat_z1_paired <- run_paired_focused_stats(data_z1_base, data_z1_post, "ZUMA1_Paired")

## =============================================================================
## 7. Plotagem de Gráficos Individuais
## =============================================================================

# 7.1 Função para Violinos (Análises Não-Pareadas)
plot_single_target_violin <- function(df_obj, stats_df, feature_name, cohort_name) {
  df <- df_obj$data
  
  if (!"response_original" %in% colnames(df)) {
    df$response_original <- "Missing"
  }
  
  plot_df <- df %>%
    select(sample, group = response_group, response_original, score = all_of(feature_name)) %>%
    filter(group != "Unknown" & group != "Missing" & !is.na(group)) %>%
    mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original)) %>% 
    drop_na(score)
  
  stat_row <- stats_df %>% filter(feature == feature_name) %>% 
    mutate(group1 = group_1, group2 = group_2) %>% drop_na(padj_BH)
  
  y_pos <- max(plot_df$score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE))
  
  p <- ggplot(plot_df, aes(x = group, y = score)) +
    geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", trim = FALSE, scale = "width") +
    geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
    geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85, stroke = 0.5) +
    scale_fill_manual(values = pal_binaria) +
    scale_color_manual(values = pal_granular) +
    theme_classic(base_size = 14) +
    labs(title = paste0(cohort_name, ": ", feature_name), x = NULL, y = expression(Log[2]~Expression)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "none",
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  if(nrow(stat_row) > 0) {
    p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  safe_filename <- gsub("[^A-Za-z0-9]", "_", feature_name)
  ggsave(paste0("results/target_bcell_individual/plots/", cohort_name, "_", safe_filename, "_Violin.pdf"), p, width = 4.5, height = 5)
}

# 7.2 Função para Gráfico de Linhas Pareadas (Slope Graph)
plot_paired_lines_genes <- function(df_base_obj, df_post_obj, stats_df, feature_name, cohort_name, id_col = "patient_id") {
  
  df_base <- df_base_obj$data
  df_post <- df_post_obj$data
  
  b_tmp <- df_base %>% select(all_of(id_col), score = all_of(feature_name)) %>% mutate(Timepoint = "Baseline")
  p_tmp <- df_post %>% select(all_of(id_col), score = all_of(feature_name)) %>% mutate(Timepoint = "Post-CART")
  
  plot_df <- bind_rows(b_tmp, p_tmp) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("Baseline", "Post-CART"))) %>%
    drop_na(score)
  
  # Garantir que só plota pacientes com os dois pontos de tempo
  complete_patients <- plot_df %>% count(!!sym(id_col)) %>% filter(n == 2) %>% pull(!!sym(id_col))
  plot_df <- plot_df %>% filter(!!sym(id_col) %in% complete_patients)
  
  stat_row <- stats_df %>% filter(feature == feature_name) %>%
    mutate(group1 = "Baseline", group2 = "Post-CART") %>% drop_na(padj_BH)
  
  y_pos <- max(plot_df$score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE))
  
  p <- ggplot(plot_df, aes(x = Timepoint, y = score)) +
    geom_line(aes(group = !!sym(id_col)), color = "gray60", alpha = 0.5, linewidth = 0.8) +
    geom_boxplot(aes(fill = Timepoint), width = 0.2, alpha = 0.7, outlier.shape = NA, color = "black") +
    geom_point(aes(fill = Timepoint), shape = 21, size = 2.5, color = "white", alpha = 0.9) +
    
    scale_fill_manual(values = pal_tempo_cr) +
    
    theme_classic(base_size = 14) +
    labs(title = paste0(cohort_name, ": ", feature_name), x = NULL, y = expression(Log[2]~Expression)) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "none",
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  if(nrow(stat_row) > 0) {
    p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  safe_filename <- gsub("[^A-Za-z0-9]", "_", feature_name)
  ggsave(paste0("results/target_bcell_individual/plots/", cohort_name, "_", safe_filename, "_Paired.pdf"), p, width = 4.5, height = 5)
}

## =============================================================================
## 8. Execução Final das Plotagens
## =============================================================================
message("\nExportando gráficos individuais (Violinos)...")
walk(data_z1_base$genes, ~plot_single_target_violin(data_z1_base, stat_z1_base, .x, "ZUMA1_Baseline"))
walk(data_z7_axi$genes,  ~plot_single_target_violin(data_z7_axi, stat_z7_axi, .x, "ZUMA7_AxiCel"))
walk(data_soc$genes,     ~plot_single_target_violin(data_soc, stat_soc, .x, "ZUMA7_SOC"))

message("Exportando gráficos pareados (Slope Graphs)...")
valid_paired_features <- intersect(data_z1_base$genes, data_z1_post$genes)
walk(valid_paired_features, ~plot_paired_lines_genes(data_z1_base, data_z1_post, stat_z1_paired, .x, "ZUMA1_Paired"))

message("\nAnálise concluída com sucesso! Verifique a pasta 'results/target_bcell_individual'.")
