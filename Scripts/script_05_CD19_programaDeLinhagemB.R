################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 05 - Disponibilidade do Alvo (CD19) e Programa de Linhagem B
# Objetivo:
#   1. Avaliar CD19 e calcular o Score do Programa B (Z-score)
#   2. Estatística e Figuras de Alta Qualidade (Violin + Boxplot)
#   3. Correlação com a Composição Celular
#   4. Análise de Sobrevida Livre de Eventos (EFS) no ZUMA-7
################################################################################

## =============================================================================
## 1. Pacotes
## =============================================================================
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(survival)
library(survminer)

setwd("/media/jean/OneDrive/PesquisaCientifica/BigData/Projetos/cart_cd19_dlbcl_062026")

options(stringsAsFactors = FALSE)
dir.create("results/target_bcell_program", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados Purificados (Fases 02 e 03)
## =============================================================================
message("Carregando matrizes e composição celular...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")

# Carregar composição celular gerada no Script 03 (NOME CORRIGIDO)
comp_z1 <- read_csv("results/cellular_composition/ZUMA1_Baseline_Manual_Scores.csv", show_col_types = FALSE)
comp_z7 <- read_csv("results/cellular_composition/ZUMA7_AxiCel_mRNA_Manual_Scores.csv", show_col_types = FALSE)

## =============================================================================
## 3. Paletas de Cores Estilo Publicação
## =============================================================================
pal_binaria <- c(
  "CR" = "#4DBBD5FF", "Non-CR" = "#E64B35FF",
  "Ongoing" = "#00A087FF", "Others" = "#3C5488FF"
)

pal_granular <- c(
  "CR" = "#4DBBD5FF", "PR" = "#F39B7FFF", "SD" = "#7E6148FF", "PD" = "#E64B35FF",
  "Ongoing Response" = "#00A087FF", "Relapsed" = "#8491B4FF", "Nonresponders" = "#3C5488FF", "Missing" = "gray80"
)

heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

## =============================================================================
## 4. Definição e Cálculo do Score de Programa B
## =============================================================================
b_program_genes <- c("CD19", "MS4A1", "CD22", "CD79A", "CD79B", "PAX5", "POU2AF1", "BTK", "BLK", "FCRL2", "ST6GAL1", "BLNK", "CD74", "HLA-DRA")

calculate_b_program_score <- function(expr_matrix) {
  genes_present <- intersect(b_program_genes, rownames(expr_matrix))
  if (length(genes_present) < 2) stop("Poucos genes B disponíveis.")
  
  # Extrair, transformar em Z-score (scale) e calcular a média por paciente
  expr_sub <- expr_matrix[genes_present, , drop = FALSE]
  expr_scaled <- t(scale(t(expr_sub)))
  score <- colMeans(expr_scaled, na.rm = TRUE)
  return(score)
}

build_target_df <- function(dataset_list, cohort_name, composition_df) {
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  
  genes_present <- intersect(b_program_genes, rownames(expr))
  
  df <- t(expr[genes_present, , drop = FALSE]) %>%
    as.data.frame() %>%
    rownames_to_column("sample")
  
  df$B_program_score <- calculate_b_program_score(expr)
  
  df <- df %>%
    left_join(clin, by = "sample") %>%
    left_join(composition_df %>% select(sample, B_lineage_proxy, immune_proxy, tumor_immune_ratio_proxy), by = "sample") %>%
    mutate(cohort = cohort_name)
  
  return(df)
}

df_z1 <- build_target_df(z1_base, "ZUMA1_Baseline", comp_z1)
df_z7 <- build_target_df(z7_mrna, "ZUMA7_mRNA", comp_z7)

## =============================================================================
## 5. Estatística Gene a Gene (Wilcoxon + FDR)
## =============================================================================
run_target_stats <- function(df, feature_cols, cohort_name) {
  valid_features <- intersect(feature_cols, colnames(df))
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>% select(group = response_group, value = all_of(feature)) %>% filter(group != "Unknown") %>% drop_na()
    tmp$group <- as.factor(tmp$group)
    if (nlevels(tmp$group) != 2) return(NULL)
    
    test <- wilcox.test(value ~ group, data = tmp)
    tibble(
      cohort = cohort_name, feature = feature,
      group_1 = levels(tmp$group)[1], group_2 = levels(tmp$group)[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(padj_BH = p.adjust(p_value, method = "BH"), label_padj = paste0("FDR = ", signif(padj_BH, 2)))
  
  write_csv(stats, paste0("results/target_bcell_program/", cohort_name, "_Bprogram_Stats.csv"))
  return(stats)
}

#features_to_test <- c("CD19", "MS4A1", "CD79A", "PAX5", "B_program_score")
features_to_test <- c(intersect(b_program_genes, colnames(df_z1)), "B_program_score")# Intersecta a nossa lista dos sonhos com os genes que realmente existem na matriz do ZUMA-7 e junta com o Score
stat_z1 <- run_target_stats(df_z1, features_to_test, "ZUMA1_Baseline")
features_to_test <- c(intersect(b_program_genes, colnames(df_z7)), "B_program_score")# Intersecta a nossa lista dos sonhos com os genes que realmente existem na matriz do ZUMA-7 e junta com o Score
stat_z7 <- run_target_stats(df_z7, features_to_test, "ZUMA7_mRNA")

## =============================================================================
## 6. Violin Plots Híbridos Individuais (Prontos para Montagem de Painel)
## =============================================================================
message("Gerando Figuras Individuais para cada Gene/Score...")

# Criar uma subpasta para manter a organização
dir.create("results/target_bcell_program/single_plots", showWarnings = FALSE, recursive = TRUE)

plot_single_target_violin <- function(df, stats_df, feature_name, cohort_name) {
  
  # Preparar dados focados apenas na feature atual
  plot_df <- df %>%
    select(sample, group = response_group, response_original, score = all_of(feature_name)) %>%
    filter(group != "Unknown") %>%
    mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original)) %>% 
    drop_na(score)
  
  # Extrair a estatística específica
  stat_row <- stats_df %>% 
    filter(feature == feature_name) %>% 
    mutate(group1 = group_1, group2 = group_2) %>% 
    drop_na(padj_BH)
  
  # Definir teto global para a barra de P-value não invadir o gráfico
  y_pos <- max(plot_df$score, na.rm = TRUE) + 0.10 * diff(range(plot_df$score, na.rm = TRUE))
  
  p <- ggplot(plot_df, aes(x = group, y = score)) +
    geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", trim = FALSE, scale = "width") +
    geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
    geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85, stroke = 0.5) +
    
    scale_fill_manual(values = pal_binaria) +
    scale_color_manual(values = pal_granular) +
    
    theme_classic(base_size = 14) +
    labs(
      title = paste0(cohort_name, "\n", feature_name), 
      x = NULL, 
      y = "Log2 Expression / Z-score"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      legend.position = "none", # Removendo a legenda para montagem limpa
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  if(nrow(stat_row) > 0) {
    p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, 
                                tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  # Limpeza do nome do arquivo
  safe_filename <- gsub("[^A-Za-z0-9]", "_", feature_name)
  
  ggsave(paste0("results/target_bcell_program/single_plots/", cohort_name, "_", safe_filename, "_Violin.pdf"), p, width = 4.5, height = 5)
  ggsave(paste0("results/target_bcell_program/single_plots/", cohort_name, "_", safe_filename, "_Violin.png"), p, width = 4.5, height = 5, dpi = 300)
}

# 6.1 Criar e Exportar a Legenda Isolada
p_legend <- ggplot(df_z7, aes(x = response_group, y = B_program_score, fill = response_group)) +
  geom_violin() +
  geom_jitter(aes(color = response_original)) +
  scale_fill_manual(values = pal_binaria, name = "Resposta Global:") +
  scale_color_manual(values = pal_granular, name = "Classificação Clínica:") +
  theme_classic() + theme(legend.position = "bottom", legend.box = "vertical")

legenda_isolada <- cowplot::get_legend(p_legend)
ggsave("results/target_bcell_program/single_plots/00_Legenda_Isolada.pdf", cowplot::ggdraw(legenda_isolada), width = 6, height = 2)

# 6.2 Executar a plotagem individual para as duas coortes

message("Exportando gráficos individuais do ZUMA-1...")
# Cria uma lista apenas com os genes que sobreviveram no ZUMA-1
features_z1 <- c(intersect(b_program_genes, colnames(df_z1)), "B_program_score")
walk(features_z1, ~plot_single_target_violin(df_z1, stat_z1, .x, "ZUMA1_Baseline"))

message("Exportando gráficos individuais do ZUMA-7...")
# Cria uma lista apenas com os genes que sobreviveram no ZUMA-7
features_z7 <- c(intersect(b_program_genes, colnames(df_z7)), "B_program_score")
walk(features_z7, ~plot_single_target_violin(df_z7, stat_z7, .x, "ZUMA7_mRNA"))
## =============================================================================
## 7. Análise de Sobrevida (EFS) - ZUMA-7 Axi-cel
## =============================================================================
message("Rodando Análise de Sobrevida de Kaplan-Meier (ZUMA-7)...")

# Mapeando os dados originais de sobrevida com os nomes exatos da matriz
surv_df <- df_z7 %>%
  mutate(
    efs_months = as.numeric(event_free_survival_months_ch1),
    efs_event = as.numeric(event_free_survival_event_ch1)
  ) %>%
  filter(!is.na(efs_months), !is.na(efs_event)) %>%
  mutate(
    CD19_Status = ifelse(CD19 >= median(CD19, na.rm = TRUE), "CD19 High", "CD19 Low"),
    B_Program_Status = ifelse(B_program_score >= median(B_program_score, na.rm = TRUE), "B-Program High", "B-Program Low"),
    CD19_Status = factor(CD19_Status, levels = c("CD19 Low", "CD19 High")),
    B_Program_Status = factor(B_Program_Status, levels = c("B-Program Low", "B-Program High"))
  )

plot_km <- function(df, group_col, title, filename) {
  
  # 1. Contorno do Bug: Criar uma coluna fixa chamada 'km_group'
  plot_df <- df %>%
    filter(!is.na(.data[[group_col]])) %>%
    mutate(km_group = .data[[group_col]])
  
  # 2. Rodar o modelo usando a coluna fixa (sem as.formula)
  fit <- survfit(Surv(efs_months, efs_event) ~ km_group, data = plot_df)
  
  # 3. Gerar o gráfico
  p <- ggsurvplot(
    fit, data = plot_df, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
    palette = c("#E64B35FF", "#4DBBD5FF"), # Vermelho (Low) e Azul (High)
    title = title, xlab = "Event-Free Survival (Months)", ylab = "Probability",
    legend.title = "Expression:", ggtheme = theme_classic()
  )
  
  # 4. Exportar o PDF
  pdf(paste0("results/target_bcell_program/", filename, ".pdf"), width = 7, height = 7, onefile = FALSE)
  print(p)
  dev.off()
}

# Agora as chamadas funcionarão perfeitamente:
plot_km(surv_df, "CD19_Status", "ZUMA-7: EFS by Baseline CD19", "ZUMA7_EFS_CD19")
plot_km(surv_df, "B_Program_Status", "ZUMA-7: EFS by B-cell Program", "ZUMA7_EFS_BProgram")

## =============================================================================
## 8. Modelos de Cox (Razão de Risco) e Forest Plot
## =============================================================================
message("Calculando Modelos de Cox e gerando Forest Plot...")

cox_cd19 <- coxph(Surv(efs_months, efs_event) ~ CD19, data = surv_df)
cox_bprog <- coxph(Surv(efs_months, efs_event) ~ B_program_score, data = surv_df)

# Juntando os resultados em uma tabela
cox_res <- bind_rows(
  broom::tidy(cox_cd19, exponentiate = TRUE, conf.int = TRUE) %>% mutate(Model = "CD19 Expression (Contínuo)"),
  broom::tidy(cox_bprog, exponentiate = TRUE, conf.int = TRUE) %>% mutate(Model = "B-Program (Z-score)")
)

write_csv(cox_res, "results/target_bcell_program/ZUMA7_Cox_Survival_Models.csv")

# 8.1 Construção do Forest Plot com ggplot2
plot_forest <- cox_res %>%
  mutate(
    # Formatando os números para ficarem bonitos no gráfico (ex: "0.85 (0.70 - 0.95)")
    HR_Label = sprintf("%.2f (%.2f - %.2f)", estimate, conf.low, conf.high),
    P_Label = ifelse(p.value < 0.001, "p < 0.001", sprintf("p = %.3f", p.value)),
    Significance = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ "")
  ) %>%
  ggplot(aes(y = fct_rev(Model))) +
  
  # 1. A linha de referência pontilhada no 1.0 (Risco Neutro)
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # 2. Barras de erro (Intervalo de Confiança de 95%)
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, x = estimate), height = 0.15, linewidth = 1) +
  
  # 3. O quadrado central (Hazard Ratio exato)
  geom_point(aes(x = estimate), shape = 15, size = 5) +
  
  # 4. Texto do HR flutuando em cima da barra
  geom_text(aes(x = estimate, label = HR_Label), vjust = -1.5, size = 4.5, fontface = "bold") +
  
  # 5. Texto do P-value alinhado à direita
  geom_text(aes(x = max(conf.high) * 1.15, label = paste0(P_Label, Significance)), hjust = 0, size = 4.5, fontface = "italic") +
  
  # 6. Escala logarítmica (Padrão para visualização de Hazard Ratio)
  scale_x_continuous(trans = "log10") +
  
  theme_classic(base_size = 14) +
  labs(
    title = "ZUMA-7: Risco de Progressão/Morte (Hazard Ratios)",
    x = "Hazard Ratio (Escala Logarítmica)",
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(face = "bold", color = "black", size = 12),
    axis.line.y = element_blank(), # Remove a linha vertical do eixo Y para ficar mais limpo
    axis.ticks.y = element_blank()
  ) +
  # Evita que o texto do p-valor seja cortado da imagem
  coord_cartesian(clip = "off", xlim = c(min(cox_res$conf.low)*0.85, max(cox_res$conf.high)*1.4))

ggsave("results/target_bcell_program/ZUMA7_Cox_Forest_Plot.pdf", plot_forest, width = 10, height = 4)
ggsave("results/target_bcell_program/ZUMA7_Cox_Forest_Plot.png", plot_forest, width = 10, height = 4, dpi = 300)

