################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 04 - Estados Funcionais do TME (Assinaturas IO360 no ZUMA-7)
# Objetivo:
#   1. Interrogar TODAS as assinaturas do painel NanoString IO360 sem viés
#   2. Comparar Ongoing vs Others (Wilcoxon + FDR)
#   3. Isolar as Top assinaturas mais discrepantes
#   4. Gerar figuras híbridas de alta qualidade (Violin + Boxplot)
################################################################################

## =============================================================================
## 1. Pacotes e Diretórios
## =============================================================================
library(tidyverse)
library(ggpubr)
library(patchwork)

options(stringsAsFactors = FALSE)
dir.create("results/io360_signatures", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados Purificados (Fase 02)
## =============================================================================
message("Carregando matriz de assinaturas IO360 do ZUMA-7...")

z7_sig <- readRDS("data/filtered/ZUMA7_AxiCel_Sigs_Clean.rds")
expr_sig <- z7_sig$expr
clin_sig <- z7_sig$clin %>% rownames_to_column("sample")

## =============================================================================
## 3. Paletas de Cores Estilo Publicação
## =============================================================================
pal_binaria <- c("Ongoing" = "#00A087FF", "Others" = "#3C5488FF")
pal_granular <- c("Ongoing Response" = "#00A087FF", "Relapsed" = "#8491B4FF", 
                  "Nonresponders" = "#3C5488FF", "Missing" = "gray80")

## =============================================================================
## 4. Estruturar o DataFrame Analítico (TODAS AS ASSINATURAS)
## =============================================================================
# Em vez de filtrar, pegamos os nomes de todas as assinaturas presentes na matriz
all_io360_features <- rownames(expr_sig)
message("Testando de forma imparcial todas as ", length(all_io360_features), " assinaturas disponíveis...")

df_io360 <- t(expr_sig) %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(clin_sig, by = "sample")

## =============================================================================
## 5. Estatística Global (Wilcoxon + FDR)
## =============================================================================
message("Calculando estatísticas (Ongoing vs Others)...")

run_signature_stats <- function(df, feature_cols) {
  
  stats <- map_dfr(feature_cols, function(feature) {
    tmp <- df %>% 
      select(group = response_group, value = all_of(feature)) %>% 
      filter(group != "Unknown") %>% 
      drop_na()
    
    tmp$group <- as.factor(tmp$group)
    if (nlevels(tmp$group) != 2) return(NULL)
    
    test <- wilcox.test(value ~ group, data = tmp)
    
    summary_group <- tmp %>%
      group_by(group) %>%
      summarise(median = median(value, na.rm = TRUE), .groups = "drop")
    
    tibble(
      feature = feature,
      group_1 = levels(tmp$group)[1], 
      group_2 = levels(tmp$group)[2],
      median_group_1 = summary_group$median[1],
      median_group_2 = summary_group$median[2],
      p_value = test$p.value
    )
  }) %>%
    mutate(
      padj_BH = p.adjust(p_value, method = "BH"), 
      label_padj = paste0("FDR = ", signif(padj_BH, 2))
    ) %>%
    arrange(p_value) # Ordena pela força do sinal estatístico
  
  write_csv(stats, "results/io360_signatures/ZUMA7_IO360_All_Stats.csv")
  return(stats)
}

stat_io360 <- run_signature_stats(df_io360, all_io360_features)

# Imprimir no console as 10 assinaturas mais fortes para interpretação imediata
message("\n--- TOP 10 Assinaturas Mais Discrepantes ---")
print(stat_io360 %>% select(feature, p_value, padj_BH) %>% head(10))

## =============================================================================
## 6. Violin Plots Híbridos Individuais (Prontos para Montagem de Painel)
## =============================================================================
message("\nGerando Figuras Individuais para cada Assinatura...")

# Criar uma subpasta para não poluir o diretório principal
dir.create("results/io360_signatures/single_plots", showWarnings = FALSE, recursive = TRUE)

plot_single_signature <- function(df, stats_df, feature_name) {
  
  # Preparar os dados apenas para a assinatura específica
  plot_df <- df %>%
    select(sample, group = response_group, response_original, score = all_of(feature_name)) %>%
    filter(group != "Unknown") %>%
    mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original)) %>% 
    drop_na(score)
  
  # Extrair a estatística específica dessa assinatura
  stat_row <- stats_df %>% 
    filter(feature == feature_name) %>% 
    mutate(group1 = group_1, group2 = group_2) %>% 
    drop_na(padj_BH)
  
  # Calcular a altura da barra do P-value
  y_pos <- max(plot_df$score) + 0.10 * diff(range(plot_df$score))
  
  p <- ggplot(plot_df, aes(x = group, y = score)) +
    
    geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", trim = FALSE, scale = "width") +
    geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
    geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85, stroke = 0.5) +
    
    scale_fill_manual(values = pal_binaria, name = "Resposta Global:") +
    scale_color_manual(values = pal_granular, name = "Classificação Clínica:") +
    
    theme_classic(base_size = 14) +
    labs(
      title = feature_name, 
      x = NULL, 
      y = "Score da Assinatura (Log2)"
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "none", # Ocultamos a legenda nos individuais para o gráfico ficar limpo
      axis.text.x = element_text(face = "bold", color = "black")
    )
  
  # Adicionar a barra de P-value apenas se a assinatura foi testada e tem resultado
  if(nrow(stat_row) > 0) {
    p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, 
                                tip.length = 0.015, size = 3.5, fontface = "bold")
  }
  
  # Limpar o nome da feature para não quebrar o salvamento do arquivo (ex: tira espaços e pontos)
  safe_filename <- gsub("[^A-Za-z0-9]", "_", feature_name)
  
  # Salvar os gráficos menores e com proporções otimizadas para painéis individuais
  ggsave(paste0("results/io360_signatures/single_plots/", safe_filename, "_Violin.pdf"), p, width = 4.5, height = 5)
  ggsave(paste0("results/io360_signatures/single_plots/", safe_filename, "_Violin.png"), p, width = 4.5, height = 5, dpi = 300)
}

# Criar uma legenda separada em um arquivo único para você colar na sua montagem final
p_legend <- ggplot(df_io360, aes(x = response_group, y = 1, fill = response_group)) +
  geom_violin() +
  geom_jitter(aes(color = response_original)) +
  scale_fill_manual(values = pal_binaria, name = "Resposta Global:") +
  scale_color_manual(values = pal_granular, name = "Classificação Clínica:") +
  theme_classic() + theme(legend.position = "bottom", legend.box = "vertical")
legenda_isolada <- cowplot::get_legend(p_legend)
ggsave("results/io360_signatures/single_plots/00_Legenda_Isolada.pdf", cowplot::ggdraw(legenda_isolada), width = 6, height = 2)

# Rodar um loop (walk) para gerar e salvar um gráfico por vez para todas as assinaturas
message("Exportando gráficos individuais (isso pode levar alguns segundos)...")
walk(all_io360_features, ~plot_single_signature(df_io360, stat_io360, .x))

cat("\n==========================================================\n")
cat("Análise Funcional do TME (IO360) Concluída!\n")
cat("Gráficos soltos e legenda salvos em: 'results/io360_signatures/single_plots/'\n")
cat("==========================================================\n")
