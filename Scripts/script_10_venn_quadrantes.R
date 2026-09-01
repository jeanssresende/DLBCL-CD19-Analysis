################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 10 - Diagrama de Venn e Correlação Transcriptômica Global (Quadrant Plot)
#
# Objetivo:
#   1. Gerar Diagrama de Venn focando na interseção direta ZUMA-1 vs ZUMA-7
#   2. Unir todos os genes das duas coortes (independente do p-valor)
#   3. Gerar um Gráfico de Quadrantes (LogFC Scatter Plot) de alta densidade
#   4. Calcular a correlação global de Pearson entre as coortes
################################################################################

# ==============================================================================
# 1. Pacotes e Diretórios
# ==============================================================================
library(tidyverse)
library(ggVennDiagram)
library(ggrepel)
library(patchwork)
# install.packages("ggpubr") # Descomente se não tiver o ggpubr instalado
library(ggpubr) 

out_dir <- "results/comparative_transcriptomics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Carregar Tabelas de DEA
# ==============================================================================
message("Carregando resultados do Limma...")
dea_z1 <- read_csv("results/dea_gsea/tables/DEA_ZUMA1_Baseline_NonCR_vs_CR_RawP.csv", show_col_types = FALSE)
dea_z7 <- read_csv("results/dea_gsea/tables/DEA_ZUMA7_AxiCel_Others_vs_Ongoing_RawP.csv", show_col_types = FALSE)

# ==============================================================================
# 3. Diagrama de Venn (Apenas ZUMA-1 vs ZUMA-7)
# ==============================================================================
message("Gerando Diagrama de Venn...")

# Filtrar DEGs pelos cortes estritos clássicos
degs_z1 <- dea_z1 %>% filter(P.Value < 0.05, abs(logFC) > 0.5, !is.na(Feature)) %>% pull(Feature)
degs_z7 <- dea_z7 %>% filter(P.Value < 0.05, abs(logFC) > 0.5, !is.na(Feature)) %>% pull(Feature)

venn_list <- list(
  "ZUMA-1 (Discovery)" = degs_z1,
  "ZUMA-7 (Validation)" = degs_z7
)

p_venn <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 1.2) +
  scale_fill_gradient(low = "#F4FAFE", high = "#3C5488FF") +
  theme(legend.position = "none") +
  labs(title = "Intersection of Resistance DEGs") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))

ggsave(file.path(out_dir, "Venn_Z1_vs_Z7.pdf"), p_venn, width = 6, height = 5)
ggsave(file.path(out_dir, "Venn_Z1_vs_Z7.png"), p_venn, width = 6, height = 5, dpi = 300)

# ==============================================================================
# 4. Gráfico de Quadrantes: Correlação Global de LogFC (COM N DE GENES)
# ==============================================================================
message("Processando dados para o Gráfico de Quadrantes (LogFC Scatter Plot)...")

# Fazer a junção de TODOS os genes encontrados em ambas as matrizes
df_merged <- inner_join(
  dea_z1 %>% select(Feature, logFC_Z1 = logFC, P_Z1 = P.Value),
  dea_z7 %>% select(Feature, logFC_Z7 = logFC, P_Z7 = P.Value),
  by = "Feature"
) %>% filter(!is.na(logFC_Z1), !is.na(logFC_Z7))

# Classificar os genes pela sua significância biológica
df_merged <- df_merged %>%
  mutate(
    Sig_Status = case_when(
      P_Z1 < 0.05 & abs(logFC_Z1) > 0.5 & P_Z7 < 0.05 & abs(logFC_Z7) > 0.5 ~ "Significant in Both",
      P_Z1 < 0.05 & abs(logFC_Z1) > 0.5 ~ "Significant in ZUMA-1",
      P_Z7 < 0.05 & abs(logFC_Z7) > 0.5 ~ "Significant in ZUMA-7",
      TRUE ~ "Not Significant"
    ),
    Sig_Status = factor(Sig_Status, levels = c("Not Significant", "Significant in ZUMA-1", "Significant in ZUMA-7", "Significant in Both"))
  )

# Identificar os genes top (que estão significativos em ambos) para a label
top_genes <- df_merged %>% filter(Sig_Status == "Significant in Both")

# ------------------------------------------------------------------------------
# CÁLCULO AUTOMÁTICO DOS "N" PARA A LEGENDA
# ------------------------------------------------------------------------------
n_total   <- nrow(df_merged)
n_not_sig <- sum(df_merged$Sig_Status == "Not Significant")
n_z1      <- sum(df_merged$Sig_Status == "Significant in ZUMA-1")
n_z7      <- sum(df_merged$Sig_Status == "Significant in ZUMA-7")
n_both    <- sum(df_merged$Sig_Status == "Significant in Both")

message(paste("Gerando Gráfico de Dispersão Global com N Total =", n_total))

p_quadrant <- ggplot(df_merged, aes(x = logFC_Z1, y = logFC_Z7)) +
  # Linhas dos eixos (Quadrantes)
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  
  # Adicionar os pontos
  geom_point(aes(color = Sig_Status, alpha = Sig_Status, size = Sig_Status)) +
  
  # Estilizar as cores e ADICIONAR O "n" AOS RÓTULOS
  scale_color_manual(
    values = c(
      "Not Significant" = "grey80",
      "Significant in ZUMA-1" = "#4DBBD5FF",
      "Significant in ZUMA-7" = "#E64B35FF",
      "Significant in Both" = "#3C5488FF"
    ),
    labels = c(
      paste0("Not Significant (n=", n_not_sig, ")"),
      paste0("Significant in ZUMA-1 (n=", n_z1, ")"),
      paste0("Significant in ZUMA-7 (n=", n_z7, ")"),
      paste0("Significant in Both (n=", n_both, ")")
    )
  ) +
  scale_alpha_manual(values = c(0.3, 0.7, 0.7, 1.0), guide = "none") +
  scale_size_manual(values = c(0.8, 1.5, 1.5, 3.0), guide = "none") +
  
  # Linha de regressão linear
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.8) +
  
  # Estatística de correlação
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  
  # Nomes dos genes significativos em ambos
  geom_text_repel(data = top_genes, aes(label = Feature), 
                  size = 4, fontface = "bold", box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf) +
  
  theme_bw(base_size = 14) +
  labs(
    title = "Global Transcriptomic Concordance",
    # Adicionando o N Total no subtítulo
    subtitle = paste0("Comparing Log2 Fold Changes (Total shared genes mapped: n=", n_total, ")"),
    x = "ZUMA-1 Log2 Fold Change (Non-CR vs CR)",
    y = "ZUMA-7 Log2 Fold Change (Others vs Ongoing)",
    color = "Gene Significance"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(file.path(out_dir, "Quadrant_LogFC_Z1_vs_Z7.pdf"), p_quadrant, width = 9, height = 7)
ggsave(file.path(out_dir, "Quadrant_LogFC_Z1_vs_Z7.png"), p_quadrant, width = 9, height = 7, dpi = 300)

# ==============================================================================
# EXTRAÇÃO DA TABELA DE DADOS DOS QUADRANTES (CSV)
# ==============================================================================
message("Exportando a tabela de dados dos quadrantes...")

# Organizar a tabela final para salvamento
table_quadrant_genes <- df_merged %>%
  select(
    Gene = Feature,
    Sig_Status,
    Log2FC_ZUMA1 = logFC_Z1,
    PValue_ZUMA1 = P_Z1,
    Log2FC_ZUMA7 = logFC_Z7,
    PValue_ZUMA7 = P_Z7
  ) %>%
  # Ordenar para que os genes "Significant in Both" apareçam primeiro,
  # seguidos pelos do ZUMA-1, ZUMA-7 e, por fim, os não significativos.
  arrange(Sig_Status, PValue_ZUMA1)

# Salvar em formato CSV na pasta de resultados comparativos
write_csv(
  table_quadrant_genes, 
  file.path(out_dir, "Transcriptomic_Concordance_Genes_Table.csv")
)

message("Tabela salva com sucesso em: ", file.path(out_dir, "Transcriptomic_Concordance_Genes_Table.csv"))

# ==============================================================================
# 5. Figura Combinada
# ==============================================================================
fig_combined <- p_venn + p_quadrant + plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = 'bold', size = 16))

ggsave(file.path(out_dir, "Figure_Comparative_Transcriptomics.pdf"), fig_combined, width = 14, height = 6)
ggsave(file.path(out_dir, "Figure_Comparative_Transcriptomics.png"), fig_combined, width = 14, height = 6, dpi = 300)

message("\n==========================================================")
message("Análise concluída! Verifique a pasta: ", out_dir)
message("==========================================================")