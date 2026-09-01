################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 12 - Validação de Biomarcadores: Axi-cel vs Standard of Care (SOC)
#
# Objetivo:
#   1. Carregar dados das coortes Axi-cel e SOC do ensaio ZUMA-7.
#   2. Testar a Sobrevida Livre de Eventos (EFS) para os genes prioritários em ambos os braços.
#   3. Gerar gráficos Side-by-Side (Axi-cel vs SOC) para distinguir 
#      biomarcadores preditivos (específicos do CAR-T) de prognósticos (gerais).
################################################################################

# ==============================================================================
# 1. Pacotes e Diretórios
# ==============================================================================
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

out_dir <- "results/zuma7_soc_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "KM_SideBySide"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Definir a Lista de Genes Validados Manualmente
# ==============================================================================
message("Carregando lista de genes prioritários (Axi-cel vs SOC)...")

target_genes <- c("DUSP5","HLA-DOB","MS4A1","KLRD1","CDH2","IRF8","LY9","WNT11",
                  "CD45RA","CD19","TNFRSF17","NCAM1","CCL20","BTLA","IFI27","LRRC32")

# ==============================================================================
# 3. Carregar Dados ZUMA-7 (Axi-cel)
# ==============================================================================
message("Carregando dados do braço Axi-cel (ZUMA-7)...")
z7_axicel <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")

expr_axi <- z7_axicel$expr
clin_axi <- z7_axicel$clin %>%
  rownames_to_column("sample") %>%
  mutate(
    EFS_time = as.numeric(event_free_survival_months_ch1),
    EFS_event = as.numeric(event_free_survival_event_ch1)
  ) %>%
  filter(!is.na(EFS_time), !is.na(EFS_event))

valid_samples_axi <- intersect(clin_axi$sample, colnames(expr_axi))
clin_axi <- clin_axi %>% filter(sample %in% valid_samples_axi)
expr_axi <- expr_axi[, clin_axi$sample, drop = FALSE]

# ==============================================================================
# 4. Carregar Dados ZUMA-7 (Standard of Care - SOC)
# ==============================================================================
message("Carregando dados do braço Standard of Care (SOC)...")
# NOTA: Certifique-se de ter este arquivo preparado com a mesma estrutura do Axi-cel
z7_soc <- readRDS("data/filtered/ZUMA7_SOC_mRNA_Clean.rds") 

expr_soc <- z7_soc$expr
clin_soc <- z7_soc$clin %>%
  rownames_to_column("sample") %>%
  mutate(
    EFS_time = as.numeric(event_free_survival_months_ch1),
    EFS_event = as.numeric(event_free_survival_event_ch1)
  ) %>%
  filter(!is.na(EFS_time), !is.na(EFS_event))

valid_samples_soc <- intersect(clin_soc$sample, colnames(expr_soc))
clin_soc <- clin_soc %>% filter(sample %in% valid_samples_soc)
expr_soc <- expr_soc[, clin_soc$sample, drop = FALSE]

# ==============================================================================
# 5. Análise de Sobrevida e Gráficos Side-by-Side
# ==============================================================================
message("Rodando modelos de sobrevida e gerando gráficos combinados...")

soc_results <- list()

for(gene in target_genes) {
  
  # Pular se o gene não estiver presente em alguma das matrizes
  if(!gene %in% rownames(expr_axi) | !gene %in% rownames(expr_soc)) next
  
  # ---------------------------------------------------------
  # PARTE 1: ZUMA-7 (Axi-cel)
  # ---------------------------------------------------------
  axi_df <- clin_axi %>%
    mutate(
      GeneExpr = as.numeric(expr_axi[gene, clin_axi$sample]),
      Strata = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm=T), "High", "Low"), levels = c("Low", "High"))
    )
  
  fit_axi <- survfit(Surv(EFS_time, EFS_event) ~ Strata, data = axi_df)
  
  p_axi <- ggsurvplot(fit_axi, data = axi_df, pval = TRUE, risk.table = TRUE, 
                      title = paste("ZUMA-7 (Axi-cel)\nEFS by", gene),
                      xlab = "Event-Free Survival (Months)", ylab = "Probability",
                      palette = c("#4DBBD5FF", "#E64B35FF"), legend = "none",
                      ggtheme = theme_classic())
  
  # ---------------------------------------------------------
  # PARTE 2: ZUMA-7 (Standard of Care)
  # ---------------------------------------------------------
  soc_df <- clin_soc %>%
    mutate(
      GeneExpr = as.numeric(expr_soc[gene, clin_soc$sample]),
      Strata = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm=T), "High", "Low"), levels = c("Low", "High"))
    )
  
  fit_soc <- survfit(Surv(EFS_time, EFS_event) ~ Strata, data = soc_df)
  surv_diff_soc <- survdiff(Surv(EFS_time, EFS_event) ~ Strata, data = soc_df)
  p_val_soc <- 1 - pchisq(surv_diff_soc$chisq, length(surv_diff_soc$n) - 1)
  
  soc_results[[gene]] <- tibble(Gene = gene, SOC_KM_PValue = p_val_soc)
  
  p_soc <- ggsurvplot(fit_soc, data = soc_df, pval = TRUE, risk.table = TRUE, 
                      title = paste("ZUMA-7 (Standard of Care)\nEFS by", gene),
                      xlab = "Event-Free Survival (Months)", ylab = "Probability",
                      palette = c("#4DBBD5FF", "#E64B35FF"), legend.title = "Expression", legend.labs = c("Low", "High"),
                      ggtheme = theme_classic())
  
  # ---------------------------------------------------------
  # PARTE 3: Unir Gráficos (Patchwork)
  # ---------------------------------------------------------
  plot_combined <- (p_axi$plot + p_soc$plot) / (p_axi$table + p_soc$table) + 
    plot_layout(heights = c(4, 1))
  
  ggsave(file.path(out_dir, "KM_SideBySide", paste0("AxiCel_vs_SOC_", gene, ".pdf")), plot_combined, width = 10, height = 6)
  ggsave(file.path(out_dir, "KM_SideBySide", paste0("AxiCel_vs_SOC_", gene, ".png")), plot_combined, width = 10, height = 6, dpi = 300)
}

# Salvar tabela com os p-valores do SOC
soc_results_df <- bind_rows(soc_results) %>% arrange(SOC_KM_PValue)
write_csv(soc_results_df, file.path(out_dir, "ZUMA7_SOC_Validation_Results.csv"))

message("\n==========================================================")
message("Comparação Axi-cel vs SOC concluída com sucesso!")
message("Verifique a pasta: ", file.path(out_dir, "KM_SideBySide"))
message("==========================================================")
