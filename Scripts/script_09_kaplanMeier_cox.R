################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 09 - Validação Prognóstica (ZUMA-7) e Intersecção de DEGs
#
# Objetivo:
#   1. Extrair DEGs do ZUMA-1 (Base), ZUMA-7 (Base) e ZUMA-1 (Pareado).
#   2. Gerar Diagrama de Venn e Tabela de intersecção dos DEGs.
#   3. Testar impacto na sobrevida (EFS) no ZUMA-7 via Kaplan-Meier e Cox Univariado.
#   4. Filtrar genes duplamente significativos (KM p<0.05 E Uni_Cox p<0.05).
#   5. Rodar Cox Multivariado ajustado (Tumor Burden SPD e Cell-of-Origin).
#   6. Gerar Forest Plot e curvas de Kaplan-Meier individuais.
################################################################################

## =============================================================================
## 1. Pacotes e Diretórios
## =============================================================================
# Se não tiver o ggVennDiagram, instale com: install.packages("ggVennDiagram")
library(tidyverse)
library(survival)
library(survminer)
library(openxlsx)
library(ggplot2)
library(RColorBrewer)
library(ggVennDiagram)

out_dir <- "results/survival_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "KM_plots"), recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## 2. Carregar Tabelas de DEA e Definir DEGs
## =============================================================================
message("Carregando resultados do Limma (Script 07)...")

dea_z1 <- read_csv("results/dea_gsea/tables/DEA_ZUMA1_Baseline_NonCR_vs_CR_RawP.csv", show_col_types = FALSE)
dea_z7 <- read_csv("results/dea_gsea/tables/DEA_ZUMA7_AxiCel_Others_vs_Ongoing_RawP.csv", show_col_types = FALSE)
dea_z1_paired <- read_csv("results/dea_gsea/tables/DEA_ZUMA1_Paired_Post_vs_Baseline_RawP.csv", show_col_types = FALSE)

# Função para extrair os genes significativos (p < 0.05 e |logFC| > 0.5)
get_sig_genes <- function(df) {
  df %>% filter(P.Value < 0.05, abs(logFC) > 0.5, !is.na(Feature)) %>% pull(Feature) %>% unique()
}

deg_list <- list(
  "ZUMA-1 Baseline" = get_sig_genes(dea_z1),
  "ZUMA-7 Baseline" = get_sig_genes(dea_z7),
  "ZUMA-1 Paired"   = get_sig_genes(dea_z1_paired)
)

## =============================================================================
## 3. Diagrama de Venn e Tabela de Intersecção
## =============================================================================
message("Gerando Diagrama de Venn dos DEGs...")

# 3.1 Plot do Venn
p_venn <- ggVennDiagram(deg_list, label_alpha = 0, edge_size = 1) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4DBBD5FF") +
  labs(title = "Overlap of Differentially Expressed Genes", fill = "Gene Count") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "Venn_Diagram_DEGs.pdf"), p_venn, width = 7, height = 6)
ggsave(file.path(out_dir, "Venn_Diagram_DEGs.png"), p_venn, width = 7, height = 6, dpi = 300)

# 3.2 Tabela de Intersecção
all_degs <- unique(unlist(deg_list))

deg_table <- tibble(Gene = all_degs) %>%
  mutate(
    In_ZUMA1_Base   = Gene %in% deg_list$`ZUMA-1 Baseline`,
    In_ZUMA7_Base   = Gene %in% deg_list$`ZUMA-7 Baseline`,
    In_ZUMA1_Paired = Gene %in% deg_list$`ZUMA-1 Paired`
  ) %>%
  mutate(
    Source = paste(
      ifelse(In_ZUMA1_Base, "ZUMA1_Base", ""),
      ifelse(In_ZUMA7_Base, "ZUMA7_Base", ""),
      ifelse(In_ZUMA1_Paired, "ZUMA1_Paired", ""),
      sep = " | "
    ) %>% str_replace_all("^ \\| | \\| $", "") %>% str_replace_all(" \\|  \\| ", " | ")
  ) %>% arrange(desc(In_ZUMA1_Base + In_ZUMA7_Base + In_ZUMA1_Paired), Gene)

# ALTERAÇÃO: Salvando em CSV
write_csv(deg_table, file.path(out_dir, "DEG_Intersection_Table.csv"))

## =============================================================================
## 4. Preparar Dados Clínicos do ZUMA-7 para Sobrevida (CORRIGIDO)
## =============================================================================
message("Preparando dados do ZUMA-7 para análise de sobrevida...")

z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
expr_z7 <- z7_mrna$expr

# Mapeamento direto com os nomes exatos do seu dataframe
clin_z7 <- z7_mrna$clin %>%
  rownames_to_column("sample") %>%
  mutate(
    # Extraindo as colunas com base no print
    EFS_time = as.numeric(event_free_survival_months_ch1),
    EFS_event = as.numeric(event_free_survival_event_ch1),
    TumorBurden = as.numeric(baseline_tumor_burden_spd_ch1),
    COO = factor(cell_of_origin_ch1)
  ) %>%
  # Limpar quem não tem dados de sobrevida registrados
  filter(!is.na(EFS_time), !is.na(EFS_event))

valid_samples <- intersect(clin_z7$sample, colnames(expr_z7))
clin_z7 <- clin_z7 %>% filter(sample %in% valid_samples)
expr_z7 <- expr_z7[, clin_z7$sample, drop = FALSE]

candidate_genes <- intersect(all_degs, rownames(expr_z7))
message(length(candidate_genes), " genes candidatos encontrados na matriz de expressão do ZUMA-7.")

## =============================================================================
## 5. Screening Inicial: Kaplan-Meier e Cox Univariado
## =============================================================================
message("Rodando Screening Univariado (KM e Cox) para todos os genes...")

run_screening <- function(gene) {
  gene_expr <- as.numeric(expr_z7[gene, clin_z7$sample])
  
  df <- clin_z7 %>%
    mutate(
      GeneExpr = gene_expr,
      GeneGroup = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm = TRUE), "High", "Low"), levels = c("Low", "High"))
    ) %>%
    filter(!is.na(GeneExpr))
  
  if (nrow(df) < 20 || n_distinct(df$GeneGroup) < 2) return(NULL)
  
  # 1. Kaplan Meier
  surv_diff <- survdiff(Surv(EFS_time, EFS_event) ~ GeneGroup, data = df)
  km_p <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  
  # 2. Cox Univariado (usando a expressão contínua)
  fit_uni <- tryCatch(coxph(Surv(EFS_time, EFS_event) ~ GeneExpr, data = df), error = function(e) NULL)
  
  uni_p <- NA; uni_hr <- NA
  if (!is.null(fit_uni)) {
    s <- summary(fit_uni)
    if("GeneExpr" %in% rownames(s$coefficients)) {
      uni_p <- s$coefficients["GeneExpr", "Pr(>|z|)"]
      uni_hr <- s$coefficients["GeneExpr", "exp(coef)"]
    }
  }
  
  tibble(Gene = gene, KM_PValue = km_p, UniCox_PValue = uni_p, UniCox_HR = uni_hr)
}

screening_results <- map_dfr(candidate_genes, run_screening)

# Filtrar o Funil: Apenas genes significativos NOS DOIS testes (< 0.05)
priority_genes <- screening_results %>%
  filter(KM_PValue < 0.05, UniCox_PValue < 0.05) %>%
  arrange(UniCox_PValue) %>%
  pull(Gene)

message("O Funil filtrou: ", length(priority_genes), " genes duplamente significativos (KM e Univariado).")

## =============================================================================
## 6. Validação Final: Cox Multivariado (Ajustado por SPD e COO)
## =============================================================================
if(length(priority_genes) > 0) {
  message("Rodando modelo Multivariado para os genes prioritários...")
  
  run_multivariate <- function(gene) {
    df <- clin_z7 %>%
      mutate(GeneExpr = as.numeric(expr_z7[gene, clin_z7$sample])) %>%
      filter(!is.na(GeneExpr), !is.na(TumorBurden), !is.na(COO), COO != "Unknown")
    
    fit_multi <- tryCatch(coxph(Surv(EFS_time, EFS_event) ~ GeneExpr + TumorBurden + COO, data = df), error = function(e) NULL)
    
    if (is.null(fit_multi)) return(NULL)
    
    s <- summary(fit_multi)
    if (!"GeneExpr" %in% rownames(s$coefficients)) return(NULL)
    
    tibble(
      Gene = gene,
      MultiCox_HR = s$coefficients["GeneExpr", "exp(coef)"],
      CI_lower = s$conf.int["GeneExpr", "lower .95"],
      CI_upper = s$conf.int["GeneExpr", "upper .95"],
      MultiCox_PValue = s$coefficients["GeneExpr", "Pr(>|z|)"]
    )
  }
  
  multi_results <- map_dfr(priority_genes, run_multivariate) %>%
    mutate(
      Risk_category = factor(ifelse(MultiCox_HR > 1, "Risk (Worse EFS)", "Protective (Better EFS)"), 
                             levels = c("Protective (Better EFS)", "Risk (Worse EFS)")),
      label = sprintf("HR %.2f; p = %.3f", MultiCox_HR, MultiCox_PValue)
    ) %>%
    arrange(MultiCox_PValue)
  
  # Salvar Tabela Mestra
  final_table <- screening_results %>% 
    left_join(multi_results, by = "Gene") %>% 
    left_join(deg_table, by = "Gene")
  
  # ALTERAÇÃO: Salvando em CSV
  write_csv(final_table, file.path(out_dir, "ZUMA7_Survival_Validation_Master.csv"))
  
  ## ===========================================================================
  ## 7. Forest Plot e Curvas KM
  ## ===========================================================================
  
  # Forest Plot
  plot_df <- multi_results %>% 
    filter(MultiCox_PValue < 0.05) %>% # Plota apenas os que sobreviveram ao ajuste multivariado
    mutate(Gene = factor(Gene, levels = rev(Gene)))
  
  if(nrow(plot_df) > 0) {
    p_forest <- ggplot(plot_df, aes(x = MultiCox_HR, y = Gene, color = Risk_category)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, linewidth = 0.8) +
      geom_point(size = 3) +
      geom_text(aes(label = label), hjust = -0.15, size = 3.5, show.legend = FALSE, color="black") +
      scale_x_log10() +
      scale_color_manual(values = c("Protective (Better EFS)" = "#2166AC", "Risk (Worse EFS)" = "#E64B35FF")) +
      theme_bw(base_size = 12) +
      labs(title = "Prognostic Validation of Candidate DEGs in ZUMA-7",
           subtitle = "Multivariate Cox proportional hazards adjusted by Tumor Burden (SPD) and Cell-of-Origin",
           x = "Hazard Ratio (Log10 Scale)", y = NULL, color = "Impact on EFS") +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
      coord_cartesian(clip = "off")
    
    ggsave(file.path(out_dir, "ForestPlot_Multivariate_ZUMA7.pdf"), p_forest, width = 9, height = max(4, nrow(plot_df)*0.4))
  }
  
  # KM Plots Individuais
  for(g in plot_df$Gene) {
    df_km <- clin_z7 %>%
      mutate(GeneExpr = as.numeric(expr_z7[g, clin_z7$sample])) %>%
      mutate(Strata = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm=T), "High", "Low"), levels = c("Low", "High"))) %>%
      filter(!is.na(GeneExpr))
    
    fit <- survfit(Surv(EFS_time, EFS_event) ~ Strata, data = df_km)
    
    p_km <- ggsurvplot(fit, data = df_km, pval = TRUE, risk.table = TRUE, title = paste("Survival by", g, "Expression"),
                       xlab = "Event-Free Survival (Months)", ylab = "Probability",
                       palette = c("#4DBBD5FF", "#E64B35FF"), legend.title = g, legend.labs = c("Low", "High"),
                       ggtheme = theme_classic())
    
    pdf(file.path(out_dir, "KM_plots", paste0("KM_", g, ".pdf")), width = 5.5, height = 6)
    print(p_km)
    dev.off()
  }
  
} else {
  message("Nenhum gene sobreviveu ao filtro Univariado (KM + Cox). O modelo Multivariado não será rodado.")
}

message("\n==========================================================")
message("Script 09 concluído! Verifique a pasta: ", out_dir)
message("==========================================================")