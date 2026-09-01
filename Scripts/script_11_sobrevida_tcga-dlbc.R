################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 11 - Validação de Biomarcadores na Coorte TCGA-DLBC (Quimioterapia)
#
# Objetivo:
#   1. Criar diretórios e validar a presença dos ficheiros do TCGA.
#   2. Carregar TPM e Sobrevida do TCGA-DLBC (Xenabrowser).
#   3. Converter Ensembl IDs para Gene Symbols.
#   4. Testar a sobrevida global (OS) para os genes validados no ZUMA-7.
#   5. Gerar gráficos Side-by-Side (ZUMA-7 EFS vs. TCGA OS).
################################################################################

# ==============================================================================
# 1. Pacotes e Diretórios de Output
# ==============================================================================
library(tidyverse)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)

out_dir <- "results/tcga_validation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "KM_SideBySide"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1.5 Preparação Automática: Download de Ficheiros do Xenabrowser (TCGA-DLBC)
# ==============================================================================
tcga_dir <- "data/TCGA"
dir.create(tcga_dir, recursive = TRUE, showWarnings = FALSE)

# URLs diretas do GDC Hub (Xenabrowser) para o TCGA-DLBC
url_surv <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-DLBC.survival.tsv.gz"
url_expr <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-DLBC.star_tpm.tsv.gz"

# Caminhos locais para salvar os ficheiros
file_surv <- file.path(tcga_dir, "TCGA-DLBC.survival.tsv")
file_expr <- file.path(tcga_dir, "TCGA-DLBC.star_tpm.tsv.gz")

# Download automático da Sobrevida (se não existir)
if (!file.exists(file_surv)) {
  message("A descarregar dados clínicos de sobrevida do TCGA-DLBC...")
  download.file(url = url_surv, destfile = file_surv, mode = "wb", quiet = FALSE)
} else {
  message("Ficheiro de sobrevida já existe localmente.")
}

# Download automático da Matriz de Expressão TPM (se não existir)
if (!file.exists(file_expr)) {
  message("A descarregar matriz de expressão TPM do TCGA-DLBC (isto pode demorar alguns minutos)...")
  # Aumentar o timeout temporariamente para ficheiros grandes
  options(timeout = max(300, getOption("timeout"))) 
  download.file(url = url_expr, destfile = file_expr, mode = "wb", quiet = FALSE)
} else {
  message("Matriz de expressão TPM já existe localmente.")
}

# Atualização na Secção 3: Ler os ficheiros diretamente (o read_tsv lida com .gz nativamente)
message("A carregar dados do TCGA-DLBC para a memória...")
tcga_clin <- read_tsv(file_surv, show_col_types = FALSE)
tcga_expr <- read_tsv(file_expr, show_col_types = FALSE)

# ==============================================================================
# 2. Definir a Lista de Genes de Interesse (Ensembl IDs)
# ==============================================================================
message("Carregando lista de genes prioritários definidos manualmente (Ensembl IDs)...")

ensg_interesse <- c("ENSG00000138166","ENSG00000241106", "ENSG00000156738", 
                    "ENSG00000134539","ENSG00000170558","ENSG00000140968",
                    "ENSG00000122224","ENSG00000085741","ENSG00000081237",
                    "ENSG00000177455","ENSG00000048462","ENSG00000149294",
                    "ENSG00000115009","ENSG00000186265","ENSG00000165949",
                    "ENSG00000137507")

# ==============================================================================
# 3. Carregar Dados Clínicos e de Expressão (TCGA-DLBC)
# ==============================================================================
message("Carregando dados do TCGA-DLBC...")

tcga_clin <- read_tsv(file_surv, show_col_types = FALSE)
tcga_expr <- read_tsv(file_expr, show_col_types = FALSE)

# ==============================================================================
# 4. Processamento dos Dados TCGA
# ==============================================================================
message("Filtrando matriz do TCGA e anotando para Gene Symbols...")

# 4.1 Limpar os Ensembl IDs da matriz do TCGA (remover a versão após o ponto, ex: .15)
expr_mat <- tcga_expr %>%
  mutate(Ensembl_Clean = str_remove(Ensembl_ID, "\\..*$"))

# 4.2 Filtrar a matriz gigante APENAS para os seus 16 genes de interesse
expr_filtered <- expr_mat %>%
  filter(Ensembl_Clean %in% ensg_interesse)

# 4.3 Traduzir os Ensembl IDs encontrados para Gene Symbol (para o gráfico ficar legível)
gene_map <- bitr(expr_filtered$Ensembl_Clean, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

# 4.4 Reduzir a matriz, trocar a coluna Ensembl pelo Gene Symbol e organizar
expr_target <- expr_filtered %>%
  inner_join(gene_map, by = c("Ensembl_Clean" = "ENSEMBL")) %>%
  select(-Ensembl_ID, -Ensembl_Clean) %>%
  # Agrupar e tirar a média (caso algum ID mapeie para dois símbolos por redundância da anotação)
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("SYMBOL")

# 4.5 Preparar dados clínicos do TCGA (Converter dias para meses para alinhar com o ZUMA-7)
tcga_clin_clean <- tcga_clin %>%
  mutate(
    OS_months = OS.time / 30.416,  # Converte dias para meses
    OS_event = OS
  ) %>%
  filter(!is.na(OS_months), !is.na(OS_event))

# Garantir que as amostras da clínica batem perfeitamente com a expressão
valid_samples_tcga <- intersect(tcga_clin_clean$sample, colnames(expr_target))
tcga_clin_clean <- tcga_clin_clean %>% filter(sample %in% valid_samples_tcga)
expr_target <- expr_target[, tcga_clin_clean$sample, drop = FALSE]

# ==============================================================================
# 5. Carregar Dados ZUMA-7 para Comparação Lado a Lado
# ==============================================================================
message("Carregando dados do ZUMA-7 para os painéis comparativos...")

z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
expr_z7 <- z7_mrna$expr
clin_z7 <- z7_mrna$clin %>%
  rownames_to_column("sample") %>%
  mutate(
    EFS_time = as.numeric(event_free_survival_months_ch1),
    EFS_event = as.numeric(event_free_survival_event_ch1)
  ) %>%
  filter(!is.na(EFS_time), !is.na(EFS_event))

valid_samples_z7 <- intersect(clin_z7$sample, colnames(expr_z7))
clin_z7 <- clin_z7 %>% filter(sample %in% valid_samples_z7)
expr_z7 <- expr_z7[, clin_z7$sample, drop = FALSE]

# ==============================================================================
# 6. Análise de Sobrevida e Gráficos Side-by-Side
# ==============================================================================
message("Correndo sobrevida para TCGA e gerando gráficos combinados...")

tcga_results <- list()

for(gene in rownames(expr_target)) {
  
  # ---------------------------------------------------------
  # PARTE 1: TCGA (Overall Survival)
  # ---------------------------------------------------------
  tcga_df <- tcga_clin_clean %>%
    mutate(
      GeneExpr = as.numeric(expr_target[gene, tcga_clin_clean$sample]),
      Strata = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm=T), "High", "Low"), levels = c("Low", "High"))
    )
  
  fit_tcga <- survfit(Surv(OS_months, OS_event) ~ Strata, data = tcga_df)
  surv_diff_tcga <- survdiff(Surv(OS_months, OS_event) ~ Strata, data = tcga_df)
  p_val_tcga <- 1 - pchisq(surv_diff_tcga$chisq, length(surv_diff_tcga$n) - 1)
  
  # Guardar o resultado estatístico do TCGA
  tcga_results[[gene]] <- tibble(Gene = gene, TCGA_KM_PValue = p_val_tcga)
  
  p_tcga <- ggsurvplot(fit_tcga, data = tcga_df, pval = TRUE, risk.table = TRUE, 
                       title = paste("TCGA-DLBC (R-CHOP)\nOS by", gene),
                       xlab = "Overall Survival (Months)", ylab = "Probability",
                       palette = c("#4DBBD5FF", "#E64B35FF"), legend = "none",
                       ggtheme = theme_classic())
  
  # ---------------------------------------------------------
  # PARTE 2: ZUMA-7 (Event-Free Survival)
  # ---------------------------------------------------------
  z7_df <- clin_z7 %>%
    mutate(
      GeneExpr = as.numeric(expr_z7[gene, clin_z7$sample]),
      Strata = factor(ifelse(GeneExpr >= median(GeneExpr, na.rm=T), "High", "Low"), levels = c("Low", "High"))
    )
  
  fit_z7 <- survfit(Surv(EFS_time, EFS_event) ~ Strata, data = z7_df)
  
  p_z7 <- ggsurvplot(fit_z7, data = z7_df, pval = TRUE, risk.table = TRUE, 
                     title = paste("ZUMA-7 (Axi-cel)\nEFS by", gene),
                     xlab = "Event-Free Survival (Months)", ylab = "Probability",
                     palette = c("#4DBBD5FF", "#E64B35FF"), legend.title = "Expression", legend.labs = c("Low", "High"),
                     ggtheme = theme_classic())
  
  # ---------------------------------------------------------
  # PARTE 3: Unir Gráficos (Patchwork)
  # ---------------------------------------------------------
  # Extrair os objetos ggplot reais de dentro da lista gerada pelo ggsurvplot
  plot_combined <- (p_z7$plot + p_tcga$plot) / (p_z7$table + p_tcga$table) + 
    plot_layout(heights = c(4, 1))
  
  ggsave(file.path(out_dir, "KM_SideBySide", paste0("Comparison_", gene, ".pdf")), plot_combined, width = 10, height = 6)
  ggsave(file.path(out_dir, "KM_SideBySide", paste0("Comparison_", gene, ".png")), plot_combined, width = 10, height = 6, dpi = 300)
}

# Guardar tabela com os p-valores do TCGA
tcga_results_df <- bind_rows(tcga_results) %>% arrange(TCGA_KM_PValue)
write_csv(tcga_results_df, file.path(out_dir, "TCGA_Validation_Results.csv"))

message("\n==========================================================")
message("Validação no TCGA concluída com sucesso!")
message("Verifique a pasta: ", file.path(out_dir, "KM_SideBySide"))
message("==========================================================")