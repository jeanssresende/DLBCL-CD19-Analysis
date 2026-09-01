################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 12b - Extração do Braço SOC a partir do RDS Processado do GEO
# Objetivo: Gerar o arquivo ZUMA7_SOC_mRNA_Clean.rds
################################################################################

library(tidyverse)

# ==============================================================================
# 1. Carregar o Ficheiro Completo (Axi-cel + SOC)
# ==============================================================================
message("A carregar o ficheiro processado do ZUMA-7...")

z7_full <- readRDS("data/ZUMA7_GSE248835_mRNA_processed.rds")

# CORREÇÃO: Usando os nomes exatos que estão dentro da lista
expr_full <- z7_full$expression_log2
clin_full <- z7_full$clinical

# ==============================================================================
# 2. Extrair o Braço Standard of Care (SOC)
# ==============================================================================
message("A isolar os pacientes do braço Standard of Care...")

clin_soc <- clin_full %>%
  rownames_to_column("sample") %>%
  # Usar o str_detect para evitar problemas com espaços invisíveis
  filter(str_detect(treatment_arm_ch1, "Standard of Care Chemotherapy")) %>% 
  column_to_rownames("sample")

# ==============================================================================
# 3. Alinhar a Matriz de Expressão e Salvar
# ==============================================================================
message("A cruzar os dados clínicos com a matriz de expressão...")

# Garantir que pegamos apenas as colunas (amostras) dos pacientes SOC
valid_samples_soc <- intersect(rownames(clin_soc), colnames(expr_full))

clin_soc_final <- clin_soc[valid_samples_soc, ]
expr_soc_final <- expr_full[, valid_samples_soc, drop = FALSE]

# Mantemos os nomes 'expr' e 'clin' aqui para que o Script 12 (que já escrevemos) funcione perfeitamente
z7_soc_clean <- list(
  expr = expr_soc_final,
  clin = clin_soc_final
)

# Salvar dentro da pasta filtered
saveRDS(z7_soc_clean, "data/filtered/ZUMA7_SOC_mRNA_Clean.rds")

message("==========================================================")
message("Sucesso! Extraímos ", length(valid_samples_soc), " pacientes do braço SOC.")
message("O ficheiro 'ZUMA7_SOC_mRNA_Clean.rds' foi criado em data/filtered/")
message("Já pode correr o Script 12 completo!")
message("==========================================================")
