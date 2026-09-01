################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19
# Script: 05 - Análise Multivariada (Gene vs Composição Celular)
# Objetivo:
#   1. Avaliar a independência da Assinatura de Resistência face ao Microambiente
#   2. Rodar modelos de Cox Multivariados: Surv(EFS) ~ Gene + Celula
#   3. Extrair Hazard Ratios (HR) e P-values ajustados
################################################################################

## =============================================================================
## 1. Pacotes
## =============================================================================
library(tidyverse)
library(survival)
library(survminer)
library(broom)

options(stringsAsFactors = FALSE)
dir.create("results/multivariate_cox", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados do ZUMA-7
## =============================================================================
message("A carregar dados do ZUMA-7...")

# Matriz purificada (do Script 02)
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
expr_matrix <- z7_mrna$expr
clin_df <- z7_mrna$clin %>% rownames_to_column("sample")

# Escores celulares (do Script 03) - Usando MCP-counter como exemplo principal
mcp_z7 <- read_csv("results/cellular_composition/ZUMA7_AxiCel_mRNA_MCPcounter_scores.csv")

## =============================================================================
## 3. Preparar o Dataframe Integrado
## =============================================================================
## =============================================================================
## 3. Preparar o Dataframe Integrado
## =============================================================================
time_col <- "event_free_survival_months_ch1"   
event_col <- "event_free_survival_event_ch1"

signature_genes <- c("CD19", "MS4A1", "TNFRSF17", "DUSP5", "BTLA", "WNT11", 
                     "CDH2", "HLA-DOB", "IRF8", "KLRD1", "CD45RA", "IFI27", 
                     "LRRC32", "LY9", "NCAM1", "CCL20")

# Extrair genes disponíveis e transpor
valid_genes <- intersect(signature_genes, rownames(expr_matrix))
expr_subset <- t(expr_matrix[valid_genes, , drop = FALSE]) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

# Populações celulares que SABEMOS que existem no seu mcp_z7
valid_mcp_cols <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", 
                    "B lineage", "NK cells", "Monocytic lineage", 
                    "Neutrophils", "Endothelial cells")

# Juntar Clínica + Expressão + Populações MCP-counter
df_model <- clin_df %>%
  select(sample, time = all_of(time_col), event = all_of(event_col)) %>%
  inner_join(expr_subset, by = "sample") %>%
  inner_join(mcp_z7 %>% select(sample, all_of(valid_mcp_cols)), by = "sample")

# Converter tempo e evento para numérico
df_model$time <- as.numeric(df_model$time)
df_model$event <- as.numeric(df_model$event)

# Limpar nomes das colunas (trocar espaços/traços por sublinhados para o Cox não dar erro)
colnames(df_model) <- gsub(" ", "_", colnames(df_model))
colnames(df_model) <- gsub("-", "_", colnames(df_model))

# Limpar também os vetores de nomes para o loop
valid_genes_clean <- gsub("-", "_", valid_genes)
covariates_clean <- gsub(" ", "_", valid_mcp_cols)

## =============================================================================
## 4. Rodar Modelos de Cox Multivariados
## =============================================================================
message("A calcular regressões multivariadas...")

results_list <- list()

for (gene in valid_genes_clean) {
  for (covar in covariates_clean) {
    
    # Criar a fórmula: Surv(time, event) ~ Gene + Celula
    fmla <- as.formula(paste("Surv(time, event) ~", gene, "+", covar))
    
    # Ajustar o modelo de Cox
    fit <- tryCatch(coxph(fmla, data = df_model), error = function(e) NULL)
    
    if (!is.null(fit)) {
      # Extrair os resultados formatados com o pacote broom
      res <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == gene) %>% # Queremos saber apenas se o GENE sobreviveu ao ajuste
        mutate(
          Target_Gene = gene,
          Adjusted_By = covar,
          Independent = ifelse(p.value < 0.05, "YES", "NO")
        ) %>%
        select(Target_Gene, Adjusted_By, HR = estimate, Conf.Low = conf.low, Conf.High = conf.high, P_Value = p.value, Independent)
      
      results_list[[paste(gene, covar, sep = "_")]] <- res
    }
  }
}

# Consolidar numa única tabela
cox_multivariate_results <- bind_rows(results_list) %>% arrange(Target_Gene, P_Value)

# Salvar o xeque-mate
write_csv(cox_multivariate_results, "results/multivariate_cox/ZUMA7_Multivariate_Cox_Independence.csv")

cat("\n==================================================================\n")
cat("Sucesso! Tabela de Cox Multivariada guardada na pasta results/multivariate_cox/\n")
cat("==================================================================\n")