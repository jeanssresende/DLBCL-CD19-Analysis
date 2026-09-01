################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 01 – Download, Processamento e Harmonização de Dados
# Datasets: GSE197977 (ZUMA-1) e GSE248835 (ZUMA-7)
################################################################################

# =============================================================================
# 1. Configurações e Pacotes
# =============================================================================
library(GEOquery)
library(Biobase)
library(tidyverse)
library(limma)
library(janitor)

# Configurações globais para performance
options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 4)

# Criar diretório de dados se não existir
if (!dir.exists("data")) dir.create("data")

# =============================================================================
# 2. Função Core de Processamento
# =============================================================================
process_zuma_dataset <- function(gse_id, gene_col_name) {
  
  message(paste("\n--- Iniciando processamento do dataset:", gse_id, "---"))
  
  # 2.1 Download do dataset do GEO
  gds <- getGEO(gse_id, GSEMatrix = TRUE, destdir = "data/")
  e <- gds[[1]]
  rm(gds)
  
  # 2.2 Extração e Limpeza das Matrizes
  gexp   <- exprs(e)
  fdat   <- fData(e)
  
  # Limpeza dos nomes das colunas clínicas para o padrão snake_case
  pdat   <- pData(e) %>% clean_names() 
  
  stopifnot(identical(rownames(gexp), rownames(fdat)))
  
  # 2.3 Definição de nomes de genes (Rownames)
  if(gene_col_name %in% colnames(fdat)) {
    # Remove NAs e strings vazias
    valid_genes <- !is.na(fdat[[gene_col_name]]) & fdat[[gene_col_name]] != ""
    gexp <- gexp[valid_genes, ]
    fdat <- fdat[valid_genes, ]
    
    # Atribuir nomes usando a coluna selecionada
    rownames(gexp) <- fdat[[gene_col_name]]
    
    # Agregação de sondas duplicadas pela média
    if(any(duplicated(rownames(gexp)))) {
      message("Aviso: Sondas múltiplas mapeando para o mesmo gene. Agregando pela média...")
      gexp <- avereps(gexp)
    }
  } else {
    stop(paste("ERRO: Coluna", gene_col_name, "não encontrada na anotação de", gse_id))
  }
  
  # 2.4 Transformação Log2 Inteligente (Evita dupla transformação)
  max_val <- max(gexp, na.rm = TRUE)
  if (max_val > 50) {
    message("Status: Valores altos detectados (Max = ", round(max_val, 2), "). Aplicando log2(x + 1)...")
    gexp_log2 <- log2(gexp + 1)
  } else {
    message("Status: Matriz parece já estar em escala logarítmica (Max = ", round(max_val, 2), "). Mantendo os valores originais...")
    gexp_log2 <- gexp
  }
  
  # 2.5 Organização do Objeto Final
  output_list <- list(
    expression_raw  = gexp,            # Matriz original como baixada
    expression_log2 = gexp_log2,       # Matriz final pronta para análise
    clinical        = pdat,
    feature_data    = fdat,
    metadata = list(
      gse = gse_id,
      platform = annotation(e),
      date_processed = Sys.Date()
    )
  )
  
  message(paste("Sucesso:", gse_id, "processado com", nrow(gexp_log2), "genes e", ncol(gexp_log2), "amostras."))
  return(output_list)
}

# =============================================================================
# 3. Execução e Tratamento Específico
# =============================================================================

# ZUMA-1: Dataset apenas com genes individuais
zuma1 <- process_zuma_dataset("GSE197977", "ORF")

# ZUMA-7: Dataset Misto (Genes + Assinaturas)
zuma7_full <- process_zuma_dataset("GSE248835", "Gene_Signature_Name")

message("\n--- Iniciando separação do dataset misto ZUMA-7 ---")

# 3.1 Identificar os índices das linhas de mRNA e Assinaturas
idx_mrna <- zuma7_full$feature_data$`Analyte Type` == "mRNA"
idx_sigs <- zuma7_full$feature_data$`Analyte Type` == "IO360 Signature"

# 3.2 Criar objeto EXCLUSIVO para mRNA (Copiando a estrutura e filtrando)
zuma7_mrna <- zuma7_full
zuma7_mrna$expression_raw  <- zuma7_full$expression_raw[idx_mrna, ]
zuma7_mrna$expression_log2 <- zuma7_full$expression_log2[idx_mrna, ]
zuma7_mrna$feature_data    <- zuma7_full$feature_data[idx_mrna, ]

# 3.3 Criar objeto EXCLUSIVO para Assinaturas
zuma7_sigs <- zuma7_full
zuma7_sigs$expression_raw  <- zuma7_full$expression_raw[idx_sigs, ]
zuma7_sigs$expression_log2 <- zuma7_full$expression_log2[idx_sigs, ]
zuma7_sigs$feature_data    <- zuma7_full$feature_data[idx_sigs, ]

message("ZUMA-7 mRNA: ", nrow(zuma7_mrna$expression_log2), " genes isolados.")
message("ZUMA-7 Sigs: ", nrow(zuma7_sigs$expression_log2), " assinaturas isoladas.")

# =============================================================================
# 4. Exportação dos Resultados
# =============================================================================

saveRDS(zuma1, file = "./data/ZUMA1_GSE197977_processed.rds")
saveRDS(zuma7_mrna, file = "./data/ZUMA7_GSE248835_mRNA_processed.rds")
saveRDS(zuma7_sigs, file = "./data/ZUMA7_GSE248835_Signatures_processed.rds")

cat("\n==========================================================\n")
cat("Processamento concluído. Objetos perfeitamente separados e prontos.\n")
cat("Arquivos salvos na pasta /data/\n")
cat("==========================================================\n")

