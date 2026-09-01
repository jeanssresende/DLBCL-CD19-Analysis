################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19
# Script: 04 - Correlação: Assinatura de Resistência vs Composição Celular
# Objetivo:
#   1. Carregar matrizes de expressão (Painéis NanoString) e escores celulares
#   2. Calcular correlação de Spearman (Assinatura vs Infiltração)
#   3. Calcular P-values e FDR para justificar independência biológica
#   4. Gerar e salvar Heatmaps de Correlação e Tabelas Estatísticas
################################################################################

## =============================================================================
## 1. Pacotes e Preparação
## =============================================================================
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

options(stringsAsFactors = FALSE)

# Criar diretório para as figuras e tabelas de correlação
dir.create("results/cellular_correlation", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Objetos de Expressão (Painéis NanoString)
## =============================================================================
message("Carregando matrizes de expressão purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")

## =============================================================================
## 3. Carregar Escores de Composição Celular (Gerados no Script 03)
## =============================================================================
message("Carregando escores de composição celular...")

# Ajuste os nomes dos arquivos CSV consoante o que foi salvo no Script 03
mcp_z1b <- read_csv("results/cellular_composition/ZUMA1_Baseline_MCPcounter_scores.csv")
mcp_z7  <- read_csv("results/cellular_composition/ZUMA7_AxiCel_mRNA_MCPcounter_scores.csv")

man_z1b <- read_csv("results/cellular_composition/ZUMA1_Baseline_Manual_Scores.csv")
man_z7  <- read_csv("results/cellular_composition/ZUMA7_AxiCel_mRNA_Manual_Scores.csv")

## =============================================================================
## 4. Definir Genes e Populações Celulares Alvo
## =============================================================================
# A sua assinatura completa
signature_genes <- c("CD19", "MS4A1", "TNFRSF17", "DUSP5", "BTLA", "WNT11", 
                     "CDH2", "HLA-DOB", "IRF8", "KLRD1", "CD45RA", "IFI27", 
                     "LRRC32", "LY9", "NCAM1", "CCL20")

# Colunas numéricas que queremos testar contra a assinatura
mcp_cols <- c("T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", 
              "B lineage", "Monocytic lineage", "Myeloid dendritic cells", 
              "Endothelial cells", "Fibroblasts")

man_cols <- c("B_lineage_proxy", "T_cell", "Cytotoxic", "T_cell_exhaustion", 
              "Myeloid", "HLA_class_I", "HLA_class_II", "IFN_gamma", 
              "immune_proxy", "tumor_immune_ratio_proxy")

## =============================================================================
## 5. Função Principal: Correlação Estatística e Heatmap
## =============================================================================
run_correlation_analysis <- function(expr_list, comp_df, cell_cols, cohort_name, analysis_type) {
  
  message(paste("\nAnalisando", cohort_name, "-", analysis_type))
  
  # 5.1 Isolar os genes da assinatura que REALMENTE existem neste painel NanoString
  expr_matrix <- expr_list$expr
  valid_genes <- intersect(signature_genes, rownames(expr_matrix))
  
  message(paste("Genes encontrados no painel:", length(valid_genes), "de", length(signature_genes)))
  if(length(valid_genes) == 0) stop("Nenhum gene da assinatura encontrado no painel!")
  
  # Preparar dataframe de expressão transposta
  expr_subset <- t(expr_matrix[valid_genes, , drop = FALSE]) %>% 
    as.data.frame() %>% 
    rownames_to_column("sample")
  
  # 5.2 Cruzar Expressão com os Escores Celulares
  # Interceptar apenas as células que existem na tabela
  valid_cells <- intersect(cell_cols, colnames(comp_df))
  
  merged_df <- expr_subset %>%
    inner_join(comp_df %>% select(sample, all_of(valid_cells)), by = "sample")
  
  # 5.3 Gerar Tabela Estatística Tidy (Gene x Célula, com P-value e FDR)
  # Isso é fundamental para responder rigorosamente aos revisores
  cor_stats <- expand_grid(Gene = valid_genes, Cell_Type = valid_cells) %>%
    mutate(
      res = map2(Gene, Cell_Type, function(g, c) {
        test <- cor.test(merged_df[[g]], merged_df[[c]], method = "spearman", exact = FALSE)
        tibble(rho = test$estimate, p_value = test$p.value)
      })
    ) %>%
    unnest(res) %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  write_csv(cor_stats, paste0("results/cellular_correlation/Stats_", cohort_name, "_", analysis_type, ".csv"))
  
  # ====================================================================
  # 5.4 Preparar Matrizes para o Heatmap (Usando o FDR exato da tabela)
  # ====================================================================
  cor_matrix <- matrix(NA, nrow = length(valid_genes), ncol = length(valid_cells),
                       dimnames = list(valid_genes, valid_cells))
  
  display_matrix <- matrix("", nrow = length(valid_genes), ncol = length(valid_cells),
                           dimnames = list(valid_genes, valid_cells))
  
  # Vamos preencher as matrizes usando a tabela cor_stats onde o FDR já está calculado
  for(i in 1:nrow(cor_stats)) {
    g <- cor_stats$Gene[i]
    c <- cor_stats$Cell_Type[i]
    
    rho_val <- cor_stats$rho[i]
    fdr_val <- cor_stats$FDR[i] # AQUI ESTÁ A CORREÇÃO: Puxando o FDR!
    
    cor_matrix[g, c] <- rho_val
    
    # Nova regra de ouro da significância BASEADA NO FDR
    stars <- ifelse(fdr_val < 0.001, "***",
                    ifelse(fdr_val < 0.01, "**",
                           ifelse(fdr_val < 0.05, "*", "")))
    
    # Junta o número arredondado com as estrelas (ex: "0.52***" ou "-0.12")
    display_matrix[g, c] <- paste0(round(rho_val, 2), stars)
  }
  
  # 5.5 Plotar Heatmap de Alta Qualidade
  # Centralizar a paleta de cores no 0 (Azul = Negativo, Branco = Neutro, Vermelho = Positivo)
  breaks_centered <- seq(-1, 1, length.out = 100)
  color_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
  
  pheatmap(cor_matrix,
           color = color_palette,
           breaks = breaks_centered,
           display_numbers = display_matrix, # AQUI PASSAMOS A NOSSA NOVA MATRIZ DE TEXTO
           fontsize_number = 8,
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           main = paste("Spearman Correlation:", cohort_name, "(", analysis_type, ")"),
           filename = paste0("results/cellular_correlation/Heatmap_", cohort_name, "_", analysis_type, ".pdf"),
           width = 10, height = 8)
  
  message(">> Tabelas e Heatmap gerados com sucesso!")
}
## =============================================================================
## 6. Execução
## =============================================================================

# ZUMA-1
run_correlation_analysis(z1_base, mcp_z1b, mcp_cols, "ZUMA1_Baseline", "MCPcounter")
run_correlation_analysis(z1_base, man_z1b, man_cols, "ZUMA1_Baseline", "ManualScores")

# ZUMA-7
run_correlation_analysis(z7_mrna, mcp_z7, mcp_cols, "ZUMA7_mRNA", "MCPcounter")
run_correlation_analysis(z7_mrna, man_z7, man_cols, "ZUMA7_mRNA", "ManualScores")

cat("\n==================================================================\n")
cat("Análise concluída! Verifique a pasta 'results/cellular_correlation'\n")
cat("==================================================================\n")
