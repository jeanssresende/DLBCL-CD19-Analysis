################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 08 - Gene Ontology (GO), KEGG Clássico e Pathview Maps
#
# Objetivo:
#   1. Extrair os DEGs (Genes Diferencialmente Expressos) do Limma
#   2. Rodar Enriquecimento Clássico (ORA) para GO (BP) e KEGG
#   3. Gerar Dotplots de alta qualidade (clusterProfiler)
#   4. Renderizar os mapas KEGG (Pathview) para vias ESPECÍFICAS escolhidas
################################################################################

# ==============================================================================
# 1. Carregar Pacotes Necessários
# ==============================================================================
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)

dir.create("results/go_kegg/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/go_kegg/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/go_kegg/figures/pathview_maps", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2. Carregar Tabelas de DEA (Limma) e Filtrar Significativos
# ==============================================================================
message("Carregando resultados do Limma (Script 07) e filtrando DEGs...")

dea_z1 <- read_csv("results/dea_gsea/tables/DEA_ZUMA1_Baseline_NonCR_vs_CR_RawP.csv")
dea_z7 <- read_csv("results/dea_gsea/tables/DEA_ZUMA7_AxiCel_Others_vs_Ongoing_RawP.csv")
dea_z1_paired <- read_csv("results/dea_gsea/tables/DEA_ZUMA1_Paired_Post_vs_Baseline_RawP.csv")

prepare_ora_lists <- function(dea_df) {
  # Filtra apenas os genes significativos (p < 0.05 e |logFC| > 0.5)
  degs <- dea_df %>% 
    filter(P.Value < 0.05, abs(logFC) > 0.5, !is.na(Feature))
  
  # Traduzir Gene Symbols para Entrez IDs (Necessário para KEGG)
  entrez_map <- bitr(degs$Feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # Preparar vetor de LogFC de TODOS os genes para o Pathview colorir o background
  logfc_all <- dea_df$logFC
  entrez_all <- bitr(dea_df$Feature, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  df_all <- dea_df %>% inner_join(entrez_all, by = c("Feature" = "SYMBOL")) %>% distinct(ENTREZID, .keep_all = TRUE)
  
  logfc_vector <- df_all$logFC
  names(logfc_vector) <- as.character(df_all$ENTREZID)
  
  return(list(
    symbol = degs$Feature,           # Usado no GO
    entrez = entrez_map$ENTREZID,    # Usado no KEGG
    logfc_vector = logfc_vector      # Usado no Pathview
  ))
}

lists_z1 <- prepare_ora_lists(dea_z1)
lists_z7 <- prepare_ora_lists(dea_z7)
lists_z1_paired <- prepare_ora_lists(dea_z1_paired)

# ==============================================================================
# 3. Executar o Enriquecimento Clássico: Gene Ontology (GO - BP)
# ==============================================================================
message("Rodando Enriquecimento para Gene Ontology (Biological Process)...")

run_go_ora <- function(gene_symbols, cohort_name) {
  if(length(gene_symbols) == 0) return(NULL)
  
  ego <- enrichGO(
    gene          = gene_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP", 
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  if(!is.null(ego) && nrow(as_tibble(ego)) > 0) {
    write_csv(as_tibble(ego), paste0("results/go_kegg/tables/GO_ORA_", cohort_name, ".csv"))
  }
  return(ego)
}

go_z1 <- run_go_ora(lists_z1$symbol, "ZUMA1_Baseline")
go_z7 <- run_go_ora(lists_z7$symbol, "ZUMA7_Baseline")
go_z1_paired <- run_go_ora(lists_z1_paired$symbol, "ZUMA1_Paired")

# ==============================================================================
# 4. Executar o Enriquecimento Clássico: KEGG Pathways
# ==============================================================================
message("Rodando Enriquecimento para vias KEGG...")

run_kegg_ora <- function(gene_entrez, cohort_name) {
  if(length(gene_entrez) == 0) return(NULL)
  
  ekegg <- enrichKEGG(
    gene         = gene_entrez,
    organism     = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if(!is.null(ekegg) && nrow(as_tibble(ekegg)) > 0) {
    ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    write_csv(as_tibble(ekegg_readable), paste0("results/go_kegg/tables/KEGG_ORA_", cohort_name, ".csv"))
  }
  return(ekegg)
}

kegg_z1 <- run_kegg_ora(lists_z1$entrez, "ZUMA1_Baseline")
kegg_z7 <- run_kegg_ora(lists_z7$entrez, "ZUMA7_Baseline")
kegg_z1_paired <- run_kegg_ora(lists_z1_paired$entrez, "ZUMA1_Paired")

# ==============================================================================
# 5. Gerar Dotplots Globais
# ==============================================================================
message("Gerando Dotplots...")

generate_dotplot <- function(enrich_obj, title, filename) {
  if(is.null(enrich_obj) || nrow(as_tibble(enrich_obj)) == 0) {
    message("Sem vias significativas para: ", title)
    return(NULL)
  }
  
  p <- dotplot(enrich_obj, showCategory = 15, font.size = 11) +
    theme_bw() +
    labs(title = title) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(paste0("results/go_kegg/figures/", filename, ".pdf"), p, width = 9, height = 7)
  ggsave(paste0("results/go_kegg/figures/", filename, ".png"), p, width = 9, height = 7, dpi = 300)
}

generate_dotplot(go_z1, "ZUMA-1 Baseline: Top GO BP", "Dotplot_GO_ZUMA1_Baseline")
generate_dotplot(go_z7, "ZUMA-7 Baseline: Top GO BP", "Dotplot_GO_ZUMA7_Baseline")
generate_dotplot(go_z1_paired, "ZUMA-1 Paired: Top GO BP", "Dotplot_GO_ZUMA1_Paired")

generate_dotplot(kegg_z1, "ZUMA-1 Baseline: Top KEGG", "Dotplot_KEGG_ZUMA1_Baseline")
generate_dotplot(kegg_z7, "ZUMA-7 Baseline: Top KEGG", "Dotplot_KEGG_ZUMA7_Baseline")
generate_dotplot(kegg_z1_paired, "ZUMA-1 Paired: Top KEGG", "Dotplot_KEGG_ZUMA1_Paired")

# ==============================================================================
# 6. MAPAS KEGG COM PATHVIEW (VIAS ESCOLHIDAS PELO USUÁRIO)
# ==============================================================================
message("Renderizando mapas do KEGG (Pathview) para as vias selecionadas...")

current_dir <- getwd()
setwd(paste0(current_dir, "/results/go_kegg/figures/pathview_maps"))

# ------------------------------------------------------------------------------
# INSIRA AQUI OS IDs DAS VIAS KEGG QUE VOCÊ QUER GERAR OS MAPAS
# (Você pode olhar o arquivo CSV do KEGG gerado na pasta tables para escolher)
# Exemplo: "hsa04630" é a via JAK-STAT, "hsa04060" é Cytokine-cytokine receptor
# ------------------------------------------------------------------------------
vias_z1_baseline <- c("hsa04060", "hsa04630") 
vias_z7_baseline <- c("hsa04060")             
vias_z1_paired   <- c("hsa04611", "hsa04630") 

run_pathview_custom <- function(pathways_list, logfc_vector, cohort_name) {
  if(length(pathways_list) == 0) return(NULL)
  
  for (path_id in pathways_list) {
    message(paste("Renderizando mapa:", path_id, "para", cohort_name))
    tryCatch({
      pathview(
        gene.data  = logfc_vector,
        pathway.id = path_id,
        species    = "hsa",
        limit      = list(gene = max(abs(logfc_vector), na.rm = TRUE), cpd = 1),
        low        = list(gene = "#313695", cpd = "blue"),  # Azul = Down
        mid        = list(gene = "gray90", cpd = "gray"),
        high       = list(gene = "#A50026", cpd = "yellow"),# Vermelho = Up
        kegg.native = TRUE,
        out.suffix = cohort_name
      )
    }, error = function(e) {
      message(paste("Erro no download do mapa", path_id))
    })
  }
}

run_pathview_custom(vias_z1_baseline, lists_z1$logfc_vector, "ZUMA1_Baseline")
run_pathview_custom(vias_z7_baseline, lists_z7$logfc_vector, "ZUMA7_Baseline")
run_pathview_custom(vias_z1_paired, lists_z1_paired$logfc_vector, "ZUMA1_Paired")

setwd(current_dir)

message("\n==========================================================")
message("Análise ORA (GO/KEGG) e Pathview concluída com sucesso!")
message("==========================================================")

