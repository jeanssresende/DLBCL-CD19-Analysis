################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 07 - Differential Expression (Limma) e ssGSEA Per-Sample
#
# Objetivo:
#   1. Construir pareamento do ZUMA-1 (Baseline vs Post) dinamicamente
#   2. DEA clássico (Z1 e Z7) e Pareado
#   3. DEA em assinaturas NanoString (ZUMA-7)
#   4. ssGSEA (fgseaMultilevel) por amostra com expressão centralizada
#   5. Dotplot de categorias biológicas baseado no mean NES
################################################################################

library(tidyverse)
library(limma)
library(ggrepel)
library(patchwork)
library(fgsea)
library(msigdbr)
library(RColorBrewer)

dir.create("results/dea_gsea/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/dea_gsea/figures", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 1. Carregar Dados Purificados e Construir Pareamento (ZUMA-1)
## =============================================================================
message("Carregando matrizes purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
z7_sig  <- readRDS("data/filtered/ZUMA7_AxiCel_Sigs_Clean.rds") # Corrigido pelo print

message("Construindo coorte pareada (Baseline vs Post) do ZUMA-1...")
if(file.exists("data/filtered/ZUMA1_Post_Clean.rds")) {
  z1_post <- readRDS("data/filtered/ZUMA1_Post_Clean.rds")
  
  # Extrair IDs dos pacientes usando Regex na coluna 'title'
  clin_base <- z1_base$clin %>% rownames_to_column("sample") %>% 
    mutate(time_point = "Baseline", patient_id = str_extract(title, "Patient \\d+"))
  
  clin_post <- z1_post$clin %>% rownames_to_column("sample") %>% 
    mutate(time_point = "Post", patient_id = str_extract(title, "Patient \\d+"))
  
  # Intersecção: Apenas pacientes que têm ambas as biópsias
  common_patients <- intersect(clin_base$patient_id, clin_post$patient_id)
  common_patients <- common_patients[!is.na(common_patients)]
  
  clin_paired <- bind_rows(
    clin_base %>% filter(patient_id %in% common_patients),
    clin_post %>% filter(patient_id %in% common_patients)
  ) %>% column_to_rownames("sample")
  
  # Montar a matriz de expressão pareada
  samples_paired <- rownames(clin_paired)
  expr_paired <- cbind(
    z1_base$expr[, intersect(samples_paired, colnames(z1_base$expr))],
    z1_post$expr[, intersect(samples_paired, colnames(z1_post$expr))]
  )
  expr_paired <- expr_paired[, rownames(clin_paired)] # Ordenar
  
  z1_paired <- list(expr = expr_paired, clin = clin_paired)
  message("✅ ZUMA-1 Pareado construído! Pacientes: ", length(common_patients), " (Total amostras: ", ncol(expr_paired), ")")
} else {
  z1_paired <- NULL
  message("⚠️ Arquivo ZUMA1_Post_Clean.rds não encontrado. Pulando análise pareada.")
}

## =============================================================================
## 2. Funções Avançadas de DEA (Limma) - USANDO P-VALOR BRUTO
## =============================================================================
run_dea_limma <- function(expr, clin, group_col, group_levels, contrast_name) {
  
  meta <- clin %>% rownames_to_column("sample") %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Unknown")
  meta[[group_col]] <- factor(meta[[group_col]], levels = group_levels)
  meta <- meta %>% filter(!is.na(.data[[group_col]]))
  
  valid_samples <- intersect(meta$sample, colnames(expr))
  meta <- meta %>% filter(sample %in% valid_samples)
  expr_sub <- expr[, meta$sample, drop = FALSE]
  
  design <- model.matrix(~ 0 + meta[[group_col]])
  colnames(design) <- make.names(group_levels)
  
  contrast_formula <- paste0(colnames(design)[2], " - ", colnames(design)[1])
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  fit <- lmFit(expr_sub, design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("Feature") %>%
    mutate(
      contrast = contrast_name,
      comparison = paste0(group_levels[2], " vs ", group_levels[1]),
      # ALTERAÇÃO: Usando P.Value bruto no lugar de adj.P.Val
      direction_RawP = case_when(
        P.Value < 0.05 & logFC > 0.5  ~ paste0("Up in ", group_levels[2]),
        P.Value < 0.05 & logFC < -0.5 ~ paste0("Up in ", group_levels[1]),
        TRUE ~ "Not significant"
      )
    )
  
  write_csv(res, paste0("results/dea_gsea/tables/DEA_", contrast_name, "_RawP.csv"))
  return(res)
}

run_paired_dea_limma <- function(dataset, contrast_name) {
  if(is.null(dataset)) return(NULL)
  
  meta <- dataset$clin %>% rownames_to_column("sample") %>% filter(!is.na(patient_id), !is.na(time_point))
  meta$patient_id <- factor(meta$patient_id)
  meta$time_point <- factor(meta$time_point, levels = c("Baseline", "Post"))
  
  valid_samples <- intersect(meta$sample, colnames(dataset$expr))
  meta <- meta %>% filter(sample %in% valid_samples)
  expr_sub <- dataset$expr[, meta$sample, drop = FALSE]
  
  design <- model.matrix(~ patient_id + time_point, data = meta)
  
  fit <- lmFit(expr_sub, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = "time_pointPost", number = Inf, sort.by = "P") %>%
    rownames_to_column("Feature") %>%
    mutate(
      contrast = contrast_name, 
      comparison = "Post vs Baseline",
      # ALTERAÇÃO: Usando P.Value bruto
      direction_RawP = case_when(
        P.Value < 0.05 & logFC > 0.5  ~ "Up in Post",
        P.Value < 0.05 & logFC < -0.5 ~ "Up in Baseline",
        TRUE ~ "Not significant"
      )
    )
  
  write_csv(res, paste0("results/dea_gsea/tables/DEA_", contrast_name, "_RawP.csv"))
  return(res)
}

## =============================================================================
## 3. Execução da DEA e Volcano Plots (Baseado em P-valor Bruto)
## =============================================================================
message("Rodando Modelos Diferenciais Limma...")

dea_z1 <- run_dea_limma(z1_base$expr, z1_base$clin, "response_group", c("CR", "Non-CR"), "ZUMA1_Baseline_NonCR_vs_CR")
dea_z7 <- run_dea_limma(z7_mrna$expr, z7_mrna$clin, "response_group", c("Ongoing", "Others"), "ZUMA7_AxiCel_Others_vs_Ongoing")
dea_z7_sigs <- run_dea_limma(z7_sig$expr, z7_sig$clin, "response_group", c("Ongoing", "Others"), "ZUMA7_Signatures_Others_vs_Ongoing")
dea_z1_paired <- run_paired_dea_limma(z1_paired, "ZUMA1_Paired_Post_vs_Baseline")

plot_volcano <- function(dea_df, title, logfc_cutoff = 0.5, p_cutoff = 0.05) {
  if(is.null(dea_df)) return(NULL)
  
  plot_df <- dea_df %>%
    mutate(
      # ALTERAÇÃO: Eixo Y agora é o log10 do P.Value bruto
      neg_log10_p = -log10(P.Value),
      sig_status = ifelse(is.na(direction_RawP), "Not significant", direction_RawP)
    )
  
  # Seleciona os top 15 genes com menor p-valor bruto para colocar o nome
  top_genes <- plot_df %>% filter(sig_status != "Not significant") %>% arrange(P.Value) %>% slice_head(n = 15)
  
  ggplot(plot_df, aes(x = logFC, y = neg_log10_p)) +
    geom_point(aes(color = sig_status), alpha = 0.7, size = 1.5) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color="gray40") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color="gray40") +
    geom_text_repel(data = top_genes, aes(label = Feature), size = 3.5, fontface="bold", box.padding = 0.5, max.overlaps = Inf) +
    scale_color_manual(values = c(
      "Up in CR" = "#4DBBD5FF", "Up in Non-CR" = "#E64B35FF",
      "Up in Ongoing" = "#00A087FF", "Up in Others" = "#3C5488FF",
      "Up in Baseline" = "#4DBBD5FF", "Up in Post" = "#E64B35FF",
      "Not significant" = "grey80"
    )) +
    theme_classic(base_size = 14) +
    # ALTERAÇÃO: Títulos dos eixos e legenda atualizados para refletir o Raw P
    labs(title = title, subtitle = "Raw P-value < 0.05 and |log2FC| > 0.5", x = "Log2 Fold Change", y = "-Log10(Raw P-value)") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom", legend.title = element_blank())
}

p_vol_z1 <- plot_volcano(dea_z1, "ZUMA-1 Baseline: Non-CR vs CR")
p_vol_z7 <- plot_volcano(dea_z7, "ZUMA-7 Axi-cel: Others vs Ongoing")
p_vol_z1_paired <- plot_volcano(dea_z1_paired, "ZUMA-1 Paired: Post vs Baseline")

if(!is.null(p_vol_z1)) ggsave("results/dea_gsea/figures/Volcano_ZUMA1_RawP.pdf", p_vol_z1, width = 6, height = 6)
if(!is.null(p_vol_z7)) ggsave("results/dea_gsea/figures/Volcano_ZUMA7_RawP.pdf", p_vol_z7, width = 6, height = 6)
if(!is.null(p_vol_z1_paired)) ggsave("results/dea_gsea/figures/Volcano_ZUMA1_Paired_RawP.pdf", p_vol_z1_paired, width = 6, height = 6)

## =============================================================================
## 4. Dicionário de Categorias Funcionais (Hallmarks MSigDB)
## =============================================================================
cat_mapping <- list(
  "Metabolic" = c("HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_HEME_METABOLISM", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_GLYCOLYSIS", "HALLMARK_XENOBIOTIC_METABOLISM", "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
  "Proliferation" = c("HALLMARK_G2M_CHECKPOINT", "HALLMARK_E2F_TARGETS", "HALLMARK_MITOTIC_SPINDLE", "HALLMARK_P53_PATHWAY", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYC_TARGETS_V1"),
  "Immune" = c("HALLMARK_COMPLEMENT", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_COAGULATION"),
  "Signaling" = c("HALLMARK_NOTCH_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_HEDGEHOG_SIGNALING", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_MTORC1_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
  "Pathway" = c("HALLMARK_APOPTOSIS", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_PROTEIN_SECRETION", "HALLMARK_HYPOXIA"),
  "DNA damage" = c("HALLMARK_UV_RESPONSE_DN", "HALLMARK_DNA_REPAIR", "HALLMARK_UV_RESPONSE_UP"),
  "Cellular component" = c("HALLMARK_APICAL_SURFACE", "HALLMARK_APICAL_JUNCTION", "HALLMARK_PEROXISOME"),
  "Development" = c("HALLMARK_PANCREAS_BETA_CELLS", "HALLMARK_ANGIOGENESIS", "HALLMARK_SPERMATOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ADIPOGENESIS", "HALLMARK_MYOGENESIS")
)

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

## =============================================================================
## 5. Módulo ssGSEA (Per Sample fgseaMultilevel) - SEM FILTRO DE P-VALOR
## =============================================================================
message("Rodando ssGSEA por amostra...")

run_ssgsea_per_sample <- function(dataset, group_col, cohort_name) {
  if(is.null(dataset)) return(NULL)
  
  expr <- dataset$expr
  meta <- dataset$clin %>% rownames_to_column("sample") %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] != "Unknown")
  
  # Centralização por Gene
  expr_centered <- expr - rowMeans(expr, na.rm = TRUE)
  
  res_list <- list()
  samples <- meta$sample[meta$sample %in% colnames(expr_centered)]
  
  # Cálculo indivídual
  for(s in samples) {
    gene_ranks <- sort(expr_centered[, s], decreasing = TRUE)
    set.seed(42)
    fgsea_res <- fgseaMultilevel(pathways = hallmark_list, stats = gene_ranks, minSize = 10, maxSize = 500)
    res_list[[s]] <- fgsea_res %>% mutate(sample = s)
  }
  
  all_ssgsea <- bind_rows(res_list) %>% left_join(meta %>% select(sample, Group = all_of(group_col)), by = "sample")
  
  # ALTERAÇÃO: Removido o p-valor. Criada a métrica de Magnitude Absoluta (Abs_Mean_NES)
  summary_df <- all_ssgsea %>%
    group_by(Group, pathway) %>%
    summarise(
      Mean_NES = mean(NES, na.rm = TRUE),
      Abs_Mean_NES = abs(mean(NES, na.rm = TRUE)), # Nova métrica para o tamanho da bolinha
      .groups = "drop"
    ) %>%
    mutate(Cohort = cohort_name)
  
  return(list(raw = all_ssgsea, summary = summary_df))
}

ssgsea_z1 <- run_ssgsea_per_sample(z1_base, "response_group", "ZUMA-1 (Discovery)")
ssgsea_z7 <- run_ssgsea_per_sample(z7_mrna, "response_group", "ZUMA-7 (Validation)")

## =============================================================================
## 6. Geração do Dotplot Funcional Categorizado (CORES E LEGENDA CORRIGIDAS)
## =============================================================================
build_dotplot <- function(ssgsea_summary1, ssgsea_summary2) {
  
  combined <- bind_rows(ssgsea_summary1, ssgsea_summary2)
  
  pathway_cat_map <- unlist(lapply(names(cat_mapping), function(cat) { setNames(rep(cat, length(cat_mapping[[cat]])), cat_mapping[[cat]]) }))
  
  plot_df <- combined %>%
    mutate(
      Category = pathway_cat_map[pathway],
      Category = factor(Category, levels = names(cat_mapping)),
      Pathway_Clean = str_remove(pathway, "HALLMARK_")
    ) %>%
    filter(!is.na(Category))
  
  plot_df$Pathway_Clean <- factor(plot_df$Pathway_Clean, levels = rev(unique(plot_df$Pathway_Clean[order(plot_df$Category)])))
  
  p <- ggplot(plot_df, aes(x = Group, y = Pathway_Clean)) +
    geom_point(aes(size = Abs_Mean_NES, fill = Mean_NES), shape = 21, color = "gray30", stroke = 0.5) +
    facet_grid(Category ~ Cohort, scales = "free", space = "free") +
    # CORREÇÃO 1: Escala travada entre -1 e 1 para cores mais vivas
    scale_fill_gradientn(colors = rev(brewer.pal(7, "RdBu")), limits = c(-1, 1), oob = scales::squish, name = "Mean NES") +
    # CORREÇÃO 2: Forçando a legenda de tamanho a aparecer corretamente
    scale_size_continuous(range = c(2, 7), name = "Magnitude\n|Mean NES|", 
                          guide = guide_legend(override.aes = list(fill = "gray70"))) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = NULL, title = "Pathway Enrichment Profiles across Clinical Responses") +
    theme(
      strip.text.y.right = element_text(angle = 0, face = "bold", size = 10),
      strip.text.x = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = "white"),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.spacing = unit(0.1, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  return(p)
}

fig_dotplot <- build_dotplot(ssgsea_z1$summary, ssgsea_z7$summary)

ggsave("results/dea_gsea/figures/Functional_Dotplot_Z1_Z7.pdf", fig_dotplot, width = 12, height = 14)
ggsave("results/dea_gsea/figures/Functional_Dotplot_Z1_Z7.png", fig_dotplot, width = 12, height = 14, dpi = 300)

## =============================================================================
## 7. GSEA Clássico e Enrichment Plots Individuais
## =============================================================================
message("Rodando GSEA Clássico e gerando Enrichment Plots...")

dir.create("results/dea_gsea/figures/enrichment_plots", showWarnings = FALSE, recursive = TRUE)

run_and_plot_classic_gsea <- function(dea_res, cohort_name) {
  if(is.null(dea_res)) return(NULL)
  
  # 1. Ranqueando os genes pela estatística T do Limma (Força da diferença entre os grupos)
  ranks <- dea_res %>%
    filter(!is.na(t), !is.na(Feature)) %>%
    arrange(desc(t)) %>%
    distinct(Feature, .keep_all = TRUE) %>%
    { setNames(.$t, .$Feature) }
  
  # 2. Rodando o GSEA Clássico
  set.seed(42)
  fgsea_res <- fgsea(pathways = hallmark_list, stats = ranks, minSize = 10, maxSize = 500)
  
  # Salvando a tabela do GSEA Clássico
  export_df <- fgsea_res %>% 
    arrange(pval) %>% 
    as_tibble() %>% 
    mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ", ", FUN.VALUE = character(1)))
  
  write_csv(export_df, paste0("results/dea_gsea/tables/Classic_GSEA_", cohort_name, ".csv"))
  
  # 3. Gerando o Gráfico de Montanha (Enrichment Plot) para as vias do nosso dicionário
  pathways_to_plot <- unlist(cat_mapping, use.names = FALSE)
  
  for (pw in pathways_to_plot) {
    pw_res <- fgsea_res %>% filter(pathway == pw)
    if(nrow(pw_res) == 0) next # Pula se a via não foi calculada
    
    # Criando o plot
    p <- plotEnrichment(hallmark_list[[pw]], ranks) +
      labs(
        title = str_remove(pw, "HALLMARK_"),
        subtitle = paste0(cohort_name, " | NES: ", round(pw_res$NES, 2), " | Raw p-val: ", signif(pw_res$pval, 3))
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(color = "gray30", size = 10)
      )
    
    # Salvando em PDF
    safe_pw <- str_remove(pw, "HALLMARK_")
    ggsave(paste0("results/dea_gsea/figures/enrichment_plots/", cohort_name, "_", safe_pw, ".pdf"), 
           p, width = 6, height = 4)
  }
}

run_and_plot_classic_gsea(dea_z1, "ZUMA1_Baseline")
run_and_plot_classic_gsea(dea_z7, "ZUMA7_AxiCel")
run_and_plot_classic_gsea(dea_z1_paired, "ZUMA1_Paired")

message("\n==========================================================")
message("Análises Completas! Verifique a pasta 'enrichment_plots'.")
message("==========================================================")
