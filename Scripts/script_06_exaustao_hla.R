################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19 (ZUMA-1 vs ZUMA-7)
# Script: 06 - Escudos Tumorais (Core Exhaustion e HLA) e Heatmaps de Alto Impacto
# Objetivo:
#   1. Estatística e Violinos (FDR INDEPENDENTE POR HIPÓTESE)
#   2. Metadados Clínicos Completos (Response original, COO, CRS, NE, EFS)
#   3. Heatmaps Complexos (Z1 Base, Z7 Axi, Z7 SOC, Z1 Pareado)
################################################################################

## =============================================================================
## 1. Pacotes e Diretórios
## =============================================================================
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
dir.create("results/exhaustion_hla/single_plots", showWarnings = FALSE, recursive = TRUE)

## =============================================================================
## 2. Carregar Dados Purificados
## =============================================================================
message("Carregando matrizes purificadas...")
z1_base <- readRDS("data/filtered/ZUMA1_Baseline_Clean.rds")
z1_post <- readRDS("data/filtered/ZUMA1_Post_Clean.rds")
z7_mrna <- readRDS("data/filtered/ZUMA7_AxiCel_mRNA_Clean.rds")
z7_soc  <- readRDS("data/filtered/ZUMA7_SOC_mRNA_Clean.rds")
z7_sig  <- readRDS("data/filtered/ZUMA7_AxiCel_Sigs_Clean.rds")

# Padronização ZUMA-7 SOC
if("ongoing_2grps_ch1" %in% colnames(z7_soc$clin) && !"response_group" %in% colnames(z7_soc$clin)) {
  z7_soc$clin$response_group <- z7_soc$clin$ongoing_2grps_ch1
}

# Paletas Originais (Mantidas para consistência)
pal_binaria <- c("CR" = "#4DBBD5", "Non-CR" = "#E64B35", "Ongoing" = "#00A087", "Others" = "#3C5488")
pal_granular <- c("CR" = "#4DBBD5", "PR" = "#F39B7F", "SD" = "#7E6148", "PD" = "#E64B35",
                  "Ongoing Response" = "#00A087", "Relapsed" = "#8491B4", "Nonresponders" = "#3C5488")
pal_tempo_cr <- c("Baseline" = "#4DBBD5", "Post-CART" = "#0072B5")

# Paleta para Subgrupos Moleculares / Cell of Origin
pal_coo <- c("GCB" = "#B3E2CD", "ABC" = "#FDCDAC", "UNCLASSIFIED" = "#CBD5E8", "Unclassified" = "#CBD5E8")

## =============================================================================
## 3. Definição dos Grupos de Genes
## =============================================================================
exhaustion_signature_groups <- list(
  "Core Exhaustion" = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "EOMES", "ENTPD1", "CXCL13")
)

hla_signature_groups <- list(
  "Classical HLA-I" = c("HLA-A", "HLA-B", "HLA-C"),
  "HLA-I processing" = c("B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5"),
  "Classical HLA-II" = c("HLA-DRA", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
  "Non-classical HLA-I" = c("HLA-E", "HLA-F", "HLA-G")
)

genes_exhaustion <- unlist(exhaustion_signature_groups, use.names = FALSE)
genes_hla <- unlist(hla_signature_groups, use.names = FALSE)
all_target_genes <- unique(c(genes_exhaustion, genes_hla))

## =============================================================================
## 4. Integração de Expressão (IO360 e ZUMA-1 Pareado)
## =============================================================================
io360_checkpoints <- c("CTLA4.IO360", "PD-L1.IO360", "PD-L2.IO360", "TIGIT.IO360", "Treg.IO360", "Exhausted CD8.IO360")
valid_io360 <- intersect(io360_checkpoints, rownames(z7_sig$expr))
expr_z7_axi_combined <- rbind(z7_mrna$expr, z7_sig$expr[valid_io360, , drop = FALSE])
genes_z7_exhaustion <- unique(c(genes_exhaustion, valid_io360))

build_plot_df <- function(dataset_list, target_genes) {
  expr <- dataset_list$expr
  clin <- dataset_list$clin %>% rownames_to_column("sample")
  if("title" %in% colnames(clin)) clin <- clin %>% mutate(patient_id = str_extract(title, "Patient \\d+"))
  genes_present <- intersect(target_genes, rownames(expr))
  df_genes <- t(expr[genes_present, , drop = FALSE]) %>% as.data.frame() %>% rownames_to_column("sample")
  return(clin %>% left_join(df_genes, by = "sample"))
}

df_z1_base <- build_plot_df(z1_base, all_target_genes)
df_z1_post <- build_plot_df(z1_post, all_target_genes)
df_z7_axi  <- build_plot_df(list(expr = expr_z7_axi_combined, clin = z7_mrna$clin), c(all_target_genes, valid_io360))
df_z7_soc  <- build_plot_df(z7_soc, all_target_genes)

## =============================================================================
## 5. Estatísticas e Violinos (Análise Individual Gene a Gene)
## =============================================================================
message("\n--- Iniciando Análise Estatística e Plotagem de Violinos ---")

run_unpaired_and_plot <- function(df, features, analysis_name) {
  valid_features <- intersect(features, colnames(df))
  if(length(valid_features) == 0) return(NULL)
  
  if (!"response_original" %in% colnames(df)) df$response_original <- "Missing"
  
  stats <- map_dfr(valid_features, function(feature) {
    tmp <- df %>% select(group = response_group, value = all_of(feature)) %>% 
      filter(!group %in% c("Unknown", "Missing", NA)) %>% drop_na()
    if (nlevels(as.factor(tmp$group)) != 2) return(NULL)
    test <- wilcox.test(value ~ as.factor(group), data = tmp)
    tibble(feature = feature, group_1 = levels(as.factor(tmp$group))[1], group_2 = levels(as.factor(tmp$group))[2], p_value = test$p.value)
  }) %>% mutate(padj_BH = p.adjust(p_value, method = "BH"), label_padj = paste0("FDR = ", signif(padj_BH, 2)))
  
  write_csv(stats, paste0("results/exhaustion_hla/", analysis_name, "_Unpaired_Stats.csv"))
  
  walk(valid_features, function(feature_name) {
    plot_df <- df %>% select(sample, group = response_group, response_original, score = all_of(feature_name)) %>%
      filter(!group %in% c("Unknown", "Missing", NA)) %>% drop_na(score) %>%
      mutate(response_original = ifelse(is.na(response_original) | response_original == "", "Missing", response_original))
    stat_row <- stats %>% filter(feature == feature_name) %>% mutate(group1 = group_1, group2 = group_2)
    y_pos <- max(plot_df$score, na.rm = TRUE) + 0.10 * diff(range(plot_df$score, na.rm = TRUE))
    
    p <- ggplot(plot_df, aes(x = group, y = score)) +
      geom_violin(aes(fill = group), alpha = 0.3, color = "gray40", scale = "width") +
      geom_boxplot(aes(fill = group), width = 0.2, outlier.shape = NA, alpha = 0.7, color = "gray30") +
      geom_jitter(aes(color = response_original), width = 0.18, size = 2, alpha = 0.85) +
      scale_fill_manual(values = pal_binaria) + scale_color_manual(values = pal_granular) +
      theme_classic(base_size = 14) + labs(title = feature_name, subtitle = analysis_name, x = NULL, y = "Expression (Log2)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 10, color="gray40"), legend.position = "none")
    if(nrow(stat_row) > 0) p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, tip.length = 0.015)
    ggsave(paste0("results/exhaustion_hla/single_plots/", analysis_name, "_", gsub("[^A-Za-z0-9]", "_", feature_name), "_Violin.pdf"), p, width = 4.5, height = 5)
  })
}

run_paired_and_plot <- function(df_base, df_post, features, analysis_name, id_col = "patient_id") {
  valid_features <- intersect(features, intersect(colnames(df_base), colnames(df_post)))
  if(length(valid_features) == 0) return(NULL)
  
  stats <- map_dfr(valid_features, function(feature) {
    b_tmp <- df_base %>% select(all_of(id_col), base_val = all_of(feature))
    p_tmp <- df_post %>% select(all_of(id_col), post_val = all_of(feature))
    paired_df <- inner_join(b_tmp, p_tmp, by = id_col) %>% drop_na()
    if (nrow(paired_df) < 3) return(NULL)
    test <- wilcox.test(paired_df$post_val, paired_df$base_val, paired = TRUE)
    tibble(feature = feature, N_pairs = nrow(paired_df), p_value = test$p.value)
  }) %>% mutate(padj_BH = p.adjust(p_value, method = "BH"), label_padj = paste0("FDR = ", signif(padj_BH, 2)))
  
  write_csv(stats, paste0("results/exhaustion_hla/", analysis_name, "_Paired_Stats.csv"))
  
  walk(valid_features, function(feature_name) {
    b_tmp <- df_base %>% select(all_of(id_col), score = all_of(feature_name)) %>% mutate(Timepoint = "Baseline")
    p_tmp <- df_post %>% select(all_of(id_col), score = all_of(feature_name)) %>% mutate(Timepoint = "Post-CART")
    plot_df <- bind_rows(b_tmp, p_tmp) %>% mutate(Timepoint = factor(Timepoint, levels = c("Baseline", "Post-CART"))) %>% drop_na(score)
    complete_patients <- plot_df %>% count(!!sym(id_col)) %>% filter(n == 2) %>% pull(!!sym(id_col))
    plot_df <- plot_df %>% filter(!!sym(id_col) %in% complete_patients)
    
    stat_row <- stats %>% filter(feature == feature_name) %>% mutate(group1 = "Baseline", group2 = "Post-CART")
    y_pos <- max(plot_df$score, na.rm = TRUE) + 0.08 * diff(range(plot_df$score, na.rm = TRUE))
    
    p <- ggplot(plot_df, aes(x = Timepoint, y = score)) +
      geom_line(aes(group = !!sym(id_col)), color = "gray60", alpha = 0.5, linewidth = 0.8) +
      geom_boxplot(aes(fill = Timepoint), width = 0.2, alpha = 0.7, outlier.shape = NA, color = "black") +
      geom_point(aes(fill = Timepoint), shape = 21, size = 2.5, color = "white", alpha = 0.9) +
      scale_fill_manual(values = pal_tempo_cr) +
      theme_classic(base_size = 14) + labs(title = feature_name, subtitle = analysis_name, x = NULL, y = "Expression (Log2)") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 10, color="gray40"), legend.position = "none")
    if(nrow(stat_row) > 0) p <- p + stat_pvalue_manual(stat_row, label = "label_padj", y.position = y_pos, tip.length = 0.015)
    ggsave(paste0("results/exhaustion_hla/single_plots/", analysis_name, "_", gsub("[^A-Za-z0-9]", "_", feature_name), "_Paired.pdf"), p, width = 4.5, height = 5)
  })
}

message("-> Avaliando Core Exhaustion...")
run_unpaired_and_plot(df_z1_base, genes_exhaustion, "ZUMA1_Baseline_Exhaustion")
run_unpaired_and_plot(df_z7_axi, genes_z7_exhaustion, "ZUMA7_AxiCel_Exhaustion")
run_unpaired_and_plot(df_z7_soc, genes_exhaustion, "ZUMA7_SOC_Exhaustion")
run_paired_and_plot(df_z1_base, df_z1_post, genes_exhaustion, "ZUMA1_Paired_Exhaustion", id_col = "patient_id")

message("-> Avaliando Maquinaria HLA...")
run_unpaired_and_plot(df_z1_base, genes_hla, "ZUMA1_Baseline_HLA")
run_unpaired_and_plot(df_z7_axi, genes_hla, "ZUMA7_AxiCel_HLA")
run_unpaired_and_plot(df_z7_soc, genes_hla, "ZUMA7_SOC_HLA")
run_paired_and_plot(df_z1_base, df_z1_post, genes_hla, "ZUMA1_Paired_HLA", id_col = "patient_id")


## =============================================================================
## 6. Preparação dos Metadados (COMPLETOS E BLINDADOS CONTRA ERROS)
## =============================================================================
message("\n--- Preparando Metadados Clínicos Completos para Heatmaps ---")

# Função de Limpeza Universal de Metadados (Versão Corrigida para Numéricos)
clean_anno_full <- function(df, resp_col, class_col=NULL, coo_col=NULL, crs_col=NULL, ne_col=NULL, efs_col=NULL, tb_col=NULL) {
  res <- df %>% rownames_to_column("sample") %>%
    mutate(Response = as.character(!!sym(resp_col))) %>%
    filter(!Response %in% c("Unknown", "Missing", NA, "")) 
  
  if(!is.null(class_col) && class_col %in% colnames(res)) res$Clinical_Class <- as.character(res[[class_col]])
  if(!is.null(coo_col) && coo_col %in% colnames(res)) res$Molecular_Subgroup <- as.character(res[[coo_col]])
  if(!is.null(efs_col) && efs_col %in% colnames(res)) res$EFS_Event <- as.character(res[[efs_col]])
  if(!is.null(crs_col) && crs_col %in% colnames(res)) res$Grade3_CRS <- as.character(res[[crs_col]])
  
  # 1. Tratar variáveis categóricas (Substituir NAs por "Missing")
  cat_cols <- c("Response", "Clinical_Class", "Molecular_Subgroup", "EFS_Event", "Grade3_CRS")
  for(col in intersect(cat_cols, colnames(res))) {
    res[[col]][is.na(res[[col]]) | res[[col]] == ""] <- "Missing"
  }
  
  # 2. Tratar variáveis numéricas (MANTER COMO NUMÉRICO para o colorRamp2)
  if(!is.null(ne_col) && ne_col %in% colnames(res)) {
    if(grepl("grade3", ne_col, ignore.case=TRUE)) {
      res$Grade3_NE <- as.character(res[[ne_col]])
      res$Grade3_NE[is.na(res$Grade3_NE) | res$Grade3_NE == ""] <- "Missing"
    } else {
      res$Neurotox_Grade <- as.numeric(str_extract(res[[ne_col]], "\\d+"))
    }
  }
  
  if(!is.null(tb_col) && tb_col %in% colnames(res)) {
    res$Tumor_Burden_Log10 <- log10(as.numeric(res[[tb_col]]) + 1)
  }
  
  return(res %>% column_to_rownames("sample"))
}

# Extraindo todas as colunas pedidas
anno_z1 <- clean_anno_full(z1_base$clin, "response_group", class_col="response_original", coo_col="molecular_subgroup_ch1", ne_col="worst_grade_of_ne_ch1", tb_col="baseline_tumour_burden_mm2_ch1") %>% select(any_of(c("Response", "Clinical_Class", "Molecular_Subgroup", "Neurotox_Grade", "Tumor_Burden_Log10")))

anno_z7_axi <- clean_anno_full(z7_mrna$clin, "response_group", class_col="response_original", coo_col="cell_of_origin_ch1", crs_col="grade3_crs_ch1", ne_col="grade3_ne_ch1", efs_col="event_free_survival_event_ch1", tb_col="baseline_tumor_burden_spd_ch1") %>% select(any_of(c("Response", "Clinical_Class", "Molecular_Subgroup", "Grade3_CRS", "Grade3_NE", "EFS_Event", "Tumor_Burden_Log10")))

anno_z7_soc <- clean_anno_full(z7_soc$clin, "response_group", coo_col="cell_of_origin_ch1", efs_col="event_free_survival_event_ch1", tb_col="baseline_tumor_burden_spd_ch1") %>% select(any_of(c("Response", "Molecular_Subgroup", "EFS_Event", "Tumor_Burden_Log10")))

# ZUMA-1 Pareado
p_base <- z1_base$clin %>% rownames_to_column("sample") %>% mutate(patient_id = str_extract(title, "Patient \\d+"))
p_post <- z1_post$clin %>% rownames_to_column("sample") %>% mutate(patient_id = str_extract(title, "Patient \\d+"))
valid_pairs <- intersect(p_base$patient_id, p_post$patient_id)

anno_z1_paired <- bind_rows(
  p_base %>% filter(patient_id %in% valid_pairs) %>% mutate(Timepoint = "Baseline"),
  p_post %>% filter(patient_id %in% valid_pairs) %>% mutate(Timepoint = "Post-CART")
) %>% mutate(Response = as.character(response_group)) %>% select(sample, Timepoint, Response) %>% column_to_rownames("sample")

## --- Função de Segurança para Paletas de Cores ---
make_safe_pal <- function(x, base_pal) {
  lvls <- unique(as.character(x))
  out <- setNames(rep("gray80", length(lvls)), lvls)
  out["Missing"] <- "white"
  for(l in lvls) {
    if(l %in% names(base_pal)) out[l] <- base_pal[l]
    if(l %in% c("0", "No", "FALSE")) out[l] <- "gray90"
    if(l %in% c("1", "Yes", "TRUE")) out[l] <- "black"
  }
  return(out)
}

# Criando paletas numéricas (Agora o R sabe que é número puro, sem "Missing" atrapalhando)
col_tb_z1 <- colorRamp2(c(0, max(anno_z1$Tumor_Burden_Log10, na.rm=TRUE)), c("gray95", "black"))
col_ne_z1 <- colorRamp2(c(0, max(anno_z1$Neurotox_Grade, na.rm=TRUE)), c("gray95", "black"))
col_tb_z7 <- colorRamp2(c(0, max(anno_z7_axi$Tumor_Burden_Log10, na.rm=TRUE)), c("gray95", "black"))
col_tb_soc <- colorRamp2(c(0, max(anno_z7_soc$Tumor_Burden_Log10, na.rm=TRUE)), c("gray95", "black"))

ha_z1 <- HeatmapAnnotation(df = anno_z1, col = list(
  Response = make_safe_pal(anno_z1$Response, pal_binaria),
  Clinical_Class = make_safe_pal(anno_z1$Clinical_Class, pal_granular),
  Molecular_Subgroup = make_safe_pal(anno_z1$Molecular_Subgroup, pal_coo),
  Neurotox_Grade = col_ne_z1, Tumor_Burden_Log10 = col_tb_z1),
  na_col = "white", border = TRUE, annotation_name_side = "left"
)

ha_z7_axi <- HeatmapAnnotation(df = anno_z7_axi, col = list(
  Response = make_safe_pal(anno_z7_axi$Response, pal_binaria),
  Clinical_Class = make_safe_pal(anno_z7_axi$Clinical_Class, pal_granular),
  Molecular_Subgroup = make_safe_pal(anno_z7_axi$Molecular_Subgroup, pal_coo),
  Grade3_CRS = make_safe_pal(anno_z7_axi$Grade3_CRS, NULL),
  Grade3_NE = make_safe_pal(anno_z7_axi$Grade3_NE, NULL),
  EFS_Event = make_safe_pal(anno_z7_axi$EFS_Event, NULL),
  Tumor_Burden_Log10 = col_tb_z7),
  na_col = "white", border = TRUE, annotation_name_side = "left"
)

ha_z7_soc <- HeatmapAnnotation(df = anno_z7_soc, col = list(
  Response = make_safe_pal(anno_z7_soc$Response, pal_binaria),
  Molecular_Subgroup = make_safe_pal(anno_z7_soc$Molecular_Subgroup, pal_coo),
  EFS_Event = make_safe_pal(anno_z7_soc$EFS_Event, NULL),
  Tumor_Burden_Log10 = col_tb_soc),
  na_col = "white", border = TRUE, annotation_name_side = "left"
)

ha_z1_paired <- HeatmapAnnotation(df = anno_z1_paired, col = list(
  Response = make_safe_pal(anno_z1_paired$Response, pal_binaria),
  Timepoint = c("Baseline" = "gray85", "Post-CART" = "black")),
  border = TRUE, annotation_name_side = "left"
)

## =============================================================================
## 7. Renderização dos Heatmaps
## =============================================================================
zscore_and_cap <- function(mat, cap = 2) {
  mat_scaled <- t(scale(t(mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled[mat_scaled > cap] <- cap
  mat_scaled[mat_scaled < -cap] <- -cap
  return(mat_scaled)
}

build_editorial_heatmap <- function(expr_mat, signature_list, anno_df, top_ha, split_col, title, prefix) {
  genes_present <- unlist(signature_list, use.names = FALSE) %>% intersect(rownames(expr_mat))
  valid_samples <- intersect(rownames(anno_df), colnames(expr_mat))
  if(length(genes_present) == 0 | length(valid_samples) == 0) return(NULL)
  
  anno_df <- anno_df[valid_samples, , drop = FALSE]
  mat_scaled <- zscore_and_cap(expr_mat[genes_present, valid_samples, drop = FALSE])
  message("✅ Desenhando ", prefix, " | Matriz: ", nrow(mat_scaled), " genes x ", ncol(mat_scaled), " amostras.")
  
  gene_group <- rep(names(signature_list), lengths(signature_list))
  names(gene_group) <- unlist(signature_list, use.names = FALSE)
  valid_groups <- intersect(names(signature_list), unique(gene_group[genes_present]))
  row_split <- factor(gene_group[genes_present], levels = valid_groups)
  
  column_split <- if(is.null(split_col)) NULL else factor(anno_df[[split_col]])
  col_fun_heatmap <- colorRamp2(c(-2, 0, 2), c("#313695", "#FFFFBF", "#A50026"))
  
  pdf(paste0("results/exhaustion_hla/", prefix, "_Heatmap.pdf"), width = 11, height = 9)
  ht <- Heatmap(
    mat_scaled, name = "Z-score", col = col_fun_heatmap, top_annotation = top_ha, 
    row_split = row_split, column_split = column_split, row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
    cluster_rows = FALSE, cluster_column_slices = FALSE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 9), row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title = title, border = TRUE
  )
  invisible(draw(ht, merge_legends = TRUE))
  dev.off()
}

message("\n--- Gerando Heatmaps Editoriais ---")

build_editorial_heatmap(z1_base$expr, exhaustion_signature_groups, anno_z1, ha_z1, "Response", "ZUMA-1 Core Exhaustion", "ZUMA1_Exhaustion")
build_editorial_heatmap(expr_z7_axi_combined, exhaustion_signature_groups, anno_z7_axi, ha_z7_axi, "Response", "ZUMA-7 Axi-Cel Core Exhaustion", "ZUMA7_AxiCel_Exhaustion")
build_editorial_heatmap(z7_soc$expr, exhaustion_signature_groups, anno_z7_soc, ha_z7_soc, "Response", "ZUMA-7 SOC Core Exhaustion", "ZUMA7_SOC_Exhaustion")
build_editorial_heatmap(expr_z1_paired, exhaustion_signature_groups, anno_z1_paired, ha_z1_paired, "Timepoint", "ZUMA-1 Paired Core Exhaustion", "ZUMA1_Paired_Exhaustion")

build_editorial_heatmap(z1_base$expr, hla_signature_groups, anno_z1, ha_z1, "Response", "ZUMA-1 HLA", "ZUMA1_HLA")
build_editorial_heatmap(z7_mrna$expr, hla_signature_groups, anno_z7_axi, ha_z7_axi, "Response", "ZUMA-7 Axi-Cel HLA", "ZUMA7_AxiCel_HLA")
build_editorial_heatmap(z7_soc$expr, hla_signature_groups, anno_z7_soc, ha_z7_soc, "Response", "ZUMA-7 SOC HLA", "ZUMA7_SOC_HLA")
build_editorial_heatmap(expr_z1_paired, hla_signature_groups, anno_z1_paired, ha_z1_paired, "Timepoint", "ZUMA-1 Paired HLA", "ZUMA1_Paired_HLA")

message("\n✅ Script 06 Finalizado com Sucesso e Pronto para Publicação!")
