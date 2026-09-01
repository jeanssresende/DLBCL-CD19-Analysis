# 1. Carregar os pacotes necessários
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(ggpubr)

# ==============================================================================
# PASSO 1: BAIXAR DADOS DE EXPRESSÃO (RNA-Seq) DO TCGA-DLBC
# ==============================================================================
cat("Baixando dados de expressão do TCGA-DLBC...\n")
query_exp <- GDCquery(
  project = "TCGA-DLBC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query_exp)
exp_data <- GDCprepare(query_exp)

# Extrair a matriz de expressão (vamos usar TPM para normalização)
matriz_tpm <- assay(exp_data, "tpm_unstrand")

# Encontrar o ID do gene WNT11
info_genes <- rowData(exp_data)
id_wnt11 <- rownames(info_genes)[info_genes$gene_name == "WNT11"]

# Criar um dataframe com a expressão do WNT11 por paciente
df_exp <- data.frame(
  Barcode_Completo = colnames(matriz_tpm),
  WNT11_TPM = as.numeric(matriz_tpm[id_wnt11, ])
) %>%
  # O TCGA tem barcodes longos. Vamos pegar só os primeiros 12 caracteres (ID do paciente)
  mutate(Patient_ID = substr(Barcode_Completo, 1, 12))

# ==============================================================================
# PASSO 2: BAIXAR DADOS DE MUTAÇÃO (SNV) DO TCGA-DLBC
# ==============================================================================
cat("Baixando dados de mutação...\n")
query_mut <- GDCquery(
  project = "TCGA-DLBC",
  data.category = "Simple Nucleotide Variation",
  access = "open",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query_mut)
mut_data <- GDCprepare(query_mut)

# Filtrar pacientes que têm mutação no gene CD79B
pacientes_mutados_cd79b <- mut_data %>% 
  filter(Hugo_Symbol == "CD79B") %>% 
  pull(Tumor_Sample_Barcode) %>% 
  substr(1, 12) %>% 
  unique()

# ==============================================================================
# PASSO 3: JUNTAR OS DADOS E CLASSIFICAR
# ==============================================================================
# Vamos cruzar os pacientes que têm RNA-seq com a lista de mutados
df_final <- df_exp %>%
  # Remove duplicatas caso haja múltiplas amostras do mesmo paciente
  distinct(Patient_ID, .keep_all = TRUE) %>%
  mutate(
    CD79B_Status = ifelse(Patient_ID %in% pacientes_mutados_cd79b, "Mutated", "Wild-Type")
  )

# Mostrar a contagem de quantos são mutados e quantos são Wild-Type
table(df_final$CD79B_Status)

# ==============================================================================
# PASSO 4: ESTATÍSTICA E GRÁFICO (COM ESCALA LOGARÍTMICA)
# ==============================================================================

# Criar uma nova coluna com a transformação log2(TPM + 1)
df_final <- df_final %>%
  mutate(WNT11_log2 = log2(WNT11_TPM + 1))

# Teste estatístico (permanece o mesmo, mas podemos rodar com a nova coluna)
teste_estatistico <- wilcox.test(WNT11_log2 ~ CD79B_Status, data = df_final)
print(teste_estatistico)

# Criar o Boxplot com a escala correta e visual limpo
grafico_log <- ggboxplot(
  df_final, 
  x = "CD79B_Status", 
  y = "WNT11_log2",
  color = "CD79B_Status", 
  palette = c("#E7B800", "#00AFBB"),
  add = "jitter", # adiciona os pontinhos dos pacientes
  ylab = "WNT11 Expression (log2 TPM + 1)", 
  xlab = "CD79B Mutation Status",
  title = "WNT11 Expression in TCGA-DLBC"
) + 
  stat_compare_means(method = "wilcox.test", label.x.npc = "center") +
  theme_minimal() +
  theme(legend.position = "none")

# Visualizar o gráfico corrigido
print(grafico_log)

# Salvar em alta resolução
ggsave("WNT11_vs_CD79B_TCGA_LogScale.pdf", plot = grafico_log, width = 6, height = 5, dpi = 300)