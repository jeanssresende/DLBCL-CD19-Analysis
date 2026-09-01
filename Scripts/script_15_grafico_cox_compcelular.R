################################################################################
# Projeto: Mecanismos de Resistência ao CAR-T anti-CD19
# Script: 06 - Figura Resumo da Regressão Multivariada
# Objetivo: Gerar um "Tile Plot" elegante com os resultados do Cox Multivariado
################################################################################

library(tidyverse)
library(ggplot2)

# 1. Carregar os dados gerados no passo anterior
cox_res <- read_csv("results/multivariate_cox/ZUMA7_Multivariate_Cox_Independence.csv")

# 2. Preparar os dados para a figura
plot_data <- cox_res %>%
  mutate(
    # Log2 do Hazard Ratio para centralizar a cor no zero (Branco)
    Log2_HR = log2(HR),
    
    # Criar os símbolos de significância clássicos
    Significance = case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01 ~ "**",
      P_Value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    
    # Limpar os nomes das populações para ficarem bonitos no eixo X
    Adjusted_By = str_replace_all(Adjusted_By, "_", " ")
  )

# 3. O Pulo do Gato: Ordenar os genes para que os "Melhores" fiquem no topo!
gene_order <- plot_data %>%
  group_by(Target_Gene) %>%
  summarise(Yes_Count = sum(Independent == "YES")) %>%
  arrange(desc(Yes_Count), Target_Gene) %>%
  pull(Target_Gene)

plot_data$Target_Gene <- factor(plot_data$Target_Gene, levels = rev(gene_order))

# 4. Criar a Figura de Alta Qualidade (Tile Plot)
p <- ggplot(plot_data, aes(x = Adjusted_By, y = Target_Gene)) +
  
  # Os quadrados coloridos com base no Hazard Ratio
  geom_tile(aes(fill = Log2_HR), color = "white", linewidth = 0.8) +
  
  # Adicionar os asteriscos (pretos para significativos, brancos/cinzas para os 'ns')
  geom_text(aes(label = Significance, color = Independent), size = 5, vjust = 0.7) +
  
  # Escala de Cores: Azul (Fator de Proteção), Branco (Neutro), Vermelho (Risco)
  scale_fill_gradient2(
    low = "#3C5488FF", mid = "white", high = "#E64B35FF", midpoint = 0,
    name = "Risco (Log2 HR)"
  ) +
  
  # Cores dos textos: Preto destaca o sucesso, Cinza esconde a falha
  scale_color_manual(values = c("YES" = "black", "NO" = "gray80"), guide = "none") +
  
  # Design Limpo (Tema minimalista)
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold.italic", color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, color = "gray30", margin = margin(b=15)),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "Independência Preditiva da Assinatura Genómica",
    subtitle = "Modelos Cox Ajustados pela Composição do Microambiente Tumoral",
    x = "População Celular Ajustada (Co-variável MCP-counter)",
    y = "Gene Alvo"
  )

# 5. Salvar a figura em formato vetorial (PDF) e PNG de alta resolução
ggsave("results/multivariate_cox/Figure_Multivariate_Summary.pdf", p, width = 11, height = 9)
ggsave("results/multivariate_cox/Figure_Multivariate_Summary.png", p, width = 11, height = 9, dpi = 300)

cat("\n==================================================================\n")
cat("Figura 'Figure_Multivariate_Summary.pdf' gerada com sucesso!\n")
cat("==================================================================\n")