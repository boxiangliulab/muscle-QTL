# Load required libraries
library(ggplot2)
library(dplyr)
library(stringr)

# Set the working directory
setwd("/Users/wenjingwang/NUS/Liulab/muscleQTL/differential_expression/Plot/")

# Load pathway enrichment data
load_pathway_data <- function(filename, description, type, denominator = 516) {
  read.delim(filename) %>%
    mutate(
      direction = type,
      des = paste0(Description, ' (', description, ')'),
      perDE = Count / denominator * 100,
      realp = qvalue
    )
}

# Load data for different types of analyses
down_GO <- load_pathway_data("go_down.txt", "GO", "Down")
down_REA <- load_pathway_data("rea_down.txt", "Reactome", "Down")
down_KEGG <- load_pathway_data("kegg_down.txt", "KEGG", "Down")
up_GO <- load_pathway_data("go_up.txt", "GO", "Up")
up_KEGG <- load_pathway_data("kegg_up.txt", "KEGG", "Up")
up_REA <- load_pathway_data("rea_up.txt", "Reactome", "Up")
up_GSEA <- load_pathway_data("gsea_up.txt", "GSEA", "Up", setSize)

# Combine all data and filter by significance
all_pathways <- bind_rows(down_GO, down_REA, down_KEGG, up_GO, up_KEGG, up_REA, up_GSEA) %>%
  filter(realp < 0.06)

# Renaming geneID column across different data sources
all_pathways <- all_pathways %>%
  rename(geneID = starts_with("geneID")) %>%
  arrange(desc(realp)) 

# Select specific pathways for the plot
selected_indices <- c(6, 8, 9, 25, 33, 34, 37, 40, 44, 55, 68, 70, 80, 102, 157, 177, 205, 299, 301, 304)
plot_data <- all_pathways[selected_indices, ]

# Assign classes to pathways
plot_data$class <- factor(
  c("Metabolism", "Skeletal Muscle", "Skeletal Muscle", "Metabolism", "Others", "Others",
    "Others", "Others", "Insulin-related Pathway", "Insulin-related Pathway", "Metabolism",
    "Diet", "Skeletal Muscle", "Others", "Insulin-related Pathway", "Insulin-related Pathway",
    "Metabolism", "Diet", "Insulin-related Pathway", "Others"),
  levels = c("Insulin-related Pathway", "Metabolism", "Skeletal Muscle", "Diet", "Energy", "Others")
)

# Save the data for future use
write.csv(plot_data, file = "enrichment_for_plot.txt", quote = TRUE, col.names = TRUE)

# Adjusting the text and finalizing the plot
enrichment_for_plot <- plot_data %>%
  mutate(des = str_to_title(des))

# Create the plot
P1 <- ggplot(enrichment_for_plot, aes(x = neg_log10_pvalue, y = des)) +
  geom_point(aes(size = perDE, shape = Significant, color = Significant)) +
  scale_shape_manual(values = c(15, 19)) +
  scale_color_manual(values = c('blue', 'red')) +
  facet_grid(class ~ ., scales = "free", space = "free") +
  theme_bw() +
  labs(x = NULL, y = NULL, shape = "", size = "% DE genes") +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(P1)

# Save the plot to a file
ggsave('Enrichment_Pathway_v3.pdf', plot = P1, width = 17, height = 11)

