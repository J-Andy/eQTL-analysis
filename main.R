# Load necessary libraries
library(dplyr)
library(ggplot2)
library(qvalue) 
library(ggrepel)

# Set the output directory (inside the container, this will be bind-mounted to the host)
output_dir <- "/usr/src/app/output"

# Load the mutation and gene KO datasets
mutations <- read.table("Mutations.tsv", header = TRUE, sep = "\t", row.names = 1)
gene_kos <- read.table("Gene_KOs.tsv", header = TRUE, sep = "\t", row.names = 1)

# Align mutation and KO data by model (cell line)
common_models <- intersect(colnames(mutations), row.names(gene_kos))
mutations <- mutations[, common_models]
gene_kos <- gene_kos[common_models, ]

# Initialize an empty results data frame
results <- data.frame()

# Loop over each mutation and gene pair
for (mutation in rownames(mutations)) {
  for (gene in colnames(gene_kos)) {
    
    # Extract mutation status and gene KO values
    mutation_status <- as.numeric(mutations[mutation, ])
    ko_values <- as.numeric(gene_kos[, gene])
    
    # Ensure there are both mutated (1) and non-mutated (0) cell lines
    if (length(unique(mutation_status)) == 2) {
      
      # Perform Welch's t-test
      test <- t.test(ko_values ~ mutation_status, var.equal = FALSE)
      
      # Collect results (mutation, gene, p-value, t-statistic)
      result <- data.frame(
        Mutation = mutation,
        Gene = gene,
        T_statistic = test$statistic,
        P_value = test$p.value
      )
      
      # Append to results data frame
      results <- rbind(results, result)
    }
  }
}

# FDR Correction using Benjamini-Hochberg procedure
results$FDR <- p.adjust(results$P_value, method = "BH")

# Filter for significant results (FDR < 0.05)
significant_results <- results %>% filter(FDR < 0.05)

# Save the significant results to a file in the output directory
write.table(significant_results, file.path(output_dir, "significant_mutation_gene_associations.tsv"), sep = "\t", row.names = FALSE)

# Print completion message
print("Statistical analysis complete. Significant associations saved to 'significant_mutation_gene_associations.tsv'.")


# Q-Q Plot for p-values
qqplot_pvalue <- ggplot(results, aes(sample = -log10(P_value))) +
  stat_qq(distribution = stats::qunif) + 
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(
    title = "Q-Q Plot of P-values",
    x = "Expected -log10(P)",
    y = "Observed -log10(P)"
  ) +
  theme_minimal()

# Save the Q-Q plot in the output directory
ggsave(file.path(output_dir, "qqplot_pvalues.png"), plot = qqplot_pvalue, width = 6, height = 6)


# Histogram of p-values 
hist_pvalue <- ggplot(results, aes(x = P_value)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black") +
  labs(
    title = "Histogram of P-values",
    x = "P-value",
    y = "Frequency"
  ) +
  theme_minimal()

# Save the histogram in the output directory
ggsave(file.path(output_dir, "histogram_pvalues.png"), plot = hist_pvalue, width = 6, height = 6)


# Volcano Plot
volcano_plot <- ggplot(results, aes(x = T_statistic, y = -log10(FDR))) +
  geom_point(aes(color = FDR < 0.05), size = 2) +  # Color points based on FDR significance
  scale_color_manual(values = c("grey", "red")) +  # Grey for non-significant, red for significant
  labs(
    title = "Volcano Plot",
    x = "T-statistic (Effect Size)",
    y = "-log10(FDR)"
  ) +
  theme_minimal() +
  geom_text_repel(data = filter(results, FDR < 0.05),  # Label only significant points
                  aes(label = paste(Mutation, Gene, sep = ":")),
                  size = 3)

# Save the volcano plot in the output directory
ggsave(file.path(output_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)


# For each significant mutation-gene pair, create a boxplot of KO values by mutation status
for (i in 1:nrow(significant_results)) {
  
  # Extract mutation and gene from the significant results
  mutation <- significant_results$Mutation[i]
  gene <- significant_results$Gene[i]
  
  # Get mutation status and KO values for the significant pair
  mutation_status <- as.factor(unlist(mutations[mutation, , drop = TRUE]))  
  ko_values <- unlist(gene_kos[, gene, drop = TRUE])
  
  # Combine into a data frame for plotting
  plot_data <- data.frame(
    Mutation_Status = mutation_status,
    KO_Value = ko_values
  )
  
  # Create boxplot
  boxplot <- ggplot(plot_data, aes(x = Mutation_Status, y = KO_Value)) +
    geom_boxplot(fill = "lightblue", color = "black") +
    geom_jitter(width = 0.2, alpha = 0.7, color = "darkblue") +
    labs(
      title = paste("KO Outcome by Mutation Status\n", mutation, "-", gene),
      x = "Mutation Status",
      y = "KO Outcome"
    ) +
    theme_minimal()
  
  # Save each boxplot as a separate file in the output directory
  ggsave(file.path(output_dir, paste0("boxplot_", mutation, "_", gene, ".png")), plot = boxplot, width = 6, height = 4)
}

print("Completed!")
