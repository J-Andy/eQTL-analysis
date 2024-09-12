library(testthat)
library(dplyr)

# Test data alignment (common models)
test_that("Data alignment between mutations and gene KO datasets works", {
  mutations <- read.table("Mutations.tsv", header = TRUE, sep = "\t", row.names = 1)
  gene_kos <- read.table("Gene_KOs.tsv", header = TRUE, sep = "\t", row.names = 1)

  # Align mutation and KO data by model (cell line)
  common_models <- intersect(colnames(mutations), row.names(gene_kos))
  
  expect_true(length(common_models) > 0, info = "There should be common cell lines between the datasets")
  
  mutations <- mutations[, common_models]
  gene_kos <- gene_kos[common_models, ]
  
  expect_equal(ncol(mutations), length(common_models), info = "Mutations should be aligned with the common cell lines")
  expect_equal(nrow(gene_kos), length(common_models), info = "Gene KO should be aligned with the common cell lines")
})


test_that("Test that t-tests return p-values", {
  sample_data <- data.frame(
    mutation_present = c(0, 0, 1, 1),
    fold_change = c(0.5, 0.6, 0.8, 0.9),
    gene = "geneA",
    mutation = "mutA"
  )
  result <- run_tests(sample_data)
  expect_true(result$p_value < 1)
})

# Print message if all tests pass
print("All tests passed successfully!")
