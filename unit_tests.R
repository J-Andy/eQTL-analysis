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



# Test that t-tests return p-values
test_that("Test that t-tests return p-values", {
  # Sample data mimicking mutation and fold change
  sample_data <- data.frame(
    mutation_present = c(0, 0, 1, 1),
    fold_change = c(0.5, 0.6, 0.8, 0.9),
    gene = "geneA",
    mutation = "mutA"
  )
  
  # Perform Welch's t-test
  t_test_result <- t.test(fold_change ~ mutation_present, data = sample_data, var.equal = FALSE)
  
  # Check that the p-value is valid (should be numeric and less than 1)
  expect_true(is.numeric(t_test_result$p.value), info = "The p-value should be numeric")
  expect_true(t_test_result$p.value < 1, info = "The p-value should be less than 1")
})


# Print message if all tests pass
print("All tests passed successfully!")
