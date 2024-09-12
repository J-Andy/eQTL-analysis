library(testthat)

test_that("Test merging works correctly", {
  merged <- inner_join(mutations_long, gene_kos_long, by = "Cell_Line")
  expect_true(nrow(merged) > 0)
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