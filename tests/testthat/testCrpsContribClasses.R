library(atmoswing)
context("CRPS contributions")

test_that("CRPS contributions are processed correctly", {
  A <- atmoswing::parseAllNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                    1, 'validation')

  crps.vect <- atmoswing::crpsVector(A, 10)
  res <- atmoswing::crpsContribClasses(crps.vect, A$target.values.raw)
  
  expect_equal(length(res$contrib), length(res$labels))
  expect_equal(sum(res$contrib), 100)
  
})
