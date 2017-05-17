library(atmoswing)
context("Optimal nb analogues")

test_that("Optimal nb of analogues number is correct", {
  A <- atmoswing::parseNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                 1, 'calibration')
  
  crps.nb.analogs <- atmoswing::crpsNbAnalogs(A, 5)
  
  res <- atmoswing::optimalNbAnalog(crps.nb.analogs$crps)
  
  expect_equal(unname(res$nb.min.first[1]), 1)
  expect_equal(unname(res$nb.min.med[1]), 5.5)
  expect_equal(unname(res$nb.min.last[1]), 10)
  
  expect_equal(unname(res$nb.min.first[4]), 1)
  expect_equal(unname(res$nb.min.med[4]), 1.5)
  expect_equal(unname(res$nb.min.last[4]), 2)
  
  expect_equal(unname(res$nb.min.first[5]), 1)
  expect_equal(unname(res$nb.min.med[5]), 1)
  expect_equal(unname(res$nb.min.last[5]), 1)
  
  res.smoothed <- atmoswing::optimalNbAnalog(crps.nb.analogs$crps.smoothed)
  
  expect_equal(unname(res.smoothed$nb.min.first[1]), 3)
  expect_equal(unname(res.smoothed$nb.min.med[1]), 5.5)
  expect_equal(unname(res.smoothed$nb.min.last[1]), 8)
  
  expect_equal(unname(res.smoothed$nb.min.first[4]), 8)
  expect_equal(unname(res.smoothed$nb.min.med[4]), 8)
  expect_equal(unname(res.smoothed$nb.min.last[4]), 8)
  
  expect_equal(unname(res.smoothed$nb.min.first[16]), 3)
  expect_equal(unname(res.smoothed$nb.min.med[16]), 3.5)
  expect_equal(unname(res.smoothed$nb.min.last[16]), 4)
})
