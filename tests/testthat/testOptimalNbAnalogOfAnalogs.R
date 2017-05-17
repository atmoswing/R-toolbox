library(atmoswing)
context("Optimal nb analogues of analogues")

test_that("Optimal nb of analogues of analogues is correct", {
  A <- atmoswing::parseNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                 1, 'calibration')
  
  crps.nb.analogs <- atmoswing::crpsNbAnalogs(A, 5)
  
  target.a.nb <- atmoswing::optimalNbAnalog(crps.nb.analogs$crps.smoothed)
  
  res <- atmoswing::optimalNbAnalogOfAnalogs(A, target.a.nb$nb.min.med)
  
  expect_equal(res[5, 2], 3)
  expect_equal(res[6, 5], 5.5)
  expect_equal(res[24, 10], 8)
})
