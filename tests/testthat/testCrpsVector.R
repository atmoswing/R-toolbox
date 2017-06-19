library(atmoswing)
context("CRPS vector")

test_that("CRPS vector is correct on all analogues", {
  A <- atmoswing::parseAllNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                    1, 'calibration')
  
  res <- atmoswing::crpsVector(A)
  val <- unname(res[14])
  ref <- atmoswing::crps(c(0.057982057, 0.091804922, 0.025769804, 0.000000000, 
                          0.000000000, 0.000000000, 0.000000000, 0.035433479, 
                          0.000000000, 0.000000000), 0.048318382)
  expect_equal(val, ref, tolerance = .00001)
  
})

test_that("CRPS vector is correct on 5 analogues", {
  A <- atmoswing::parseAllNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                    1, 'calibration')
  
  res <- atmoswing::crpsVector(A, 5)
  val <- unname(res[14])
  ref <- atmoswing::crps(c(0.057982057, 0.091804922, 0.025769804, 0.000000000, 
                           0.000000000), 0.048318382)
  expect_equal(val, ref, tolerance = .00001)
  
})

test_that("CRPS vector is correct on too many analogues", {
  A <- atmoswing::parseAllNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                    1, 'calibration')

  res <- atmoswing::crpsVector(A, 20)
  val <- unname(res[14])
  ref <- atmoswing::crps(c(0.057982057, 0.091804922, 0.025769804, 0.000000000, 
                           0.000000000, 0.000000000, 0.000000000, 0.035433479, 
                           0.000000000, 0.000000000), 0.048318382)
  expect_equal(val, ref, tolerance = .00001)
  
})
