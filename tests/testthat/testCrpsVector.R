library(atmoswing)
context("CRPS vector")

test_that("CRPS vector is correct", {
  A <- atmoswing::parseNcOutputs(file.path('test_files', 'optim', '1', 'results'),
                                 1, 'calibration')
  
  # On all analogues
  res1 <- atmoswing::crpsVector(A)
  val1 <- unname(res1[14])
  ref1 <- atmoswing::crps(c(0.057982057, 0.091804922, 0.025769804, 0.000000000, 
                            0.000000000, 0.000000000, 0.000000000, 0.035433479, 
                            0.000000000, 0.000000000), 0.048318382)
  expect_equal(val1, ref1, tolerance = .00001)
  
  # On 5 analogues
  res1 <- atmoswing::crpsVector(A, 5)
  val2 <- unname(res1[14])
  ref2 <- atmoswing::crps(c(0.057982057, 0.091804922, 0.025769804, 0.000000000, 
                            0.000000000), 0.048318382)
  expect_equal(val2, ref2, tolerance = .00001)
  
})
