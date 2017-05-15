library(atmoswing)
context("CRPS nb analogues")

test_that("CRPS per analogue number is correct", {
  A <- atmoswing::parseNcOutputs(file.path('test_files', 'optimizer-outputs', '1', 'results'),
                                 1, 'calibration')
  
  res <- atmoswing::crpsNbAnalogs(A, 5)
  
  val1 <- unname(res$crps[23, 8])
  ref1 <- atmoswing::crps(c(0.2593, 0.2335, 0.2190, 0.1401, 0.2609, 0.0161, 
                            0.1611, 0.1933), 0.1675)
  expect_equal(val1, ref1, tolerance = .00001)
  
  val2 <- unname(res$crps[37, 10])
  ref2 <- atmoswing::crps(c(0.0660, 0.1111, 0.0000, 0.1482, 0.0000, 0.2964, 
                            0.1836, 0.1578, 0.1305, 0.0322), 0.1417)
  expect_equal(val1, ref1, tolerance = .00001)
  
  val1.ma <- unname(res$crps.smoothed[23, 8])
  expect_equal(val1.ma, 0.028172, tolerance = .00001)
})
