library(atmoswing)
context("Parsing outputs")

test_that("NetCDF outputs are correctly parsed for 1 level of analogy, calibration", {
  expect_true(is.dir(file.path('files', 'optimizer-outputs', '1', 'results')))
  A <- atmoswing::parseNcOutputs(file.path('files', 'optimizer-outputs', '1', 'results'),
                                 1, 'calibration')
  
  expect_equal(A$analog.dates.MJD[1,1], 53777.0)
  expect_equal(A$analog.dates.MJD[19,1], 48960.0)
  expect_equal(A$analog.dates.MJD[23,5], 47875.0)
  
  expect_equal(A$analog.criteria[1,1], 52.8974, tolerance = .0001)
  expect_equal(A$analog.criteria[21,1], 50.9112, tolerance = .0001)
  expect_equal(A$analog.criteria[16,7], 55.3957, tolerance = .0001)
  
  expect_equal(A$analog.values.norm[1,1], 0, tolerance = .0001)
  expect_equal(A$analog.values.norm[4,3], 0.0805, tolerance = .0001)
  expect_equal(A$analog.values.norm[14,1], 0.0580, tolerance = .0001)
  expect_equal(A$analog.values.norm[4,10], 0.0725, tolerance = .0001)

  expect_equal(A$analog.values.raw[1,1], 0, tolerance = .0001)
  expect_equal(A$analog.values.raw[23,1], 16.1000, tolerance = .0001)
  expect_equal(A$analog.values.raw[19,9], 13.6000, tolerance = .0001)
  
  expect_equal(A$target.dates.MJD[6], 54837.0)
  expect_equal(A$target.dates.UTC[6], as.Date('2009-01-06'))
  
  expect_equal(A$target.values.norm[19], 0.3286, tolerance = .0001)
  
  expect_equal(A$target.values.raw[19], 20.4, tolerance = .0001)
  
  expect_equal(A$predict.score[7], 0.9594, tolerance = .0001)
  
})

test_that("NetCDF outputs are correctly parsed for 1 level of analogy, validation", {
  expect_true(is.dir(file.path('files', 'optimizer-outputs', '1', 'results')))
  A <- atmoswing::parseNcOutputs(file.path('files', 'optimizer-outputs', '1', 'results'),
                                 1, 'validation')
  
  expect_equal(A$analog.dates.MJD[1,1], 49666.0)
  expect_equal(A$analog.dates.MJD[19,1], 54466.0)
  expect_equal(A$analog.dates.MJD[23,5], 53406.0)
  
  expect_equal(A$analog.criteria[1,1], 53.6841, tolerance = .0001)
  expect_equal(A$analog.criteria[21,1], 44.9280, tolerance = .0001)
  expect_equal(A$analog.criteria[16,7], 46.1833, tolerance = .0001)
  
  expect_equal(A$analog.values.norm[1,1], 0.0177, tolerance = .0001)
  expect_equal(A$analog.values.norm[4,3], 0.0000, tolerance = .0001)
  expect_equal(A$analog.values.norm[14,1], 0.1369, tolerance = .0001)
  expect_equal(A$analog.values.norm[4,10], 0.0966, tolerance = .0001)
  
  expect_equal(A$analog.values.raw[1,1], 1.1, tolerance = .0001)
  expect_equal(A$analog.values.raw[23,1], 0.0, tolerance = .0001)
  expect_equal(A$analog.values.raw[20,8], 3.4, tolerance = .0001)
  
  expect_equal(A$target.dates.MJD[6], 55202.0)
  expect_equal(A$target.dates.UTC[6], as.Date('2010-01-06'))
  
  expect_equal(A$target.values.norm[8], 0.0435, tolerance = .0001)
  
  expect_equal(A$target.values.raw[8], 2.7, tolerance = .0001)
  
  expect_equal(A$predict.score[7], 0.4730, tolerance = .0001)
  
})

test_that("NetCDF outputs are correctly parsed for the 2nd level of analogy, calibration", {
  expect_true(is.dir(file.path('files', 'optimizer-outputs', '2', 'results')))
  A <- atmoswing::parseNcOutputs(file.path('files', 'optimizer-outputs', '2', 'results'),
                                 1, 'calibration', 2)
  
  expect_equal(A$analog.dates.MJD[1,1], 46076.0)
  expect_equal(A$analog.dates.MJD[19,1], 45320.0)
  expect_equal(A$analog.dates.MJD[23,5], 54139.0)
  
  expect_equal(A$analog.criteria[1,1], 9.7133, tolerance = .0001)
  expect_equal(A$analog.criteria[21,1], 46.9275, tolerance = .0001)
  expect_equal(A$analog.criteria[17,8], 140.9442, tolerance = .0001)
  
  expect_equal(A$analog.values.norm[1,1], 0, tolerance = .0001)
  expect_equal(A$analog.values.norm[2,3], 0.0564, tolerance = .0001)
  expect_equal(A$analog.values.norm[14,1], 0.0435, tolerance = .0001)
  expect_equal(A$analog.values.norm[13,10], 0.0403, tolerance = .0001)
  
  expect_equal(A$analog.values.raw[1,1], 0, tolerance = .0001)
  expect_equal(A$analog.values.raw[23,1], 4.3, tolerance = .0001)
  expect_equal(A$analog.values.raw[19,9], 24.6, tolerance = .0001)
  
  expect_equal(A$target.dates.MJD[6], 54837.0)
  expect_equal(A$target.dates.UTC[6], as.Date('2009-01-06'))
  
  expect_equal(A$target.values.norm[19], 0.3286, tolerance = .0001)
  
  expect_equal(A$target.values.raw[19], 20.4, tolerance = .0001)
  
  expect_equal(A$predict.score[7], 0.9725, tolerance = .0001)
  
})

test_that("NetCDF outputs are correctly parsed for the 2nd level of analogy, validation", {
  expect_true(is.dir(file.path('files', 'optimizer-outputs', '2', 'results')))
  A <- atmoswing::parseNcOutputs(file.path('files', 'optimizer-outputs', '2', 'results'),
                                 1, 'validation', 2)
  
  expect_equal(A$analog.dates.MJD[1,1], 51216.0)
  expect_equal(A$analog.dates.MJD[19,1], 50821.0)
  expect_equal(A$analog.dates.MJD[23,5], 49703.0)
  
  expect_equal(A$analog.criteria[1,1], 112.7534, tolerance = .0001)
  expect_equal(A$analog.criteria[21,1], 9.0414, tolerance = .0001)
  expect_equal(A$analog.criteria[16,7], 131.3492, tolerance = .0001)
  
  expect_equal(A$analog.values.norm[1,1], 0.0532, tolerance = .0001)
  expect_equal(A$analog.values.norm[4,3], 0.0290, tolerance = .0001)
  expect_equal(A$analog.values.norm[16,1], 0.0564, tolerance = .0001)
  
  expect_equal(A$analog.values.raw[1,1], 3.3, tolerance = .0001)
  expect_equal(A$analog.values.raw[23,1], 0.0, tolerance = .0001)
  expect_equal(A$analog.values.raw[20,8], 2.2, tolerance = .0001)
  
  expect_equal(A$target.dates.MJD[6], 55202.0)
  expect_equal(A$target.dates.UTC[6], as.Date('2010-01-06'))
  
  expect_equal(A$target.values.norm[8], 0.0435, tolerance = .0001)
  
  expect_equal(A$target.values.raw[8], 2.7, tolerance = .0001)
  
  expect_equal(A$predict.score[7], 0.7440, tolerance = .0001)
  
})
