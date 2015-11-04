library('RUnit')

source('C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/distance_sampling_analysis.R')

test.suite <- defineTestSuite("example",
                              dirs = file.path("C:/Users/Lisa/Documents/phd/aerial survey/R/code/tests"),
                              testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)
