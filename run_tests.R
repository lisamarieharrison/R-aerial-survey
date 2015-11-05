library('RUnit')

pathnames <- list.files(pattern="[.]R$", path="C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/functions/", full.names=TRUE)
invisible(sapply(pathnames, FUN = source))

test.suite <- defineTestSuite("example",
                              dirs = file.path("C:/Users/Lisa/Documents/phd/aerial survey/R/code/R-aerial-survey/tests"),
                              testFileRegexp = 'test_distance_sampling_analysis')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)

