library('RUnit')
 
source('morrislib.R')
 
test.suite <- defineTestSuite("genome",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')
 
test.result <- runTestSuite(test.suite)
 
printTextProtocol(test.result)

