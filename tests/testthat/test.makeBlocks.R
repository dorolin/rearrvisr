context("computeRearrs - makeBlocks")
library(rearrvisr)

## makeBlocks<-function(blocksL,blstart,blend,blelem,leafelem,
##                      iter,tomask=logical(0),OneIter=FALSE)


## test sets
allelem<-c(1,2,3,4,3,2,4,5,6)
leafelem<-c(1,1,1,1,1,1,1,1,1)
blockori<-c(1,1,-1,-1,1,1,1,1,1)

## test results
blocksL.out<-list(rbind(c(1,   2,   1,   2,   2,   1,   1),
                        c(3,   3,   3,   3,   1,   9,   3),
                        c(4,   4,   4,   4,   1,   9,   4),
                        c(5,   5,   3,   3,   1,   9,   3),
                        c(6,   6,   2,   2,   1,   9,   2),
                        c(7,   9,   4,   6,   3,   1,   5)),
                  rbind(c(1,   2,   1,   1,   2,   9,   1),
                        c(3,   4,   3,   4,   2,   1,   3),
                        c(5,   6,   3,   2,   2,  -1,   2),
                        c(7,   9,   5,   5,   3,   9,   4)),
                  rbind(c(1,   2,   1,   1,   2,   9,   1),
                        c(3,   6,   3,   2,   4,  -1,   2),
                        c(7,   9,   4,   4,   3,   9,   3)),
                  rbind(c(1,   9,   1,   3,   9,   1,   1)))


## tests
test_that("function output of makeBlocks with test set 1", {
    expect_identical(makeBlocks(vector("list",0),1:length(allelem),
                                1:length(allelem),allelem,leafelem,1,
                                tomask=logical(0),OneIter=FALSE), blocksL.out)
})


