context("computeRearrs - assignOri3")
library(rearrvisr)

## assignOri3<-function(blocksL,maskL,allelem,elemrows,
##                      leafelem,leaves,testorientation)


## test sets
## -- 1 --
allelem1<-c(1,2,5,4,7,8,3,6)
dupli1<-numeric()
leafelem1<-c(1,1,1,1,1,0,1,1)
testorientation1<-c(1,1,1,1,1,1,1,1,1)
leaves1<-c(1,1,1,1,1,0,0,1,1)
elemrows1<-rbind(c(1,2,3,4,5,6,8,9),
                 c(1,2,3,4,5,7,8,9))
blocksL1<-makeBlocks(vector("list",0),1:length(allelem1),
                     1:length(allelem1),allelem1,leafelem1,1,
                     tomask=logical(0),OneIter=FALSE)
maskL1<-setMasks(blocksL1,allelem1,dupli1,tomask=logical(0))
## -- 2 --
allelem2<-c(2,1,5,4,6,7,3)
dupli2<-numeric()
leafelem2<-c(1,1,1,1,1,0,1)
testorientation2<-c(1,1,1,1,1,1,1,1)
leaves2<-c(1,1,1,1,1,0,0,1)
elemrows2<-rbind(c(1,2,3,4,5,6,8),
                 c(1,2,3,4,5,7,8))
blocksL2<-makeBlocks(vector("list",0),1:length(allelem2),
                     1:length(allelem2),allelem2,leafelem2,1,
                     tomask=logical(0),OneIter=FALSE)
maskL2<-setMasks(blocksL2,allelem2,dupli2,tomask=logical(0))
## -- 3 --
allelem3<-c(2,1)
dupli3<-numeric()
leafelem3<-c(1,1)
testorientation3<-c(1,1)
leaves3<-c(1,1)
elemrows3<-rbind(c(1,2),
                 c(1,2))
blocksL3<-makeBlocks(vector("list",0),1:length(allelem3),
                     1:length(allelem3),allelem3,leafelem3,1,
                     tomask=logical(0),OneIter=FALSE)
maskL3<-setMasks(blocksL3,allelem3,dupli3,tomask=logical(0))
## -- 4 --
allelem4<-c(2,1)
dupli4<-numeric()
leafelem4<-c(1,1)
testorientation4<-c(-1,-1)
leaves4<-c(1,1)
elemrows4<-rbind(c(1,2),
                 c(1,2))
blocksL4<-makeBlocks(vector("list",0),1:length(allelem4),
                     1:length(allelem4),allelem4,leafelem4,1,
                     tomask=logical(0),OneIter=FALSE)
maskL4<-setMasks(blocksL4,allelem4,dupli4,tomask=logical(0))
## -- 5 --
allelem5<-c(7,8,9,4,5,6,1,2,3,10,11,12,14,13)
dupli5<-numeric()
leafelem5<-c(1,1,1,1,1,0,1,0,1,1,1,1,1,1)
testorientation5<-c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
leaves5<-c(1,1,1,1,1,0,0,1,0,0,0,0,1,1,1,1,1,1)
elemrows5<-rbind(c(1,2,3,4,5,6,8,9,13,14,15,16,17,18),
                 c(1,2,3,4,5,7,8,12,13,14,15,16,17,18))
blocksL5<-makeBlocks(vector("list",0),1:length(allelem5),
                     1:length(allelem5),allelem5,leafelem5,1,
                     tomask=logical(0),OneIter=FALSE)
maskL5<-setMasks(blocksL5,allelem5,dupli5,tomask=logical(0))
## -- 6 --
allelem6<-c(3,2,1,3)
dupli6<-3
leafelem6<-c(0,1,1,0)
testorientation6<-c(1,1,1,1)
leaves6<-c(0,1,1,0)
elemrows6<-rbind(c(1,2,3,4),
                 c(1,2,3,4))
blocksL6<-makeBlocks(vector("list",0),1:length(allelem6),
                     1:length(allelem6),allelem6,leafelem6,1,
                     tomask=logical(0),OneIter=FALSE)
maskL6<-setMasks(blocksL6,allelem6,dupli6,tomask=logical(0))
## -- 7 --
allelem7<-c(2,3,2,1,3)
dupli7<-c(2,3)
leafelem7<-c(0,0,0,1,0)
testorientation7<-c(1,1,1,1,1)
leaves7<-c(0,0,0,1,0)
elemrows7<-rbind(c(1,2,3,4,5),
                 c(1,2,3,4,5))
blocksL7<-makeBlocks(vector("list",0),1:length(allelem7),
                     1:length(allelem7),allelem7,leafelem7,1,
                     tomask=logical(0),OneIter=FALSE)
maskL7<-setMasks(blocksL7,allelem7,dupli7,tomask=logical(0))
## -- 8 --
allelem8<-c(4,1,2,3,1,2,3)
dupli8<-c(1,2,3)
leafelem8<-c(1,0,0,0,0,0,0)
testorientation8<-c(1,1,1,1,1,1,1,1)
leaves8<-c(1,0,0,0,0,0,0,0)
elemrows8<-rbind(c(1,2,3,4,5,6,8),
                 c(1,2,3,4,5,7,8))
blocksL8<-makeBlocks(vector("list",0),1:length(allelem8),
                     1:length(allelem8),allelem8,leafelem8,1,
                     tomask=logical(0),OneIter=FALSE)
maskL8<-setMasks(blocksL8,allelem8,dupli8,tomask=logical(0))




## tests
test_that("function output of assignOri3 with test set 1", {
    expect_identical(assignOri3(blocksL1,maskL1,allelem1,elemrows1,
                                leafelem1,leaves1,testorientation1), 1)
})
test_that("function output of assignOri3 with test set 2", {
    expect_identical(assignOri3(blocksL2,maskL2,allelem2,elemrows2,
                                leafelem2,leaves2,testorientation2), -1)
})
test_that("function output of assignOri3 with test set 3", {
    expect_identical(assignOri3(blocksL3,maskL3,allelem3,elemrows3,
                                leafelem3,leaves3,testorientation3), 1)
})
test_that("function output of assignOri3 with test set 4", {
    expect_identical(assignOri3(blocksL4,maskL4,allelem4,elemrows4,
                                leafelem4,leaves4,testorientation4), -1)
})
test_that("function output of assignOri3 with test set 5", {
    expect_identical(assignOri3(blocksL5,maskL5,allelem5,elemrows5,
                                leafelem5,leaves5,testorientation5), 1)
})
test_that("function output of assignOri3 with test set 6", {
    expect_identical(assignOri3(blocksL6,maskL6,allelem6,elemrows6,
                                leafelem6,leaves6,testorientation6), -1)
})
test_that("function output of assignOri3 with test set 7", {
    expect_identical(assignOri3(blocksL7,maskL7,allelem7,elemrows7,
                                leafelem7,leaves7,testorientation7), -1)
})
test_that("function output of assignOri3 with test set 8", {
    expect_identical(assignOri3(blocksL8,maskL8,allelem8,elemrows8,
                                leafelem8,leaves8,testorientation8), 9)
})

