context("computeRearrs - tagTP2")
library(rearrvisr)

## tagTP2<-function(synt,allelem,tmprows,elemrows,TPelem,n,node,leaves,
##                  testorientation,preMasks,splitnodes,remWgt=0.05)


## common definitions
elemrows<-rbind(c(1,2,3,4,5,6,8,9,13,14,15),
                c(1,2,3,4,5,7,8,12,13,14,15))
tmprows<-(1:15)+2
leaves<-c(1,1,1,1,1,0,0,1,0,0,0,0,1,1,1)
testorientation<-c(1,1,1,1,-1,-1,1,1,1,-1,1,-1,-1,1,-1)
n<-3
nhier<-4
node<-"Q"
preMasks<-list(A=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                   FALSE,FALSE,FALSE),
               D=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                   FALSE,FALSE,FALSE))
splitnodes<-TRUE
synt<-list(TLWC=matrix(0,ncol=0,nrow=length(leaves)+2),
           IV=matrix(0,ncol=0,nrow=length(leaves)+2),
           TLWCbS=matrix(0,ncol=0,nrow=length(leaves)+2),
           TLWCbE=matrix(0,ncol=0,nrow=length(leaves)+2),
           IVbS=matrix(0,ncol=0,nrow=length(leaves)+2),
           IVbE=matrix(0,ncol=0,nrow=length(leaves)+2),
           nodeori=matrix(NA,ncol=nhier,nrow=length(leaves)+2),
           blockori=matrix(NA,ncol=nhier,nrow=length(leaves)+2),
           blockid=matrix(NA,ncol=nhier,nrow=length(leaves)+2),
           premask=matrix(NA,ncol=nhier,nrow=length(leaves)+2),
           subnode=matrix(0,ncol=nhier,nrow=length(leaves)+2))

## test sets
## -- 1 -- (no dupl)
allelem1<-c(10,9,5,4,3,7,6,8,1,2,11)
TPelem1<-matrix(0,ncol=11,nrow=0)
## -- 2 -- (7 dupl 2x)
allelem2<-c(1,2,5,4,3,7,6,7,10,9,11)
TPelem2<-matrix(c(0,0,0,0,0,1,2,1,0,0,0),nrow=1)
## -- 3 -- (7 dupl 3x)
allelem3<-c(1,2,4,7,3,7,6,5,7,9,11)
TPelem3<-rbind(c(0,0,0,1,2,1,0,0,0,0,0),c(0,0,0,0,0,1,2,2,1,0,0))
## -- 4 -- (4 dupl 3x)
allelem4<-c(4,2,1,4,3,4,7,6,5,8,9)
TPelem4<-rbind(c(1,2,2,1,0,0,0,0,0,0,0),c(0,0,0,1,2,1,0,0,0,0,0))
## -- 5 -- (multiple dupl)
allelem5<-c(1,2,3,4,3,2,4,5,6,7,8)
TPelem5<-rbind(c(0,1,2,2,2,1,0,0,0,0,0),c(0,0,1,2,1,0,0,0,0,0,0),
               c(0,0,0,1,2,2,1,0,0,0,0))
## -- 6 -- (multiple dupl)
allelem6<-c(2,3,2,3,4,5,8,7,6,3,9)
TPelem6<-rbind(c(1,2,1,0,0,0,0,0,0,0,0),c(0,1,2,1,0,0,0,0,0,0,0),
               c(0,0,0,1,2,2,2,2,2,1,0))


## test results
## -- 1 --
TLWC1<-rbind(c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.05, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.05, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.95, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.05, 0.00),
             c(0.05, 0.00, 0.00, 0.00, 0.00, 0.95),
             c(0.00, 0.95, 0.00, 0.00, 0.00, 0.00))

IV1<-rbind(c(0,   0,   0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0,   0,   0),
           c(0,   0,   1,   0,   0,   0,   0),
           c(0,   0,   0,   1,   0,   0,   0),
           c(1,   1,   0,   0,   1,   0,   0),
           c(1,   1,   0,   0,   0,   1,   0),
           c(1,   1,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0,   0,   1),
           c(0,   0,   0,   0,   0,   0,   0))

## -- 2 --
TLWC2<-rbind(c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00),
             c(0.00, 0.95, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.95, 0.00),
             c(0.00, 0.00, 0.00, 0.05),
             c(0.00, 0.00, 0.00, 0.00))

IV2<-rbind(c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(1,   1,   0,   0,   0),
           c(1,   0,   1,   0,   0),
           c(1,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   1,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   1))

## -- 3 --
TLWC3<-rbind(c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.05, 0.00, 0),
             c(0.00, 0.00, 1),
             c(0.00, 0.95, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0),
             c(0.00, 0.00, 0))

IV3<-rbind(c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   0),
           c(0,   1,   0,   0,   0),
           c(1,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0),
           c(1,   0,   1,   0,   0),
           c(1,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0),
           c(1,   0,   0,   0,   0),
           c(0,   0,   0,   1,   0),
           c(0,   0,   0,   0,   0),
           c(0,   0,   0,   0,   1))

## -- 4 --
TLWC4<-rbind(c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.95, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.05, 0.5, 0.0, 0.00, 0.00),
             c(0.00, 0.05, 0.0, 0.5, 0.00, 0.00),
             c(0.00, 0.05, 0.0, 0.0, 0.05, 0.00),
             c(0.00, 0.05, 0.0, 0.0, 0.00, 0.95),
             c(0.00, 0.05, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.05, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00),
             c(0.00, 0.00, 0.0, 0.0, 0.00, 0.00))

IV4<-rbind(c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   1,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(1,   0,   1,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   1))

## -- 5 --
TLWC5<-rbind(c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.95, 0.00, 0.00, 0.00),
             c(0.95, 0.00, 0.00, 0.00),
             c(0.00, 0.05, 0.95, 0.00),
             c(0.00, 0.05, 0.00, 0.05),
             c(0.00, 0.05, 0.00, 0.05),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00))

IV5<-rbind(c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(1,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   0,   0),
           c(0,   1,   0),
           c(0,   0,   0),
           c(0,   0,   1))

## -- 6 --
TLWC6<-rbind(c(0.00, 0.00, 0.00, 0.00),
             c(0.00, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00),
             c(0.05, 0.00, 0.00, 0.00),
             c(0.00, 0.95, 0.00, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.05, 0.00),
             c(0.00, 0.00, 0.00, 0.95),
             c(0.00, 0.00, 0.00, 0.00))

IV6<-rbind(c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   1,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   0),
           c(1,   0,   1,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(1,   0,   0,   0),
           c(0,   0,   0,   0),
           c(0,   0,   0,   1))


## tests
test_that("function output of tagTP2 with test set 1 (no dupl)", {
    expect_identical(tagTP2(synt,allelem1,tmprows,elemrows,TPelem1,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC1)
    expect_identical(tagTP2(synt,allelem1,tmprows,elemrows,TPelem1,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV1)
})
test_that("function output of tagTP2 with test set 2 (7 dupl 2x)", {
    expect_identical(tagTP2(synt,allelem2,tmprows,elemrows,TPelem2,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC2)
    expect_identical(tagTP2(synt,allelem2,tmprows,elemrows,TPelem2,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV2)
})
test_that("function output of tagTP2 with test set 3 (7 dupl 3x)", {
    expect_identical(tagTP2(synt,allelem3,tmprows,elemrows,TPelem3,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC3)
    expect_identical(tagTP2(synt,allelem3,tmprows,elemrows,TPelem3,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV3)
})
test_that("function output of tagTP2 with test set 4 (4 dupl 3x)", {
    expect_identical(tagTP2(synt,allelem4,tmprows,elemrows,TPelem4,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC4)
    expect_identical(tagTP2(synt,allelem4,tmprows,elemrows,TPelem4,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV4)
})
test_that("function output of tagTP2 with test set 5 (multiple dupl)", {
    expect_identical(tagTP2(synt,allelem5,tmprows,elemrows,TPelem5,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC5)
    expect_identical(tagTP2(synt,allelem5,tmprows,elemrows,TPelem5,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV5)
})
test_that("function output of tagTP2 with test set 6 (multiple dupl)", {
    expect_identical(tagTP2(synt,allelem6,tmprows,elemrows,TPelem6,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$TLWC, TLWC6)
    expect_identical(tagTP2(synt,allelem6,tmprows,elemrows,TPelem6,n,node,
                            leaves,testorientation,preMasks,splitnodes,
                            remWgt=0.05)$IV, IV6)
})



