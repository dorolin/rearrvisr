## ------------------------------------------------------------------------
## remove columns that contain only zeros from a matrix
## ------------------------------------------------------------------------

removeZeros<-function(mydata,mpos,markernames,ordfocal,s){
    if(length(mpos)>1){
        if(is.vector(mydata)){
            mydata<-as.matrix(mydata,ncol=1)
        }
        mydata<-mydata[,colSums(mydata)!=0]
        if(is.vector(mydata)){
            mydata<-as.matrix(mydata,ncol=1)
        }
    }else if(length(mpos)==1){
        keep<-mydata[mydata>0]
        mydata<-matrix(keep,nrow=1)
        rownames(mydata)<-markernames
    }else{
        stop(paste("Scaffold",ordfocal[s],"has no marker"))
    }
    return(mydata)
}

## ------------------------------------------------------------------------

