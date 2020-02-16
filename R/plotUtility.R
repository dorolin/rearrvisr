## ------------------------------------------------------------------------
## set up plot area given the number of blocks and rearrangements
## ------------------------------------------------------------------------

setRearrPlot<-function(BLOCKS,ordfocal,space,blockwidth,mar,
                       remThld=0.05,pad=0.7,y0pad=5){

    ## checks done in genomeRearrPlot

    nlev<-0
    nblocks<-0
    nrearr<-0

    mysubset<-match(ordfocal,names(BLOCKS))

    for(s in mysubset){
        nlev<-nlev+(ncol(BLOCKS[[s]]$blocks)-5)/9+1
        nblocks<-max(nblocks,nrow(BLOCKS[[s]]$blocks))
        allrearrs<-cbind(BLOCKS[[s]]$NM1,BLOCKS[[s]]$NM2,BLOCKS[[s]]$SM,
                         BLOCKS[[s]]$IV,BLOCKS[[s]]$IVsm)
        if(ncol(allrearrs)>0){
            newrearrs<-sum(apply(allrearrs,2,function(x) sum(x>remThld)>0))
        }else{
            newrearrs<-0
        }
        nrearr<-nrearr+newrearrs
    }

    myxlim<-c(-pad,nblocks*blockwidth+pad)
    ytmp<-length(mysubset)*(space$nmark+space$carid+space$stat)+
        (length(mysubset)-1)*space$scaff+
            (nlev-length(mysubset))*space$blkid+
                (nlev-length(mysubset))*space$ndori+
                    (nlev-length(mysubset))*space$blkori+
                        (nlev-length(mysubset))*space$elem+
                            nrearr*space$rearr

    myylim<-c(-(pad+y0pad),ytmp+pad)

    ## calculate height and width for the graphic
    myheight<-(myylim[2]-myylim[1])*0.2+sum(mar[c(1,3)])*0.2
    mywidth<-(myxlim[2]-myxlim[1])*0.5+sum(mar[c(2,4)])*0.2
    ## mar*0.2 converts units in lines to inches

    return(list(x=myxlim,y=myylim,height=myheight,width=mywidth))
}

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## get vector that has approximate position of breakpoints in bp
##  from precomputed start and end positions in SYNT matrix
##  (matrix has to be a subset of one scaffold);
##  approximate position is midpoint between markers provided
##  by startpos and endpos
## ------------------------------------------------------------------------

getBreakpnts2BP<-function(bptS,bptE,startpos,endpos,remThld=0){

    if(!is.matrix(bptS) | !is.matrix(bptE)){
        stop("Require matrix as input")
    }
    if((nrow(bptS) != nrow(bptE)) |
       (length(startpos) != length(endpos)) |
       (nrow(bptS) != length(startpos))){
        stop("Breakpoint matrices and marker positions cannot be matched")
    }
    brptsS<-matrix(NA,nrow=0,ncol=3)
    brptsE<-matrix(NA,nrow=0,ncol=3)
    if(ncol(bptS)>0){
        allS<-apply(bptS,1,max)
        posS<-which(allS>remThld)
        if(length(posS)){
            ## tag starting breakpoints
            for(n in 1:length(posS)){
                if(posS[n]==1){
                    ## breakpnt is at start of scaffold
                    brptsS<-rbind(brptsS,c(startpos[posS[n]],
                                           startpos[posS[n]],
                                           allS[posS[n]]))
                }else{
                    ## breakpnt in middle of scaffold
                    brptsS<-rbind(brptsS,c(endpos[posS[n]-1],
                                           startpos[posS[n]],
                                           allS[posS[n]]))
                }
            }
        }
    }
    if(ncol(bptE)>0){
        allE<-apply(bptE,1,max)
        posE<-which(allE>remThld)
        if(length(posE)){
            ## tag ending breakpoints
            for(n in 1:length(posE)){
                if(posE[n]==length(endpos)){
                    ## breakpnt is at end of scaffold
                    brptsE<-rbind(brptsE,c(endpos[posE[n]],
                                           endpos[posE[n]],
                                           allE[posE[n]]))
                }else{
                    ## breakpnt in middle of scaffold
                    brptsE<-rbind(brptsE,c(endpos[posE[n]],
                                           startpos[posE[n]+1],
                                           allE[posE[n]]))
                }
            }
        }
    }
    colnames(brptsS)<-c("start","end","valS")
    colnames(brptsE)<-c("start","end","valE")

    ## get max val for overlaps of start and end breakpoints
    tmp<-merge(brptsS,brptsE,by=c("start","end"))
    tmp<-as.matrix(cbind(tmp[,1:2],apply(tmp[,3:4],1,max)))
    ## bind all together
    if(nrow(tmp)>0){
        brptsS<-brptsS[-(unique(match(tmp[,1],brptsS[,1]),
                                match(tmp[,2],brptsS[,2]))),,drop=FALSE]
        brptsE<-brptsE[-(unique(match(tmp[,1],brptsE[,1]),
                                match(tmp[,2],brptsE[,2]))),,drop=FALSE]
    }
    brpts<-rbind(brptsS,brptsE,tmp)
    brpts<-brpts[order(brpts[,1],brpts[,2]),]

    ## get midpoint in base pairs
    brpts<-data.frame(bptmid=apply(brpts[,1:2],1,mean),
                      bptmin=brpts[,1],bptmax=brpts[,2],maxtagval=brpts[,3])

    return(brpts)
}

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## summarize tags (adjust values according to number of genes per
##  tagged event, and only keep max value per gene)
## ------------------------------------------------------------------------

summarizeTags<-function(mydata,mpos,markernames,ordfocal,s,remThld=0){
    mydata0<-mydata
    if(remThld>0){
        mydata0[mydata<=remThld]<-0
    }
    ## remove columns that are zero
    tmp<-removeZeros(mydata0,mpos,markernames,ordfocal,s)
    ## modify tags according to the number of markers per event
    if(ncol(tmp)>0){
        for(i in 1:ncol(tmp)){
            tags<-which(tmp[,i]>0)
            tmp[tags,i]<-tmp[tags,i]/length(tags)
        }
        tmp<-apply(tmp,1,max)
    }else{
        tmp<-0
    }
    return(tmp)
}

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## simplify BLOCKS tags for NM2 and SM (remove duplicate flank tags,
##  and remove flanks that are already covered by 1.0 tags)
## ------------------------------------------------------------------------

simplifyBlockTags<-function(BLOCKS,remThld=0){

    ## both flanks are tagged with 0.5 if flanks are of the same size,
    ##  or if both flanks are already covered by inserts

    for(s in 1:length(names(BLOCKS))){

        NM2blocks<-BLOCKS[[s]]$NM2
        SMblocks<-BLOCKS[[s]]$SM

        ## (1) remove duplicate flanks
        ## ---------------------------
        ##  (if not a flank it should never be duplicated)
        ## NM2
        toRm<-apply(NM2blocks,2,function(x) sum(x==0.5)==sum(x>0)) &
            duplicated(NM2blocks,MARGIN=2)
        if(sum(toRm)>0){
            NM2blocks<-NM2blocks[,!toRm,drop=FALSE]
        }

        ## SM
        toRm<-apply(SMblocks,2,function(x) sum(x==0.5)==sum(x>0)) &
            duplicated(SMblocks,MARGIN=2)
        if(sum(toRm)>0){
            SMblocks<-SMblocks[,!toRm,drop=FALSE]
        }

        ## (2a) remove 0.5 flanks already covered by 1.0 tags (i.e., inserts)
        ## ------------------------------------------------------------------
        ## NM2
        tmp<-apply(NM2blocks,2,function(x) sum(x==0.5)==sum(x>0))
        if(sum(tmp)>0 & sum(!tmp)>0){
            ## putative flanks and putative inserts exist
            flanks<-getGapsFlanks(NM2blocks[,tmp,drop=FALSE])$Flanks
            ## (can be other things than flanks)
            inserts<-NM2blocks[,!tmp,drop=FALSE]
            ## everything that is not a flank
            ##  (can be other things than inserts)
            tmpToRm<-numeric()
            for(f in 1:length(flanks)){
                myflank<-flanks[[f]]
                if(ncol(myflank)>1){
                    ## only flank if at least two parts involved
                    isCovered<-numeric(ncol(myflank))
                    for(p in 1:ncol(myflank)){
                        tmptag<-numeric(nrow(NM2blocks))
                        tmptag[(myflank[2,p]):(myflank[3,p])]<-1
                        if(duplicated(cbind(inserts,tmptag),MARGIN=2)[ncol(inserts)+1]){
                            isCovered[p]<-1
                        }
                        ## checks for inserts, which are tagged with 1.0
                        ##  (assuming remWgt > 0)
                    }
                    if(sum(isCovered)==ncol(myflank)){
                        ## all flanks are individually covered by 1.0 tags
                        tmpToRm<-c(tmpToRm,f)
                    }
                }
            }
            if(length(tmpToRm)>0){
                ## remove column from NM2
                toRm<-which(tmp==TRUE)[tmpToRm]
                NM2blocks<-NM2blocks[,-toRm,drop=FALSE]
            }
        }

        ## SM
        tmp<-apply(SMblocks,2,function(x) sum(x==0.5)==sum(x>0))
        if(sum(tmp)>0 & sum(!tmp)>0){
            ## putative flanks and putative inserts exist
            flanks<-getGapsFlanks(SMblocks[,tmp,drop=FALSE])$Flanks
            ## (can be other things than flanks)
            inserts<-SMblocks[,!tmp,drop=FALSE]
            ## everything that is not a flank
            ##  (can be other things than inserts)
            tmpToRm<-numeric()
            for(f in 1:length(flanks)){
                myflank<-flanks[[f]]
                if(ncol(myflank)>1){
                    ## only flank if at least two parts involved
                    isCovered<-numeric(ncol(myflank))
                    for(p in 1:ncol(myflank)){
                        tmptag<-numeric(nrow(SMblocks))
                        tmptag[(myflank[2,p]):(myflank[3,p])]<-1
                        if(duplicated(cbind(inserts,tmptag),MARGIN=2)[ncol(inserts)+1]){
                            isCovered[p]<-1
                        }
                        ## checks for inserts, which are tagged with 1.0
                        ##  (assuming remWgt > 0)
                    }
                    if(sum(isCovered)==ncol(myflank)){
                        ## all flanks are individually covered by 1.0 tags
                        tmpToRm<-c(tmpToRm,f)
                    }
                }
            }
            if(length(tmpToRm)>0){
                ## remove column from SM
                toRm<-which(tmp==TRUE)[tmpToRm]
                SMblocks<-SMblocks[,-toRm,drop=FALSE]
            }
        }

        ## (2b) remove (1-remThld) flanks already covered by 1.0 tags
        ## ----------------------------------------------------------
        ## NM2
        tmp<-NM2blocks
        tmp[tmp<=remThld]<-0
        tmp[tmp>=1-remThld]<-1
        toRm<-apply(NM2blocks,2,function(x) sum(x==remThld | x==(1-remThld))==sum(x>0) &
                        (sum(x==remThld)>0 & sum(x==(1-remThld))>0)) &
                            duplicated(tmp,MARGIN=2)
        ## flanks and the tagged flank also exists as insert
        if(sum(toRm)>0){
            NM2blocks<-NM2blocks[,!toRm,drop=FALSE]
        }

        ## SM
        tmp<-SMblocks
        tmp[tmp<=remThld]<-0
        tmp[tmp>=1-remThld]<-1
        toRm<-apply(SMblocks,2,function(x) sum(x==remThld | x==(1-remThld))==sum(x>0) &
                        (sum(x==remThld)>0 & sum(x==(1-remThld))>0)) &
                            duplicated(tmp,MARGIN=2)
        ## flanks and the tagged flank also exists as insert
        if(sum(toRm)>0){
            SMblocks<-SMblocks[,!toRm,drop=FALSE]
        }


        BLOCKS[[s]]$NM2<-NM2blocks
        BLOCKS[[s]]$SM<-SMblocks

    }

    return(BLOCKS)
}

## ------------------------------------------------------------------------



