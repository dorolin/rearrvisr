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
        allrearrs<-cbind(BLOCKS[[s]]$TLBS,BLOCKS[[s]]$TLWS,BLOCKS[[s]]$TLWC,
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
        if(remThld>0){
            allS[allS<=remThld]<-0
            allS[allS>=1-remThld]<-1
        }
        posS<-which(allS>0)
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
        if(remThld>0){
            allE[allE<=remThld]<-0
            allE[allE>=1-remThld]<-1
        }
        posE<-which(allE>0)
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
    brpts<-data.frame(bpt=apply(brpts[,1:2],1,mean),val=brpts[,3])

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
## simplify BLOCKS tags for TLWS and TLWC (remove duplicate flank tags,
##  and remove flanks that are already covered by 1.0 tags)
## ------------------------------------------------------------------------

simplifyBlockTags<-function(BLOCKS,remThld=0){

    ## both flanks are tagged with 0.5 if flanks are of the same size,
    ##  or if both flanks are already covered by inserts

    for(s in 1:length(names(BLOCKS))){

        TLWSblocks<-BLOCKS[[s]]$TLWS
        TLWCblocks<-BLOCKS[[s]]$TLWC

        ## (1) remove duplicate flanks
        ## ---------------------------
        ##  (if not a flank it should never be duplicated)
        ## TLWS
        toRm<-apply(TLWSblocks,2,function(x) sum(x==0.5)==sum(x>0)) &
            duplicated(TLWSblocks,MARGIN=2)
        if(sum(toRm)>0){
            TLWSblocks<-TLWSblocks[,!toRm,drop=FALSE]
        }

        ## TLWC
        toRm<-apply(TLWCblocks,2,function(x) sum(x==0.5)==sum(x>0)) &
            duplicated(TLWCblocks,MARGIN=2)
        if(sum(toRm)>0){
            TLWCblocks<-TLWCblocks[,!toRm,drop=FALSE]
        }

        ## (2a) remove 0.5 flanks already covered by 1.0 tags (i.e., inserts)
        ## ------------------------------------------------------------------
        ## TLWS
        tmp<-apply(TLWSblocks,2,function(x) sum(x==0.5)==sum(x>0))
        if(sum(tmp)>0 & sum(!tmp)>0){
            ## putative flanks and putative inserts exist
            flanks<-getGapsFlanks(TLWSblocks[,tmp,drop=FALSE])$Flanks
            ## (can be other things than flanks)
            inserts<-TLWSblocks[,!tmp,drop=FALSE]
            ## everything that is not a flank
            ##  (can be other things than inserts)
            tmpToRm<-numeric()
            for(f in 1:length(flanks)){
                myflank<-flanks[[f]]
                if(ncol(myflank)>1){
                    ## only flank if at least two parts involved
                    isCovered<-numeric(ncol(myflank))
                    for(p in 1:ncol(myflank)){
                        tmptag<-numeric(nrow(TLWSblocks))
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
                ## remove column from TLWS
                toRm<-which(tmp==TRUE)[tmpToRm]
                TLWSblocks<-TLWSblocks[,-toRm,drop=FALSE]
            }
        }

        ## TLWC
        tmp<-apply(TLWCblocks,2,function(x) sum(x==0.5)==sum(x>0))
        if(sum(tmp)>0 & sum(!tmp)>0){
            ## putative flanks and putative inserts exist
            flanks<-getGapsFlanks(TLWCblocks[,tmp,drop=FALSE])$Flanks
            ## (can be other things than flanks)
            inserts<-TLWCblocks[,!tmp,drop=FALSE]
            ## everything that is not a flank
            ##  (can be other things than inserts)
            tmpToRm<-numeric()
            for(f in 1:length(flanks)){
                myflank<-flanks[[f]]
                if(ncol(myflank)>1){
                    ## only flank if at least two parts involved
                    isCovered<-numeric(ncol(myflank))
                    for(p in 1:ncol(myflank)){
                        tmptag<-numeric(nrow(TLWCblocks))
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
                ## remove column from TLWC
                toRm<-which(tmp==TRUE)[tmpToRm]
                TLWCblocks<-TLWCblocks[,-toRm,drop=FALSE]
            }
        }

        ## (2b) remove (1-remThld) flanks already covered by 1.0 tags
        ## ----------------------------------------------------------
        ## TLWS
        tmp<-TLWSblocks
        tmp[tmp<=remThld]<-0
        tmp[tmp>=1-remThld]<-1
        toRm<-apply(TLWSblocks,2,function(x) sum(x==remThld | x==(1-remThld))==sum(x>0) &
                        (sum(x==remThld)>0 & sum(x==(1-remThld))>0)) &
                            duplicated(tmp,MARGIN=2)
        ## flanks and the tagged flank also exists as insert
        if(sum(toRm)>0){
            TLWSblocks<-TLWSblocks[,!toRm,drop=FALSE]
        }

        ## TLWC
        tmp<-TLWCblocks
        tmp[tmp<=remThld]<-0
        tmp[tmp>=1-remThld]<-1
        toRm<-apply(TLWCblocks,2,function(x) sum(x==remThld | x==(1-remThld))==sum(x>0) &
                        (sum(x==remThld)>0 & sum(x==(1-remThld))>0)) &
                            duplicated(tmp,MARGIN=2)
        ## flanks and the tagged flank also exists as insert
        if(sum(toRm)>0){
            TLWCblocks<-TLWCblocks[,!toRm,drop=FALSE]
        }


        BLOCKS[[s]]$TLWS<-TLWSblocks
        BLOCKS[[s]]$TLWC<-TLWCblocks

    }

    return(BLOCKS)
}

## ------------------------------------------------------------------------



