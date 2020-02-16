## ------------------------------------------------------------------------
## This file contains various helper functions for 'computeRearrs.R'
## ------------------------------------------------------------------------


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## NOTES

## term 'synteny' in this script used in sense of conserved order
## term 'duplicated' for TP elements refers to identical node IDs
##  occurring at different positions, not actual a gene duplication

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## potential further improvements:
##   - decision inversions versus transpositions
##   - add 'best hits' for flank/insert decision for P-nodes
##   - with duplicated elements: could take possibility of
##     inversions into account for decision which part has been
##     rearranged
##   - with duplicated elements: follow up on the idea of defining
##     "premasks" before making blocks [!!! having changed to use
##     'checkAdjComplexRearr()' function in 'tagTP2()' is likely to
##     be incompatible with using premasks without adding
##     appropriate adjustements]
##   - 'assignOri3()' and 'checkAdjAscend()' and 'checkAdjDescend()'
##     could be simplified when using new 'checkAdjComplexRearr()'
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


## ------------------------------------------------------------------------


## check whether elements in subtree are leaves or nodes
tagLeaf<-function(subtree,tmprows,n){
    if(ncol(subtree)<=(3+(n-1)*2)){ ## all are leaves
        tag<-rep(1,length(tmprows))
    }else{
        tag<-rep(0,length(tmprows))
        tag[is.na(subtree[tmprows,3+n*2])]<-1
    }
    return(tag)
}

## ------------------------------------------------------------------------


## rank function that handles ties so that ties get identical
##  ranks with no gaps to higher ranks
myRank<-function(x){
    rawrank<-rank(x,ties.method="min")
    ## adjust that no ranks are missing due to ties
    tmp<-sort(setdiff(1:length(rawrank),rawrank))
    y<-rep(0,length(rawrank))
    if(length(tmp)>0){
        for(i in 1:length(tmp)){
            y[which(rawrank>tmp[i])]<-y[which(rawrank>tmp[i])]+1
        }
    }
    newrank<-rawrank-y
    return(newrank)
}

## ------------------------------------------------------------------------


## check tree for nodes containing only one marker
## tree has to be in original order
getSingles<-function(tree){

    singlemarkers<-rep(0,nrow(tree))

    ## compare each row to preceding and subsequent row
    ##  (exceptions for first and last row)
    for(i in 1:nrow(tree)){
        pos<-which(!is.na(tree[i,-c(1:2)]))
        pos<-head(pos+2,n=-2L) ## +2 to adjust for -c(1:2) above
        curr<-paste0(tree[i,pos],collapse="-")
        if(i==1){
            if(curr!=paste0(tree[i+1,pos],collapse="-")){
                singlemarkers[i]<-1
            }
        }else if(i==nrow(tree)){
            if(curr!=paste0(tree[i-1,pos],collapse="-")){
                singlemarkers[i]<-1
            }
        }else{
            if(curr!=paste0(tree[i-1,pos],collapse="-") &
               curr!=paste0(tree[i+1,pos],collapse="-")){
                singlemarkers[i]<-1
            }
        }
    }
    return(singlemarkers)
}

## ------------------------------------------------------------------------


## additional processing transpositions (TPs)
tagTP2<-function(synt,allelem,tmprows,elemrows,TPelem,n,node,leaves,
                 testorientation,preMasks,splitnodes,remWgt=0.05,
                 testlim=100){

    ## >>> don't use preMasks, this will require new adjustements
    ##     after adding 'checkAdjComplexRearr()' function

    ## synt has dimension of subtree
    ## leaves, testorientation, tmprows have dimension of markers
    ## allelem, elemrows, TPelem have dimension of elements
    ## baseID contains name of higher-level node

    ## check for orientation of leaf markers
    ##  (only possible when markers were doubled)
    ## >>> markers that are single for a (Q-)node in the original,
    ##     unreduced tree have unresolved orientation to my understanding;
    ##     those are now set to NA;
    ## >>> orientation is not relevant for P-nodes to my
    ##     understanding. Only consider orientation for Q-nodes

    ## get right data dimensions (count number of leaf markers per element)
    leafelem<-rep(0,ncol(elemrows))
    for(i in 1:ncol(elemrows)){
        leafelem[i]<-sum(leaves[(elemrows[1,i]):(elemrows[2,i])])
    }

    if(node=="P"){
        ## for P-nodes:
        ##  either flanking elements or the inserted element
        ##  might have been the source of transposition

        subtagsmark<-integer(length(tmprows))
        elemids<-integer(length(tmprows))

        if(nrow(TPelem)>0){
            ## tagging of P-nodes
            ##      roughly comparable to what's done with CARs
            ##      during post-processing step
            ## filter out large inserts and tag adjacent flanks with 0.5
            ##  instead, otherwise, tag inserts with 1.0
            flankelem<-unique(allelem[apply(TPelem,2,
                                            function(x) sum(x==1)>0)])
            TPflanks<-matrix(0,ncol=0,nrow=length(tmprows))
            for(k in 1:length(flankelem)){
                ## bind flanks of same element together
                flankpos<-which(allelem==flankelem[k])
                tmpflank<-rep(0,length(tmprows))
                for(l in 1:length(flankpos)){
                    pos1<-elemrows[1,flankpos[l]]
                    pos2<-elemrows[2,flankpos[l]]
                    tmpflank[pos1:pos2]<-rep(1,length(pos1:pos2))
                    elemids[pos1:pos2]<-rep(allelem[flankpos[l]],
                                            length(pos1:pos2))
                }
                TPflanks<-cbind(TPflanks,tmpflank)
            }
            ## determine flank and gap sizes (use function originally
            ##  written for postprocessing at CAR level)
            GapFlank<-getGapsFlanks(TPflanks)
            ## (returns $Gaps and $Flanks)
            ## filter out large inserts and tag adjacent flanks instead,
            ##  otherwise, tag inserts only
            TPtags<-makeTPtags(TPflanks,GapFlank$Gaps,GapFlank$Flanks,
                               nelem=length(tmprows),remWgt=remWgt)
            ## returns $TPfilt and $subtags
            ## (matrix with tags for inserts/gaps and subnode IDs)

            ## append columns to TLWC
            tp<-matrix(0,nrow=nrow(synt$TLWC),ncol=ncol(TPtags$TPfilt))
            tp[tmprows,]<-TPtags$TPfilt
            synt$TLWC<-cbind(synt$TLWC,tp)

            ## store subnode IDs for markers that have received a TP tag
            subtagsmark<-TPtags$subtags

            ## obtain breakpoints for TPtags$TPfilt
            bpt<-getBreakpntsSE(TPtags$TPfilt)
            ## returns $bptS and $bptE
            ## append columns to TLWCbS/TLWCbE
            bptS<-matrix(0,nrow=nrow(synt$TLWCbS),ncol=ncol(bpt$bptS))
            bptS[tmprows,]<-bpt$bptS
            synt$TLWCbS<-cbind(synt$TLWCbS,bptS)
            bptE<-matrix(0,nrow=nrow(synt$TLWCbE),ncol=ncol(bpt$bptE))
            bptE[tmprows,]<-bpt$bptE
            synt$TLWCbE<-cbind(synt$TLWCbE,bptE)
        } ## close nrow(TPelem)>0
        synt$subnode[tmprows,n]<-subtagsmark
        synt$blockid[tmprows,n]<-elemids

    }else if(node=="Q"){
        ## for Q-nodes:
        ##  check for correct order of elements and potential inversions;
        ##  don't want to tag TPelems for Q-nodes,
        ##  as this is mostly redundant with checking for correct
        ##  ordering (dependent on which checks are performed);
        ##  in addition, TPelems might arise from inversions that
        ##  cut through an element, so this should be handled differently;
        ##  however, in some cases duplicated elements might be hidden
        ##  in otherwise perfect blocks/block-adjencies (possibly only
        ##  if more than one element duplicated)

        ## identify blocks of elements that are in order
        rankelem<-myRank(allelem)

        ## store splitblockid (if TP in 'checkOriAscend'/'checkOriDescend')
        splitblockid<-rep("",length(allelem))

        ## get vector of duplicated elements (TPelem=1)
        dupli<-numeric()
        ## in addition, bind flanks of same TP element together
        TPflanks<-matrix(0,ncol=0,nrow=length(allelem))

        if(nrow(TPelem)>0){
            dupli<-unique(allelem[apply(TPelem,2,function(x) is.element(1,x))])
            for(k in 1:length(dupli)){
                ## bind flanks of same element together
                tmpflank<-numeric(length(allelem))
                tmpflank[which(allelem==dupli[k])]<-1
                TPflanks<-cbind(TPflanks,tmpflank)
            }
        }


        ## >>>> removed the distinction between having/not having duplis
        ##      after adding 'checkAdjComplexRearr()' function
        ## if(length(dupli)>0){
        ##     ## in the presence of duplications, determine node order
        ##     ##  based on pre-determined masked / non-masked elements
        ## }else{
        ##     ## no duplications: order determination using function
        ##     ## (functions nevertheless account for duplicated elements
        ##     ##   as inherited from code history)
        ## }


        blocksL<-vector("list",0)
        blocksL<-makeBlocks(blocksL,1:length(allelem),1:length(allelem),
                            allelem,leafelem,1)
        maskL<-setMasks(blocksL,allelem,dupli)

        ## store blockid and blockori
        ##  (blockid might be further modified below)
        ## expand block over elements to markers
        for(k in 1:nrow(blocksL[[1]])){
            ## get columns in elemrows
            tmp<-(blocksL[[1]][k,1]):(blocksL[[1]][k,2])
            ## get positions of markers
            tmp2<-numeric()
            for(j in tmp){
                tmp2<-c(tmp2,(elemrows[1,j]):(elemrows[2,j]))
            }
            synt$blockori[tmprows[tmp2],n]<-rep(blocksL[[1]][k,6],
                                                length(tmp2))
            synt$blockid[tmprows[tmp2],n]<-rep(blocksL[[1]][k,7],
                                               length(tmp2))
        }

        ## storage for block element orientation
        ##  (can be overwritten in functions below)
        blockoriL<-vector("list",length(blocksL))
        for(z in 1:length(blocksL)){
            blockoriL[[z]]<-blocksL[[z]][,6]
        }

        ## check for rearrangements:

        ## - given final node orientation, check whether in last level
        ##   of blocks ([[length(blocksL)]] or [[length(blocksL)-1]])
        ##   adjacencies are correct, and whether orientation of blocks
        ##   is correct

        ## - if no final block combination exists,
        ##   - use 'keepBlocks()' function to exclude elements that
        ##     hinder clustering to final block
        ##   - use 'checkAdjComplexRearr()' function to cluster
        ##     remaining blocks, make rearrangement tags for them,
        ##     and to tag excluded elements
        ## - if final block combination exists,
        ##   - tag wrong adjacencies as TPs (adjacencies should always
        ##     be correct if final combination exists, except when
        ##     orientation was switched)
        ## - for wrong orientation, tag as IV or TP
        ## - if wrong orientation includes dupli, tag temporarily
        ##   (in extra matrix?) to potentially solve which part is
        ##   TP/IV when considering lower-level nodes; if one of
        ##   the dupli's is not bound to other elements but causes
        ##   a problem, tag it
        ##   >>> not done yet (in order to keep it simple);
        ##       also, decision on TP/IV is not perfect and
        ##       could be improved, but might become cumbersome
        ##       given the hierarchy of P- and Q-nodes
        ## - for blocks that were not further combined at higher level
        ##   (thus block orientation is 9), pass orientation from next
        ##   higher level

        ## - step through the levels of blocks and repeat



        ## take into account that there might be just a single
        ##  block of one element (nrow(blocksL[[1]])==1 & ori==9)
        ##  or a single block of several elements
        ##  (nrow(blocksL[[1]])==1 & ori!=9); in both cases,
        ##  all elements will be in order
        if(nrow(blocksL[[1]])==1){
            ori<-assignOri3(blocksL,maskL,allelem,elemrows,
                            leafelem,leaves,testorientation)
            tmptp<-matrix(NA,ncol=0,nrow=length(allelem))
            tmpiv<-matrix(NA,ncol=0,nrow=length(allelem))
            blocklevel<-1
        }else if(nrow(blocksL[[length(blocksL)]])>1){
            ## no final block; exclude smallest blocks
            ##   until final clustering possible
            tokeep<-keepBlocks(blocksL,maskL,elemrows,
                               allelem,leafelem,testlim=testlim)

            ## check adjacencies for block clustering done
            ##   with non-excluded blocks only, and make tags
            ##   for incorrect adjacencies; also make tags
            ##   for excluded blocks; ori is identified from
            ##   clustering of non-excluded blocks
            xxx<-checkAdjComplexRearr(blocksL,maskL,allelem,elemrows,
                                      leafelem,leaves,testorientation,
                                      tokeep,remWgt)
            ori<-xxx$ori
            tmptp<-xxx$tmptp
            tmpiv<-xxx$tmpiv
            invelem<-xxx$invelem
            splitblockid<-xxx$splitblockid
            expectedOri<-xxx$expectedOri
            blocklevel<-length(blocksL)
        }else{
            ori<-assignOri3(blocksL,maskL,allelem,elemrows,
                            leafelem,leaves,testorientation)
            ## final block exists, which can be skipped
            blocklevel<-length(blocksL)-1
            if(ori==1){
                tmptp<-checkAdjAscend(blocksL[[blocklevel]],
                                      maskL[[blocklevel]],
                                      length(allelem))
            }else if(ori== -1){
                tmptp<-checkAdjDescend(blocksL[[blocklevel]],
                                       maskL[[blocklevel]],
                                       length(allelem))
            }else{
                stop("Node has unexpected orientation")
            }
            ## if all block elements received a tp tag, adjust so that
            ##  largest block won't get a tag
            if(sum(rowSums(tmptp)>=1)==nrow(tmptp)){
                tmptp<-adjustTPtags(tmptp,blocksL[[blocklevel]],
                                    maskL[[blocklevel]],
                                    blocksL[[1]],elemrows,remWgt)
            }
            tmpiv<-matrix(NA,ncol=0,nrow=length(allelem))
        }
        synt$nodeori[tmprows,n]<-ori
        synt$premask[tmprows,n]<-0


        ## === check orientation of blocks (top level) ===
        if(nrow(blocksL[[length(blocksL)]])>1){
            ## continue with special treatment of blocksL[[length(blocksL)]]
            ##   after pre-processing with 'checkAdjComplexRearr()' above

            ## run checkOriAscend/checkOriDescend separately for each
            ##   element, depending on orientation in expectedOri
            for(k in 1:(nrow(blocksL[[blocklevel]]))){
                if(expectedOri[k]==9){
                    ## pass orientation from higher level
                    expectedOri[k]<-ori
                }
                ## get the ones with opposite of expected orientation
                if(expectedOri[k]==1){
                    xxx<-checkOriAscend(blocksL,blocklevel,allelem,dupli,
                                        elemrows,leaves,testorientation,
                                        blockoriL,subbl=k,splitblockid,
                                        remWgt=remWgt)
                    ## returns tpElA, ivElA, (modified) blockoriL,
                    ##  splitblockid
                    tmptp<-cbind(tmptp,xxx$tpElA)
                    tmpiv<-cbind(tmpiv,xxx$ivElA)
                    blockoriL<-xxx$blockoriL
                    splitblockid<-xxx$splitblockid
                    ## keep track of expected orientation of leaves
                    if(ncol(xxx$ivElA)>0){
                        for(i in 1:ncol(xxx$ivElA)){
                            newinv<-which(xxx$ivElA[,i]==1 & invelem==0)
                            backinv<-which(xxx$ivElA[,i]==1 & invelem==1)
                            if(length(newinv)>0){
                                invelem[newinv]<-rep(1,length(newinv))
                            }
                            if(length(backinv)>0){
                                invelem[backinv]<-rep(0,length(backinv))
                            }
                        }
                    }
                }else if(expectedOri[k]== -1){
                    xxx<-checkOriDescend(blocksL,blocklevel,allelem,dupli,
                                         elemrows,leaves,testorientation,
                                         blockoriL,subbl=k,splitblockid,
                                         remWgt=remWgt)
                    ## returns tpElD, ivElD, (modified) blockoriL,
                    ##  splitblockid
                    tmptp<-cbind(tmptp,xxx$tpElD)
                    tmpiv<-cbind(tmpiv,xxx$ivElD)
                    blockoriL<-xxx$blockoriL
                    splitblockid<-xxx$splitblockid
                    ## keep track of expected orientation of leaves
                    if(ncol(xxx$ivElD)>0){
                        for(i in 1:ncol(xxx$ivElD)){
                            newinv<-which(xxx$ivElD[,i]==1 & invelem==0)
                            backinv<-which(xxx$ivElD[,i]==1 & invelem==1)
                            if(length(newinv)>0){
                                invelem[newinv]<-rep(1,length(newinv))
                            }
                            if(length(backinv)>0){
                                invelem[backinv]<-rep(0,length(backinv))
                            }
                        }
                    }
                }
                ## for proper function with tests for remaining levels
                ##   below, ensure that blockoriL[[blocklevel]][k]!=9
                if(blockoriL[[blocklevel]][k]==9){
                    ## pass orientation from higher level
                    blockoriL[[blocklevel]][k]<-expectedOri[k]
                }
            } ## close loop over rows k

        }else{ ## nrow(blocksL[[length(blocksL)]])==1
            if(ori==1){
                xxx<-checkOriAscend(blocksL,blocklevel,allelem,dupli,
                                    elemrows,leaves,testorientation,blockoriL,
                                    subbl=1:nrow(blocksL[[blocklevel]]),
                                    splitblockid,remWgt=remWgt)
                ## returns tpElA, ivElA, (modified) blockoriL, splitblockid
                tmptp<-cbind(tmptp,xxx$tpElA)
                tmpiv<-cbind(tmpiv,xxx$ivElA)
                blockoriL<-xxx$blockoriL
                splitblockid<-xxx$splitblockid
                ## keep track of expected orientation of leaves
                invelem<-rep(0,length(allelem))
                if(ncol(xxx$ivElA)>0){
                    for(i in 1:ncol(xxx$ivElA)){
                        newinv<-which(xxx$ivElA[,i]==1 & invelem==0)
                        backinv<-which(xxx$ivElA[,i]==1 & invelem==1)
                        if(length(newinv)>0){
                            invelem[newinv]<-rep(1,length(newinv))
                        }
                        if(length(backinv)>0){
                            invelem[backinv]<-rep(0,length(backinv))
                        }
                    }
                }
            }else if(ori == -1){
                xxx<-checkOriDescend(blocksL,blocklevel,allelem,dupli,
                                     elemrows,leaves,testorientation,blockoriL,
                                     subbl=1:nrow(blocksL[[blocklevel]]),
                                     splitblockid,remWgt=remWgt)
                ## returns tpElD, ivElD, (modified) blockoriL,splitblockid
                tmptp<-cbind(tmptp,xxx$tpElD)
                tmpiv<-cbind(tmpiv,xxx$ivElD)
                blockoriL<-xxx$blockoriL
                splitblockid<-xxx$splitblockid
                ## keep track of expected orientation of leaves
                invelem<-rep(1,length(allelem))
                if(ncol(xxx$ivElD)>0){
                    for(i in 1:ncol(xxx$ivElD)){
                        newinv<-which(xxx$ivElD[,i]==1 & invelem==0)
                        backinv<-which(xxx$ivElD[,i]==1 & invelem==1)
                        if(length(newinv)>0){
                            invelem[newinv]<-rep(1,length(newinv))
                        }
                        if(length(backinv)>0){
                            invelem[backinv]<-rep(0,length(backinv))
                        }
                    }
                }
            }else{ ## ori==9
                invelem<-rep(NA,length(allelem))
            }
        }
        ## ===

        blocklevel<-blocklevel-1

        ## === check orientation of blocks (remaining levels) ===
        ## loop through remaining levels, separately for each higher-level
        ##  block (adjacencies will always be correct by definition)
        if(blocklevel>0){
            for(z in blocklevel:1){
                for(k in 1:(nrow(blocksL[[z+1]]))){
                    ## get elements to consider
                    tmp<-which(blocksL[[z]][,1]>=blocksL[[z+1]][k,1] &
                                   blocksL[[z]][,2]<=blocksL[[z+1]][k,2])
                    passedOri<-9
                    if(blockoriL[[z+1]][k]==9){
                        ## pass orientation from higher level
                        y<-z+2
                        while(y<=length(blocksL) & passedOri==9){
                            if(nrow(blocksL[[y]])==1){
                                ## avoid taking original orientation
                                ##  but take assigned ori instead
                                break
                            }
                            tmp2<-which(blocksL[[y]][,1]<=blocksL[[z+1]][k,1] &
                                            blocksL[[y]][,2]>=blocksL[[z+1]][k,2])
                            passedOri<-blockoriL[[y]][tmp2]
                            y<-y+1
                        }
                        ## if still 9, use ori (if ori would be 9 too then
                        ##  this loop would never have been entered)
                        if(passedOri==9){
                            if(nrow(blocksL[[length(blocksL)]])>1){
                                stop("Unexpected block orientation")
                                ## it should be ensured above that this
                                ##   never happens
                            }
                            passedOri<-ori
                        }
                        blockoriL[[z+1]][k]<-passedOri
                    }
                    ## get the ones with opposite of expected orientation
                    if(blockoriL[[z+1]][k]==1){
                        xxx<-checkOriAscend(blocksL,z,allelem,dupli,
                                            elemrows,leaves,testorientation,
                                            blockoriL,subbl=tmp,splitblockid,
                                            remWgt=remWgt)
                        ## returns tpElA, ivElA, (modified) blockoriL,
                        ##  splitblockid
                        tmptp<-cbind(tmptp,xxx$tpElA)
                        tmpiv<-cbind(tmpiv,xxx$ivElA)
                        blockoriL<-xxx$blockoriL
                        splitblockid<-xxx$splitblockid
                        ## keep track of expected orientation of leaves
                        if(ncol(xxx$ivElA)>0){
                            for(i in 1:ncol(xxx$ivElA)){
                                newinv<-which(xxx$ivElA[,i]==1 & invelem==0)
                                backinv<-which(xxx$ivElA[,i]==1 & invelem==1)
                                if(length(newinv)>0){
                                    invelem[newinv]<-rep(1,length(newinv))
                                }
                                if(length(backinv)>0){
                                    invelem[backinv]<-rep(0,length(backinv))
                                }
                            }
                        }
                    }else if(blockoriL[[z+1]][k]== -1){
                        xxx<-checkOriDescend(blocksL,z,allelem,dupli,
                                             elemrows,leaves,testorientation,
                                             blockoriL,subbl=tmp,splitblockid,
                                             remWgt=remWgt)
                        ## returns tpElD, ivElD, (modified) blockoriL,
                        ##  splitblockid
                        tmptp<-cbind(tmptp,xxx$tpElD)
                        tmpiv<-cbind(tmpiv,xxx$ivElD)
                        blockoriL<-xxx$blockoriL
                        splitblockid<-xxx$splitblockid
                        ## keep track of expected orientation of leaves
                        if(ncol(xxx$ivElD)>0){
                            for(i in 1:ncol(xxx$ivElD)){
                                newinv<-which(xxx$ivElD[,i]==1 & invelem==0)
                                backinv<-which(xxx$ivElD[,i]==1 & invelem==1)
                                if(length(newinv)>0){
                                    invelem[newinv]<-rep(1,length(newinv))
                                }
                                if(length(backinv)>0){
                                    invelem[backinv]<-rep(0,length(backinv))
                                }
                            }
                        }
                    }
                } ## close loop over rows
            } ## close loop over blocklevels
        } ## close blocklevel>0


        ## make subnode IDs for elements that have received a TP tag
        ##  >>>> could be interesting to store this for IVs too, but
        ##       separately, to check whether IV including a dupli
        ##       would be recovered at next lower level
        ##  note that here, unlike for flanks of cars or
        ##   P-nodes, only the tagged piece of a two-element-
        ##   translocation, or the tagged parts when segment of
        ##   largest block has been removed, gets a subnode ID
        subtags<-integer(length(allelem))
        if(ncol(tmptp)>0){
            tagmat<-tmptp
            ## exclude removed parts of two-element-translocation
            ##  or when segment of largest block has been untagged
            if(remWgt>0){
                tagmat[tagmat<=remWgt]<-0
                tagmat[tagmat>=0.5]<-1
            }
            ## check if some columns have gaps in tags, as each flank
            ##  will need separate ID
            ## >>>> with new 'checkAdjComplexRearr()' function,
            ##      gaps in tags are indeed possible, but I don't think
            ##      that they should get separate tags as they are part
            ##      of the same rearrangement (otherwise, there shouldn't
            ##      be any gaps and computing multags was already removed)
            ## multags<-apply(tagmat,2,function(x)
            ##     length(which(diff(which(x>0))>1))>0)
            ## if(sum(multags)>0){
            ##     cat("Warning: did not expect to find gaps in tags\n")
            ##     ## >>>> don't think that it will possible for tmptp,
            ##     ##      even if yes, I don't think that flanks should
            ##     ##      get separate tags (?)
            ##     ## for(i in which(multags*1==1)){
            ##     ##     tmp<-which(tagmat[,i]>0)
            ##     ##     tmp2<-which(diff(tmp)>1)
            ##     ##     tagstart<-c(1,tmp2+1)
            ##     ##     tagend<-c(tmp2,length(tmp))
            ##     ##     for(j in 1:length(tagstart)){
            ##     ##         tagmat[tmp[(tagstart[j]):(tagend[j])],i]<-j
            ##     ##     }
            ##     ## }
            ## }
            tmptags<-apply(tagmat,1,function(x) paste0(x,collapse="-"))
            unitags<-unique(tmptags)
            ## exclude markers that are untagged
            notag<-paste0(rep(0,ncol(tagmat)),collapse="-")
            unitags<-unitags[unitags!=notag]
            for(i in 1:length(unitags)){
                tagpos<-which(tmptags==unitags[i])
                subtags[tagpos]<-rep(i,length(tagpos))
            }
        }


        ## expand all tags to markers
        ## get positions of markers
        tpmark<-matrix(NA,ncol=ncol(tmptp),nrow=length(tmprows))
        ivmark<-matrix(NA,ncol=ncol(tmpiv),nrow=length(tmprows))
        invmark<-rep(NA,length(tmprows))
        subtagsmark<-rep(NA,length(tmprows))
        splitblockmark<-rep(NA,length(tmprows))
        for(j in 1:length(allelem)){
            tmp2<-(elemrows[1,j]):(elemrows[2,j])
            if(ncol(tpmark)>0){
                for(i in tmp2){
                    tpmark[i,]<-tmptp[j,]
                }
            }
            if(ncol(ivmark)>0){
                for(i in tmp2){
                    ivmark[i,]<-tmpiv[j,]
                }
            }
            for(i in tmp2){
                invmark[i]<-invelem[j]
                subtagsmark[i]<-subtags[j]
                splitblockmark[i]<-splitblockid[j]
            }
        }


        ## for leaf markers and markers doubled, check that orientation
        ##  of markers within lowest-level block is correct
        if(splitnodes==TRUE & length(tmprows)<2){
            ## (added test for length(tmprows)>1 after adding
            ##   to work with subnode ids, as this can result
            ##   in single-marker subnodes that don't have a
            ##   'NA' for testorientation in check below)
            syntoritagsIV<-rep(0,length(tmprows))
        }else if(sum(!is.na(testorientation))>0){
            ## check that tags for inversions are -1, and outside 1,
            ##  but only for leafmarkers
            syntoritagsIV<-rep(0,length(tmprows))

            ## tag remaining conflicts as single-marker inversions
            syntoritagsIV[invmark==1 & !is.na(testorientation) &
                            leaves==1 & testorientation==1]<-1
            syntoritagsIV[invmark==0 & !is.na(testorientation) &
                            leaves==1 & testorientation== -1]<-1
        }else{
            syntoritagsIV<-rep(0,length(tmprows))
        }
        ## append columns to ivmark (these are single-marker
        ##  inversions, so each needs its own column)
        if(sum(syntoritagsIV>0)){
            for(i in which(syntoritagsIV==1)){
                tmptags<-rep(0,length(tmprows))
                tmptags[i]<-1
                ivmark<-cbind(ivmark,tmptags)
            }
        }

        ## get positions in tree and bind tags to synt,
        ##  and obtain breakpoints
        if(ncol(tpmark)>0){
            tptree<-matrix(0,nrow=nrow(synt$TLWC),ncol=ncol(tpmark))
            tptree[tmprows,]<-tpmark
            synt$TLWC<-cbind(synt$TLWC,tptree)
            ## obtain breakpoints for tpmark
            bpt<-getBreakpntsSE(tpmark)
            ## returns $bptS and $bptE
            bptS<-matrix(0,nrow=nrow(synt$TLWCbS),ncol=ncol(bpt$bptS))
            bptS[tmprows,]<-bpt$bptS
            synt$TLWCbS<-cbind(synt$TLWCbS,bptS)
            bptE<-matrix(0,nrow=nrow(synt$TLWCbE),ncol=ncol(bpt$bptE))
            bptE[tmprows,]<-bpt$bptE
            synt$TLWCbE<-cbind(synt$TLWCbE,bptE)
        }
        if(ncol(ivmark)>0){
            ivtree<-matrix(0,nrow=nrow(synt$IV),ncol=ncol(ivmark))
            ivtree[tmprows,]<-ivmark
            synt$IV<-cbind(synt$IV,ivtree)
            ## obtain breakpoints for ivmark
            bpt<-getBreakpntsSE(ivmark)
            ## returns $bptS and $bptE
            bptS<-matrix(0,nrow=nrow(synt$IVbS),ncol=ncol(bpt$bptS))
            bptS[tmprows,]<-bpt$bptS
            synt$IVbS<-cbind(synt$IVbS,bptS)
            bptE<-matrix(0,nrow=nrow(synt$IVbE),ncol=ncol(bpt$bptE))
            bptE[tmprows,]<-bpt$bptE
            synt$IVbE<-cbind(synt$IVbE,bptE)
        }
        synt$subnode[tmprows,n]<-subtagsmark
        ## modify tags for blockid
        synt$blockid[tmprows,n]<-paste0(synt$blockid[tmprows,n],splitblockmark)

    } ## close node=="Q"
    return(synt)

}

## ------------------------------------------------------------------------


## check if CAR - scaffold assignments indicate translocations
tagTLcar2<-function(markers,tree,myscafs,mycars,chromLev=0){
    ## - for each CAR, only certain number/combination of boundaries to
    ##   other CARs are allowed
    ## - for all sets of CARs, if boundaries are spread on different
    ##   scaffolds, these are not allowed to result in conflicts
    ## - when considering not scaffolds but full chromosomes
    ##   the number and positions of allowed boundaries becomes
    ##   more restrictive (>>>> not tested)

    ## >>>> TODO: clean up a bunch of unused if-statements

    TLbetween<-matrix(NA,ncol=0,nrow=nrow(markers))
    rownames(TLbetween)<-tree$marker
    TLwithin<-TLbetween

    scafcarviolB<-matrix(0,nrow=length(myscafs),ncol=length(mycars))
    rownames(scafcarviolB)<-myscafs
    colnames(scafcarviolB)<-mycars
    scafcarviolW<-scafcarviolB

    ## require at least two CARs
    if(length(mycars)<2){
        return(list(TLbetween=TLbetween,TLwithin=TLwithin,
                    scafcarviolB=scafcarviolB))
    }

    ## determine size of scaffolds
    scafsize<-integer(length(myscafs))
    for(i in 1:length(myscafs)){
        scafsize[i]<-sum(markers$scaff==myscafs[i])
    }

    ## count boundaries
    ## (i) internal boundaries (not aligning with scaffold ends)
    scafcar<-matrix(0,nrow=length(myscafs),ncol=length(mycars))
    rownames(scafcar)<-myscafs
    colnames(scafcar)<-mycars
    ## (ii) end boundaries (aligning with scaffold ends)
    endbnd<-scafcar

    for(s in 1:length(myscafs)){
        scaffpos<-which(markers$scaff==myscafs[s])
        tmpcars<-sort(unique(tree$car[scaffpos]))
        if(length(tmpcars)>1){ ## at least two cars -> min one boundary
            for(i in tmpcars){
                carpos<-which(tree$car[scaffpos]==i)
                bound<-2*sum(diff(carpos)!=1) ## gaps count two boundaries
                bound<-bound+(carpos[1]!=1) ## CAR starts in middle of scaffold
                bound<-bound+(tail(carpos,n=1L)!=length(scaffpos))
                ## CAR starts at end of scaffold
                endb<-1*(carpos[1]==1)
                endb<-endb+(tail(carpos,n=1L)==length(scaffpos))
                ## should also work if CAR or scaffold has only one marker
                scafcar[s,mycars==i]<-scafcar[s,mycars==i]+bound
                endbnd[s,mycars==i]<-endbnd[s,mycars==i]+endb
            }
        }else if(length(tmpcars)==1){ ## CARs that fully span a scaffold
            endbnd[s,mycars==tmpcars]<-endbnd[s,mycars==tmpcars]+2
        }
    }

    ## tag CARs that have too many boundaries
    for(i in 1:length(mycars)){
        if(sum(scafcar[,i] + endbnd[,i]) > 2){
            if(sum(scafcar[,i] + endbnd[,i] > 0)>1){
                ## multiple scaffolds involved -> might be TLcar
                if(chromLev==1){
                    ## chromosome-level assembly of focal species
                    ## -> always TLcar
                    scafcarviolB[(scafcar[,i] + endbnd[,i]) > 0,i]<-1
                }else{
                    if(sum(scafcar[,i])==2 & sum(scafcar[,i]>1)==0){
                        ## two internal boundaries, each on different scaffold
                        ## might not be a TLcar, but needs further testing (*,$)
                        ## -> done below

                    }else if(sum(scafcar[,i])>=2 & sum(scafcar[,i]>0)==1){
                        ## >= 2 internal boundaries on single scaffold
                        ## might not be a TLcar, but needs further testing (#)
                        intscaf<-which(scafcar[,i]>=2)
                        if(endbnd[intscaf,i]==0){
                            fullscaf<-which(endbnd[,i]==2)
                            ## get number of markers of intscaf
                            nm<-sum(markers$scaff==rownames(scafcar)[intscaf] &
                                        tree$car==as.integer(colnames(scafcar)[i]))
                            if(nm < max(scafsize[fullscaf])){
                                ## changed from '<=' to '<'
                                scafcarviolB[c(intscaf,fullscaf),i]<-1
                            }
                        }
                    }else if(sum(scafcar[,i])>2 & sum(scafcar[,i]>0)==2){
                        ## > 2 internal boundaries across two scaffolds
                        ## might not be a TLcar, but needs further testing
                        ##  (both scaffs with internal boundaries also
                        ##   need end boundaries, and (*,$); -> would have
                        ##   to modify current end boundary test below)
                        ## -> as a simplification, call it TLcar for now
                        scafcarviolB[(scafcar[,i] + endbnd[,i]) > 0,i]<-1
                    }else if(sum(scafcar[,i]>0)>2){
                        ## > 2 scaffolds with internal boundaries
                        ## -> always TLcar
                        scafcarviolB[(scafcar[,i] + endbnd[,i]) > 0,i]<-1
                    }else if(sum(scafcar[,i])==1){
                        ## 1 internal boundary
                        ## -> not a TLcar (perhaps test $)

                    }else if(sum(scafcar[,i])==0){
                        ## only end boundaries
                        ## -> not a TLcar

                    }else{
                        ## check what this can be
                        stop("Unexpected arrangement of CAR boundaries")
                    }
                }
            }
            if(sum(scafcar[,i] + endbnd[,i] > 2)>0){
                ## too many boundaries on scaffold(s)
                ## -> always TL
                scafcarviolW[(scafcar[,i] + endbnd[,i] > 2),i]<-1
            }
        }
    }
    ## (*) test that no circles are needed to fit end pieces together
    ## ($) test that ends can be joined at next lower hierarchy level
    ## (#) depends on relative size



    ## check sets of CARs for locations of boundaries
    ##   need more than only pairwise checks, e.g., violation can
    ##   also result with special arrangement of three (or more)
    ##   CARs on three (or more) scaffolds (i.e., they form circle)

    if(length(myscafs)>1){
        ## identify candidates: these are CARs that have each one part on
        ##  the end of one scaff and on the end of another scaff
        ##  (they should not yet be tagged)
        ##  - two internal boundaries on different scaffs
        tmp1<-which(apply(scafcar,2,function(x) sum(x==1)==2)==TRUE)
        ##  - two end boundaries on different scaffs
        tmp2<-which(apply(endbnd,2,function(x) sum(x==1)==2)==TRUE)
        tmp<-intersect(tmp1,tmp2)
        ##  - exclude those that are already tagged
        ToRm<-which(colSums(scafcarviolB)>0)
        tmp<-setdiff(tmp,ToRm)
        ##  - check that internal and end boundaries are on same two scaffolds
        if(length(tmp)>1){
            setcand<-tmp[(colSums(scafcar[,tmp]*endbnd[,tmp])==2)]
        }else{
            setcand<-integer(0)
        }

        if(length(setcand)>1){
            for(i in 1:length(setcand)){
                ci<-setcand[i]
                si<-which(scafcar[,ci]>0)
                if(length(si)!=2){
                    stop("tagTLcar: wrong number of scaffolds in boundary check")
                }
                if(sum(scafcarviolB[,ci])==0){ ## haven't been tagged already
                    subcand<-setcand[-i]
                    ## also exclude subcand that have tags already
                    if(length(subcand)>1){
                        tmp<-which(colSums(scafcarviolB[,subcand])>0)
                        if(length(tmp)>0){
                            subcand<-subcand[-tmp]
                        }
                    }else if(length(subcand)==1){
                        if(sum(scafcarviolB[,subcand])>0){
                            subcand<-integer()
                        }
                    }
                    ## traverse candidates from one scaff in focal CAR
                    ##  until other scaff is hit, or no more match
                    ##  between the candidates
                    hit<-0
                    tmp<-which(scafcar[si[1],subcand]>0)
                    sx<-si[1]
                    traverse<-tmp

                    while(hit==0 & length(tmp)>0 &
                              length(subcand)>=length(traverse)){
                        ## get scaffolds
                        sy<-which(scafcar[,subcand[tmp]]>0)
                        sy<-setdiff(sy,sx) ## exclude previous scaff
                        if(sy==si[2]){ ## hit
                            hit<-1
                        }else{
                            tmp<-which(scafcar[sy,subcand]>0)
                            tmp<-setdiff(tmp,traverse) ## excl. prev. CARs
                            if(length(tmp)==0){ ## no more match
                                break
                            }
                            sx<-sy
                            traverse<-c(traverse,tmp)
                        }
                    }
                    if(hit==1){
                        sxi<-which((scafcar[,ci]+endbnd[,ci])>0)
                        scafcarviolB[sxi,ci]<-1
                        for(cj in subcand[traverse]){
                            sxj<-which((scafcar[,cj]+endbnd[,cj])>0)
                            scafcarviolB[sxj,cj]<-1
                        }
                    }
                }
            }
        }
    }


    ## ## finally, remove all tags for CARs that fully span a scaffold
    ## ## (>>>> this is for scaffold-level assemblies, where one scaffold might
    ## ##       be placed in middle of other scaffold currently padded with Ns)
    ## scafcarviolB[endbnd==2 & scafcar==0]<-0


    if(sum(scafcarviolB)>0){
        ## go through all cars
        for(j in 1:ncol(scafcarviolB)){
            tl<-rep(0,nrow(markers))
            for(i in 1:nrow(scafcarviolB)){
                ## go through all scaffolds
                if(scafcarviolB[i,j]>0){
                    tmp<-which(markers$scaff==myscafs[i] & tree$car==mycars[j])
                    tl[tmp]<-rep(1,length(tmp))
                }
            }
            ## append column to TLbetween
            if(sum(tl)>0){
                TLbetween<-cbind(TLbetween,tl)
            }
        }
    }

    if(sum(scafcarviolW)>0){
        ## go through all cars
        for(j in 1:ncol(scafcarviolW)){
            tl<-rep(0,nrow(markers))
            for(i in 1:nrow(scafcarviolW)){
                ## go through all scaffolds
                if(scafcarviolW[i,j]>0){
                    tmp<-which(markers$scaff==myscafs[i] & tree$car==mycars[j])
                    tl[tmp]<-rep(1,length(tmp))
                }
            }
            ## append column to TLwithin
            if(sum(tl)>0){
                TLwithin<-cbind(TLwithin,tl)
            }
        }
    }

    return(list(TLbetween=TLbetween,TLwithin=TLwithin))

}

## ------------------------------------------------------------------


## find column in tagMat that corresponds to a particular gap
## >>>> should tagMat always have a single block per scaffold? if yes,
##      test should better be c(1:(insS-1),c(insE+1):nrow(tagMat))
findInsert<-function(tagMat,insS,insE){
    if(!is.matrix(tagMat)){
        stop("Require matrix as input")
    }else if(ncol(tagMat)==0){
        tpcol<-integer(0)
    }else if(insS>insE | insS<2 | insE>=nrow(tagMat)){
        stop("Invalid insert start or end position")
    }else if(ncol(tagMat)==1){
        vtagMat<-as.vector(tagMat)
        if(sum(vtagMat[insS:insE])==length(insS:insE) &
           sum(vtagMat[c(insS-1,insE+1)])==0){
            tpcol<-1
        }else{
            tpcol<-integer(0)
        }
    }else if(insS==insE){
        tpcol<-which(tagMat[insS,]==1 &
                         colSums(tagMat[c(insS-1,insE+1),])==0)
    }else{
        tpcol<-which(colSums(tagMat[insS:insE,])==
                         length(insS:insE) &
                             colSums(tagMat[c(insS-1,insE+1),])==0)
    }
    return(tpcol)
    ## returns integer(0) if no matches
}

## ------------------------------------------------------------------


## obtain size and positions of flanks and gaps from matrix that
##  has flank positions for each duplicated element in a column
getGapsFlanks<-function(FlankMat){

    if(!is.matrix(FlankMat)){
        stop("Require matrix as input")
    }else if(ncol(FlankMat)==0){

        ## Gaps<-list(matrix(0,nrow=3,ncol=0,
        ##                   dimnames=list(c("gapsize","gapstart","gapend"))))
        ## Flanks<-list(matrix(0,nrow=3,ncol=0,
        ##                     dimnames=list(c("flanksize","flankstart",
        ##                         "flankend"))))
        ## return(list(Gaps=Gaps,Flanks=Flanks))

        stop("Require at least one column in FlankMat")

    }else if(ncol(FlankMat)>0){

        Gaps<-vector("list",ncol(FlankMat))
        Flanks<-vector("list",ncol(FlankMat))

        for(i in 1:ncol(FlankMat)){
            tmp<-as.vector(which(FlankMat[,i]>0))
            if(length(tmp)==0){
                stop("Function 'getGapsFlanks' requires non-zero values")
            }
            tmpgaps<-which(diff(tmp)>1)
            ## save gap size
            gapsize<-tmp[tmpgaps+1]-tmp[tmpgaps]-1
            ## get gap positions
            gp<-as.vector(which(FlankMat[,i]==0))
            gp<-gp[gp>min(tmp) & gp<max(tmp)]
            gapend<-gp[c(which(diff(gp)>1),length(gp))]
            if(length(gp)>0){
                gapstart<-gp[c(1,which(diff(gp)>1)+1)]
            }else{
                gapstart<-integer(0)
            }
            Gaps[[i]]<-rbind(gapsize,gapstart,gapend)
            ## get flank positions and sizes
            flankend<-tmp[c(tmpgaps,length(tmp))]
            flankstart<-tmp[c(1,tmpgaps+1)]
            flanksize<-flankend-flankstart+1
            Flanks[[i]]<-rbind(flanksize,flankstart,flankend)
        }

        for(i in 1:ncol(FlankMat)){
            if(ncol(Flanks[[i]])!=ncol(Gaps[[i]])+1){
                stop(paste0("FlankMat in column ",i," has unexpected flank - gap count"))
            }
        }

        return(list(Gaps=Gaps,Flanks=Flanks))

    }else{
        stop("Function 'getGapsFlanks' requires a matrix")
    }
}

## ------------------------------------------------------------------


## filter out large inserts and tag adjacent flanks instead,
##  otherwise, tag inserts only
makeTPtags<-function(TPflanks,TPelemGaps,TPelemFlanks,nelem,remWgt=0.05){


    ## storage for all kept tags of TPflanks
    TPfiltF<-matrix(NA,ncol=0,nrow=nelem) ## flanks (tmp)
    TPkeptF<-matrix(NA,ncol=0,nrow=nelem) ## flanks
    TPremF<-matrix(NA,ncol=0,nrow=nelem) ## removed flanks
    TPkeptI<-matrix(NA,ncol=0,nrow=nelem) ## inserts
    TPremI<-matrix(NA,ncol=0,nrow=nelem) ## removed inserts (tmp)

    ## keep flanks if
    ## - only if flank small relative to insert
    ## keep insert otherwise

    for(i in 1:ncol(TPflanks)){

        if(ncol(TPelemFlanks[[i]])==2){

            if((TPelemGaps[[i]][1,1]>
                    TPelemFlanks[[i]][1,1]) |
               (TPelemGaps[[i]][1,1]>
                    TPelemFlanks[[i]][1,2])){
                ## tag flanks temporarily
                tag<-numeric(nelem)
                tag[(TPelemFlanks[[i]][2,1]):(TPelemFlanks[[i]][3,1])]<-rep(-0.5,TPelemFlanks[[i]][1,1])
                tag[(TPelemFlanks[[i]][2,2]):(TPelemFlanks[[i]][3,2])]<-rep(0.5,TPelemFlanks[[i]][1,2])
                TPfiltF<-cbind(TPfiltF,tag)
                ## store unused insert temporarily
                tag<-numeric(nelem)
                tag[(TPelemGaps[[i]][2,1]):(TPelemGaps[[i]][3,1])]<-rep(1,TPelemGaps[[i]][1,1])
                TPremI<-cbind(TPremI,tag)
            }else{
                ## tag insert
                tag<-numeric(nelem)
                tag[(TPelemGaps[[i]][2,1]):(TPelemGaps[[i]][3,1])]<-rep(1,TPelemGaps[[i]][1,1])
                TPkeptI<-cbind(TPkeptI,tag)
            }

        }else if(ncol(TPelemFlanks[[i]])>2){

            ## >>> this is a bit more complicated so moved elements
            ##     of older function here that might be a bit
            ##     more cumbersome but have been tested already
            ## tmpFlanks<-matrix(0,nrow=nrow(TPflanks),ncol=0)
            tmpGaps<-matrix(0,nrow=nrow(TPflanks),
                            ncol=ncol(TPelemGaps[[i]]))
            GapsToRm<-numeric()

            for(j in 1:ncol(TPelemGaps[[i]])){
                tmpGaps[(TPelemGaps[[i]][2,j]):(TPelemGaps[[i]][3,j]),j]<-rep(1,TPelemGaps[[i]][1,j])
            }

            ## check for large inserts; if any, then for each,
            ##  remove temporary insert tag, and tag both flanks
            ##  temporarily
            lrg<-integer()
            for(j in 1:ncol(TPelemGaps[[i]])){
                if((TPelemGaps[[i]][1,j]>
                        sum(TPelemFlanks[[i]][1,1:j])) |
                   (TPelemGaps[[i]][1,j]>
                        sum(TPelemFlanks[[i]][1,(j+1):ncol(TPelemFlanks[[i]])]))){
                    lrg<-c(lrg,j)
                }
            }
            ## ensures that many smaller inserts
            ##  will kept tagged as inserts

            if(length(lrg)>0){
                for(j in 1:length(lrg)){
                    ## tag flanks temporarily
                    tag<-numeric(nelem)
                    for(k in 1:(lrg[j])){
                        tag[(TPelemFlanks[[i]][2,k]):(TPelemFlanks[[i]][3,k])]<-rep(-0.5,TPelemFlanks[[i]][1,k])
                    }
                    for(k in (lrg[j]+1):ncol(TPelemFlanks[[i]])){
                        tag[(TPelemFlanks[[i]][2,k]):(TPelemFlanks[[i]][3,k])]<-rep(0.5,TPelemFlanks[[i]][1,k])
                    }
                    TPfiltF<-cbind(TPfiltF,tag)
                    ## store insert tag to be removed
                    tmp<-findInsert(tmpGaps,TPelemGaps[[i]][2,lrg[j]],
                                    TPelemGaps[[i]][3,lrg[j]])
                    if(length(tmp)!=1){
                        stop("Unique insert could not be identified")
                    }
                    GapsToRm<-c(GapsToRm,tmp)
                }
            }

            ## finally, remove inserts that should not be tagged from
            ##  tmpGaps and add to TPremI, and others to TPkeptI
            GapsToRm<-unique(as.vector(GapsToRm))
            if(length(GapsToRm)>0){
                TPkeptI<-cbind(TPkeptI,tmpGaps[,-GapsToRm])
                TPremI<-cbind(TPremI,tmpGaps[,GapsToRm])
            }else{
                TPkeptI<-cbind(TPkeptI,tmpGaps)
            }


        }else{ ## close ncol(TPelemFlanks[[i]])>2
            stop("Expected at least two flanks")
        }
    } ## close loop over ncol(TPflanks)

    ## decide which part of flank to keep
    if(ncol(TPfiltF)>0){
        TPkeptF<-matrix(0,ncol=ncol(TPfiltF),nrow=nelem)
        TPremF<-matrix(0,ncol=ncol(TPfiltF),nrow=nelem)
        for(i in 1:ncol(TPfiltF)){
            ## uncovered left of insert
            tmpL<-which(TPfiltF[,i]== -0.5)
            uncL<-sum(rowSums(TPkeptI[tmpL,,drop=FALSE])<1)
            ## uncovered right of insert
            tmpR<-which(TPfiltF[,i]==0.5)
            uncR<-sum(rowSums(TPkeptI[tmpR,,drop=FALSE])<1)
            if(uncL>0 & uncR>0){
                ## decide which flank to keep
                ##  - always keep smaller one
                ##    (>>>> this is different from algorithm at
                ##          the CAR level, where decision also
                ##          depends on whether one flank is contained
                ##          in the other for a Q-node at the next
                ##          level of hierarchy, assuming here
                ##          that P-nodes are rare and commonly
                ##          contain single-marker Q-nodes)

                ## keep smaller flank
                if(length(tmpL)<length(tmpR)){
                    TPkeptF[tmpL,i]<-rep(1,length(tmpL))
                    TPremF[tmpR,i]<-rep(1,length(tmpR))
                }else if(length(tmpR)<length(tmpL)){
                    TPkeptF[tmpR,i]<-rep(1,length(tmpR))
                    TPremF[tmpL,i]<-rep(1,length(tmpL))
                }else{
                    ## keep both flanks tagged with 0.5
                    TPkeptF[c(tmpR,tmpL),i]<-rep(0.5,length(c(tmpR,tmpL)))
                }
            }else if(uncL>0 & uncR==0){
                ## keep right flank
                TPkeptF[tmpR,i]<-rep(1,length(tmpR))
                TPremF[tmpL,i]<-rep(1,length(tmpL))
            }else if(uncL==0 & uncR>0){
                ## keep left flank
                TPkeptF[tmpL,i]<-rep(1,length(tmpL))
                TPremF[tmpR,i]<-rep(1,length(tmpR))
            }else{
                ## keep both flanks tagged with 0.5, or remove both
                ##  (make sure that when removing, matrix dimensions
                ##   need to stay the same for TPkeptF, TPremF, TPremI,
                ##   and TPfiltF, as same dimensions are assumed below)
                TPkeptF[c(tmpR,tmpL),i]<-rep(0.5,length(c(tmpR,tmpL)))
                ##TPremF[c(tmpR,tmpL),i]<-rep(1,length(c(tmpR,tmpL)))
            }
        } ## close loop over ncol(TPfiltF)
    } ## close if(ncol(TPfiltF)>0){

    ## test for immediate adjacencies of flanks with 1.0 tags
    ##  (probably rare event >>>> would need some more testing)
    if(ncol(TPkeptF)>1){
        ## get columns that have 1.0 tags
        totest<-which(apply(TPkeptF, 2, function(x) sum(x==1)>0)==TRUE)
        if(length(totest)>1){
            toRm<-integer()
            toKeep<-integer()
            for(m in 1:(length(totest)-1)){
                for(w in (m+1):length(totest)){
                    mgaps<-sum(diff(which(TPkeptF[,totest[m]]==1))>1)
                    wgaps<-sum(diff(which(TPkeptF[,totest[w]]==1))>1)
                    tmp<-pmax(TPkeptF[,totest[m]],TPkeptF[,totest[w]])
                    tgaps<-sum(diff(which(tmp==1))>1)
                    ## require that no gap will be removed in combined
                    ##  vector, requiring that one or both single vectors
                    ##  have zero gaps (otherwise things might cancel out)
                    if((mgaps + wgaps == tgaps) & (mgaps==0 | wgaps==0)){
                        ## test that there is no overlap of flank tags
                        if(sum(rowSums(TPkeptF[,totest[c(m,w)]])>1)==0){
                            ## check if exactly one has removed insert tag
                            ##  (automatically takes into account that zero
                            ##   gaps are allowed)
                            tmp1<-TPremI[,totest[m]]*TPkeptF[,totest[w]]
                            tmp2<-TPremI[,totest[w]]*TPkeptF[,totest[m]]
                            if(sum(tmp1>0)==sum(TPremI[,totest[w]]>0) &
                               sum(tmp2)==0){
                                toRm<-c(toRm,totest[m])
                                toKeep<-c(toKeep,totest[w])
                            }else if(sum(tmp2>0)==sum(TPremI[,totest[m]]>0) &
                                     sum(tmp1)==0){
                                toRm<-c(toRm,totest[w])
                                toKeep<-c(toKeep,totest[m])
                            }
                        }
                    }
                }
            }
            ## make sure no flanks will be removed that are the ones
            ##  that should be kept in another pair
            toRm<-toRm[!is.element(toRm,toKeep) & !is.element(toKeep,toRm)]
            if(length(toRm)>0){
                for(i in unique(toRm)){
                    ngaps<-sum(diff(which(TPkeptF[,i]!=0))>1)
                    if(ngaps>0){
                        stop("Removed flank contained a gap")
                    }
                    TPremF[,i]<-TPremF[,i]+TPkeptF[,i]
                    TPkeptF[,i]<-0
                }
                TPfiltF<-TPfiltF[,-unique(toRm),drop=FALSE]
            }
        }
    }

    ## test for immediate adjacencies of inserts (always 1.0 tags)
    ##  (probably rare event >>>> would need some more testing)
    if(ncol(TPkeptI)>1){
        toRm<-integer()
        toKeep<-integer()
        toKeepBoth<-integer()
        for(m in 1:(ncol(TPkeptI)-1)){
            for(w in (m+1):ncol(TPkeptI)){
                mgaps<-sum(diff(which(TPkeptI[,m]==1))>1)
                wgaps<-sum(diff(which(TPkeptI[,w]==1))>1)
                tmp<-pmax(TPkeptI[,m],TPkeptI[,w])
                tgaps<-sum(diff(which(tmp==1))>1)
                ## require that zero gaps in combined and both
                ##  single vectors, and no overlap of tags
                if((mgaps + wgaps + tgaps)==0 &
                   sum(rowSums(TPkeptI[,c(m,w)])>1)==0){
                    ## check that there is no other insert directly
                    ##  adjacent to these two to make sure that
                    ##  no inserts will be removed that
                    ##  might need to be kept in another pair
                    insertends<-range(which(tmp==1))
                    a1<-which(TPkeptI[insertends[1],]==0 &
                                  TPkeptI[insertends[1]-1,]==1)
                    a2<-which(TPkeptI[insertends[2],]==0 &
                                  TPkeptI[insertends[2]+1,]==1)
                    if(length(c(a1,a2))==0){
                        ## get size of inserts, keep smaller insert
                        nm<-sum(TPkeptI[,m]==1)
                        nw<-sum(TPkeptI[,w]==1)
                        if(nm > nw){
                            toRm<-c(toRm,m)
                            toKeep<-c(toKeep,w)
                        }else if(nm < nw){
                            toRm<-c(toRm,w)
                            toKeep<-c(toKeep,m)
                        }else{
                            toKeepBoth<-c(toKeepBoth,c(m,w))
                        }
                    }
                }
            }
        }
        toRm<-unique(toRm)
        toKeep<-unique(toKeep)
        toKeepBoth<-unique(toKeepBoth)
        if(length(unique(c(toRm,toKeep,toKeepBoth)))<
           length(c(toRm,toKeep,toKeepBoth))){
            stop("Removed/kept inserts overlap")
        }
        ## adjust weight for removed/kept inserts
        ##  (will work if column indices are numeric())
        TPkeptI[,toRm]<-TPkeptI[,toRm]*remWgt
        TPkeptI[,toKeep]<-TPkeptI[,toKeep]*(1-remWgt)
        TPkeptI[,toKeepBoth]<-TPkeptI[,toKeepBoth]*0.5
    }


    ## bind TPkeptI and TPkeptF together, adjusting for weight
    ##  of removed flanks
    TPkeptF[TPkeptF==1]<-1-remWgt ## don't adjust 0.5 tags
    TPremF[TPremF==1]<-remWgt
    TPfilt<-cbind(TPkeptI,TPkeptF+TPremF)

    ## make subnode IDs
    subtags<-integer(nelem)
    if(ncol(TPkeptI)+ncol(TPfiltF)>0){
        tagmat<-matrix(NA,ncol=0,nrow=nelem)
        ## flanks
        if(ncol(TPfiltF)>0){
            for(i in 1:ncol(TPfiltF)){
                ## left of insert
                tmpL<-which(TPfiltF[,i]== -0.5)
                ## right of insert
                tmpR<-which(TPfiltF[,i]==0.5)
                ## make different tags
                tag<-numeric(nelem)
                tag[tmpL]<-rep(1,length(tmpL))
                tag[tmpR]<-rep(2,length(tmpR))
                tagmat<-cbind(tagmat,tag)
            }
        }
        ## inserts
        tmpmat<-TPkeptI
        if(remWgt>0){
            tmpmat[tmpmat<=remWgt]<-0
            tmpmat[tmpmat>=0.5]<-1
            ## unlike for flanks, setting all 0.5 to 1 is okay
            ##  as they are always in different columns
        }
        tagmat<-cbind(tagmat,tmpmat)

        tmptags<-apply(tagmat,1,function(x) paste0(x,collapse="-"))
        unitags<-unique(tmptags)
        ## exclude markers that are untagged
        notag<-paste0(rep(0,ncol(tagmat)),collapse="-")
        unitags<-unitags[unitags!=notag]
        for(i in 1:length(unitags)){
            tagpos<-which(tmptags==unitags[i])
            subtags[tagpos]<-rep(i,length(tagpos))
        }
    }

    return(list(TPfilt=TPfilt,subtags=subtags))
    ## matrix with tags for inserts/gaps and subnode IDs
}

## ------------------------------------------------------------------


## main loop for going through all scaffolds to identify rearrangements
tagRearr<-function(SYNT,markers,tree,orientation,nhier,myscafs,splitnodes,
                   remWgt=0.05,testlim=100){

    if(remWgt<0 | remWgt>=0.5){
        stop("Tagging weight for removed flanks needs to be [0,0.5)")
    }

    for(s in 1:length(myscafs)){
        myset<-which(markers$scaff==myscafs[s])
        subtree<-tree[myset,]
        suborientation<-orientation[myset]

        ## go through levels of hierarchy and test for synteny compliance
        TLWC<-matrix(NA,ncol=0,nrow=nrow(subtree))
        rownames(TLWC)<-subtree$marker
        IV<-matrix(NA,ncol=0,nrow=nrow(subtree))
        rownames(IV)<-subtree$marker
        nodeori<-matrix(NA,nrow=nrow(subtree),ncol=nhier)
        rownames(nodeori)<-subtree$marker
        subnode<-matrix(0,nrow=nrow(subtree),ncol=nhier)
        rownames(subnode)<-subtree$marker
        ## pass subnode ids from CAR level
        subnode[,1]<-SYNT$subnode[myset,1]

        synt<-list(TLWC=TLWC,IV=IV,TLWCbS=TLWC,TLWCbE=TLWC,IVbS=IV,
                   IVbE=IV,nodeori=nodeori,blockori=nodeori,
                   blockid=nodeori,premask=nodeori,subnode=subnode)
        rm(TLWC,IV,nodeori,subnode)

        ## ----
        ## CAR level already done
        ## ----

        ## ----
        ## Next levels past CAR level:
        ## determine rearrangements, taking orientation into account
        ## ---------------------------------------------------------
        for(n in 2:nhier){
            mycol<-2+(n-1)*2 ## column in subtree with node type
            if(sum(!is.na(subtree[,mycol]))==0){
                break
            }
            hiercol<-seq(3,mycol-1,2) ## make unique tags for preceding hierarchy
            ## * make vector with all preceding hierarchies
            ## * go through unique nodes/genes in that hierarchy
            ## * identify rows in subtree for each, and subset accordingly
            if(splitnodes==TRUE){
                ## take subnode ids into account
                splitids<-matrix(NA,ncol=0,nrow=length(myset))
                for(d in 1:(n-1)){
                    splitids<-cbind(splitids,
                                    paste(subtree[,hiercol[d]],
                                          synt$subnode[,d],sep="."))
                }
                allprev<-apply(splitids,1,
                               function(x) paste0(x,collapse="-"))
            }else{
                allprev<-apply(as.matrix(subtree[,hiercol]),1,
                               function(x) paste0(x,collapse="-"))
            }
            uniprev<-unique(allprev[!is.na(subtree[,mycol])])
            if(length(uniprev)==0){ ## only NAs left
                next
            }
            for(p in 1:length(uniprev)){
                tmprows<-which(allprev==uniprev[p])
                tmptree<-subtree[tmprows,]

                node<-"NA"
                currelem<-tmptree[1,mycol+1]
                if(is.na(currelem)){currelem<-numeric()}
                allelem<-currelem
                elemrows<-matrix(1,nrow=2,ncol=length(allelem))
                TPelem<-matrix(0,nrow=0,ncol=length(allelem))
                for(i in 1:nrow(tmptree)){
                    ## the following is for either new preceding blocks or
                    ##  continuing existing preceding blocks
                    node<-tmptree[i,mycol]
                    if(!is.element(node,c("P","Q","NA")) & !is.na(node)){
                        stop("Unrecognized node type ",node)
                    }
                    if(is.na(tmptree[i,mycol+1])){
                        currelem<-numeric()
                        allelem<-numeric()
                        elemrows<-matrix(1,nrow=2,ncol=0)
                        TPelem<-matrix(0,nrow=0,ncol=0)
                        ## next
                    }else if(tmptree[i,mycol+1]==currelem){ ## same level
                        elemrows[2,ncol(elemrows)]<-i ## save new maximum row
                    }else{ ## new or already existing level
                        currelem<-subtree[tmprows[i],mycol+1]
                        allelem<-c(allelem,currelem)
                        elemrows<-cbind(elemrows,c(i,i))
                        TPelem<-cbind(TPelem,matrix(0,ncol=1,nrow=nrow(TPelem)))
                        if(!is.element(tmptree[i,mycol+1],allelem[-length(allelem)])){
                            ## -> new level, nothing to be done here
                        }else{ ## level already present
                            ## identify problematic element and tag with 1,
                            ##  and tag inserted element with 2 (accounting
                            ##  for potential nestedness by adding rows)
                            ##  -> needs additional processing
                            TPelem<-rbind(TPelem,matrix(0,nrow=1,
                                                        ncol=ncol(TPelem)))
                            tpcol<-which(allelem==currelem)
                            ## only take last two occurences
                            tpcol<-tail(tpcol,n=2L)
                            tprow<-nrow(TPelem)
                            TPelem[tprow,tpcol]<-c(1,1)
                            TPelem[tprow,(tpcol[1]+1):(tpcol[2]-1)]<-2
                        }
                    }
                }
                ## summarize block
                if(!is.na(node) & node!="NA"){
                    leaves<-tagLeaf(subtree,tmprows,n)
                    testorientation<-suborientation[tmprows]
                    if(node=="Q" & nrow(TPelem)>0){
                        ## ## requires additional processing
                        ## preMasks<-makePreMasks(tmptree,allelem,elemrows,
                        ##                        TPelem,n,nhier)
                        ## >>>> preMasks are currently unused in tagTP2()
                        ##      -> don't compute to save time
                        preMasks<-list(A=rep(FALSE,length(allelem)),
                                       D=rep(FALSE,length(allelem)))
                    }else{
                        preMasks<-list(A=rep(FALSE,length(allelem)),
                                       D=rep(FALSE,length(allelem)))
                    }
                    synt<-tagTP2(synt,allelem,tmprows,elemrows,TPelem,n,node,
                                 leaves,testorientation,preMasks,
                                 splitnodes,remWgt=remWgt,testlim=testlim)
                }
            } ## end loop over p in 1:length(uniprev)
        } ## end loop over n in 2:nhier

        ## adjust column dimensions for TLWC and IV and their
        ##  breakpoints and append data to SYNT
        SYNT$TLWC<-appendData(SYNT$TLWC,synt$TLWC)
        SYNT$IV<-appendData(SYNT$IV,synt$IV)
        SYNT$TLWCbS<-appendData(SYNT$TLWCbS,synt$TLWCbS)
        SYNT$TLWCbE<-appendData(SYNT$TLWCbE,synt$TLWCbE)
        SYNT$IVbS<-appendData(SYNT$IVbS,synt$IVbS)
        SYNT$IVbE<-appendData(SYNT$IVbE,synt$IVbE)

        ## add nodeori, blockori, blockid, premask
        SYNT$nodeori<-rbind(SYNT$nodeori,synt$nodeori)
        SYNT$blockori<-rbind(SYNT$blockori,synt$blockori)
        SYNT$blockid<-rbind(SYNT$blockid,synt$blockid)
        SYNT$premask<-rbind(SYNT$premask,synt$premask)

        ## add subnode ID (without overwriting IDs for CARs)
        SYNT$subnode[myset,2:nhier]<-synt$subnode[,2:nhier]

    } ## close loop over s in 1:length(myscafs)

    return(SYNT)
}

## ------------------------------------------------------------------


## filter out large inserts at CAR level
filterCars3<-function(SYNT,markers,tree,nhier,myscafs,mycars,TLL,
                      scafcarbest,bestHitOpt=1,remWgt=0.05){
    ## obtain size and positions of pieces that are inserted
    ##  note: there might only be a single piece of a fragmented CAR
    ##  on one scaffold (will return matrix with 0 ncol)
    ## obtain size and positions of flanks

    if(!is.matrix(TLL$TLbetween) | !is.matrix(TLL$TLwithin)){
        stop("Require matrix as input")
    }
    if(ncol(TLL$TLbetween)==0){
        SYNT$TLBS<-matrix(NA,ncol=0,nrow=nrow(markers))
        rownames(SYNT$TLBS)<-rownames(TLL$TLbetween)
    }
    if(ncol(TLL$TLwithin)==0){
        SYNT$TLWS<-matrix(NA,ncol=0,nrow=nrow(markers))
        rownames(SYNT$TLWS)<-rownames(TLL$TLwithin)
    }
    if(ncol(TLL$TLbetween)==0 & ncol(TLL$TLwithin)==0){
        return(SYNT)
    }
    if(!is.element(bestHitOpt,1:3)){
        stop(paste("Invalid bestHitOpt",bestHitOpt))
    }
    if(remWgt<0 | remWgt>=0.5){
        stop("Tagging weight for removed flanks needs to be [0,0.5)")
    }

    ## go through all scaffolds
    ## ------------------------
    for(s in 1:length(myscafs)){
        myset<-which(markers$scaff==myscafs[s])

        if(length(myset)==0){
            stop(paste("Scaffold",myscafs[s],"has no marker"))
        }
        TLWS<-TLL$TLwithin[myset,,drop=FALSE]
        TLBS<-TLL$TLbetween[myset,,drop=FALSE]

        ## remove columns that are zero in TLWS
        tmp<-which(colSums(TLWS)!=0)
        TLWS<-TLWS[,tmp,drop=FALSE]

        ## storage for all kept tags TLWS
        TLWfiltF<-matrix(NA,ncol=0,nrow=length(myset)) ## flanks (tmp)
        TLWkeptF<-matrix(NA,ncol=0,nrow=length(myset)) ## flanks
        TLWremF<-matrix(NA,ncol=0,nrow=length(myset)) ## removed flanks
        TLWkeptI<-matrix(NA,ncol=0,nrow=length(myset)) ## inserts
        TLWremI<-matrix(NA,ncol=0,nrow=length(myset)) ## removed inserts (tmp)
        ## to modify in TLBS
        toMod<-integer()

        ## Translocations within scaffolds
        ## -------------------------------
        ## go through all flanks and make tags for either flanks or inserts
        if(isTRUE(ncol(TLWS)>0)){

            GapFlank<-getGapsFlanks(TLWS)
            ## returns $Gaps and $Flanks

            ## keep flanks only if flank small relative to insert
            ## keep insert otherwise

            for(i in 1:ncol(TLWS)){

                flankCar<-unique(tree$car[myset][TLWS[,i]>0])
                if(length(flankCar)!=1){
                    stop("Flanking CAR is not unique")
                }
                if(ncol(GapFlank$Flanks[[i]])==2){

                    if((GapFlank$Gaps[[i]][1,1]>
                            GapFlank$Flanks[[i]][1,1]) |
                       (GapFlank$Gaps[[i]][1,1]>
                            GapFlank$Flanks[[i]][1,2])){
                        ## tag flanks temporarily
                        tag<-numeric(length(myset))
                        tag[(GapFlank$Flanks[[i]][2,1]):(GapFlank$Flanks[[i]][3,1])]<-rep(-0.5,GapFlank$Flanks[[i]][1,1])
                        tag[(GapFlank$Flanks[[i]][2,2]):(GapFlank$Flanks[[i]][3,2])]<-rep(0.5,GapFlank$Flanks[[i]][1,2])
                        TLWfiltF<-cbind(TLWfiltF,tag)
                        ## store unused insert temporarily
                        tag<-numeric(length(myset))
                        tag[(GapFlank$Gaps[[i]][2,1]):(GapFlank$Gaps[[i]][3,1])]<-rep(1,GapFlank$Gaps[[i]][1,1])
                        TLWremI<-cbind(TLWremI,tag)
                    }else{
                        ## tag insert
                        tag<-numeric(length(myset))
                        tag[(GapFlank$Gaps[[i]][2,1]):(GapFlank$Gaps[[i]][3,1])]<-rep(1,GapFlank$Gaps[[i]][1,1])
                        TLWkeptI<-cbind(TLWkeptI,tag)
                    }

                }else if(ncol(GapFlank$Flanks[[i]])>2){

                    ## >>> this is a bit more complicated so moved elements
                    ##     of older function here that might be a bit
                    ##     more cumbersome but have been tested already
                    ## tmpFlanks<-matrix(0,nrow=nrow(TLWS),ncol=0)
                    tmpGaps<-matrix(0,nrow=nrow(TLWS),
                                    ncol=ncol(GapFlank$Gaps[[i]]))
                    GapsToRm<-numeric()

                    for(j in 1:ncol(GapFlank$Gaps[[i]])){
                        tmpGaps[(GapFlank$Gaps[[i]][2,j]):(GapFlank$Gaps[[i]][3,j]),j]<-rep(1,GapFlank$Gaps[[i]][1,j])
                    }

                    ## check for large inserts; if any, then for each,
                    ##  remove temporary insert tag, and tag both flanks
                    ##  temporarily
                    lrg<-integer()
                    for(j in 1:ncol(GapFlank$Gaps[[i]])){
                        if((GapFlank$Gaps[[i]][1,j]>
                                sum(GapFlank$Flanks[[i]][1,1:j])) |
                           (GapFlank$Gaps[[i]][1,j]>
                                sum(GapFlank$Flanks[[i]][1,(j+1):ncol(GapFlank$Flanks[[i]])]))){
                            lrg<-c(lrg,j)
                        }
                    }
                    ## ensures that many smaller inserts
                    ##  will kept tagged as inserts

                    if(length(lrg)>0){
                        for(j in 1:length(lrg)){
                            ## tag flanks temporarily
                            tag<-numeric(length(myset))
                            for(k in 1:(lrg[j])){
                                tag[(GapFlank$Flanks[[i]][2,k]):(GapFlank$Flanks[[i]][3,k])]<-rep(-0.5,GapFlank$Flanks[[i]][1,k])
                            }
                            for(k in (lrg[j]+1):ncol(GapFlank$Flanks[[i]])){
                                tag[(GapFlank$Flanks[[i]][2,k]):(GapFlank$Flanks[[i]][3,k])]<-rep(0.5,GapFlank$Flanks[[i]][1,k])
                            }
                            TLWfiltF<-cbind(TLWfiltF,tag)
                            ## store insert tag to be removed
                            tmp<-findInsert(tmpGaps,GapFlank$Gaps[[i]][2,lrg[j]],
                                            GapFlank$Gaps[[i]][3,lrg[j]])
                            if(length(tmp)!=1){
                                stop("Unique insert could not be identified")
                            }
                            GapsToRm<-c(GapsToRm,tmp)
                        }
                    }

                    ## finally, remove inserts that should not be tagged from
                    ##  tmpGaps and add to TLWremI, and others to TLWkeptI
                    GapsToRm<-unique(as.vector(GapsToRm))
                    if(length(GapsToRm)>0){
                        TLWkeptI<-cbind(TLWkeptI,tmpGaps[,-GapsToRm])
                        TLWremI<-cbind(TLWremI,tmpGaps[,GapsToRm])
                    }else{
                        TLWkeptI<-cbind(TLWkeptI,tmpGaps)
                    }


                }else{ ## close ncol(GapFlank$Flanks[[i]])>2
                    stop("Expected at least two flanks")
                }
            } ## close loop over ncol(TLWS)

            ## decide which part of flank to keep
            if(ncol(TLWfiltF)>0){
                TLWkeptF<-matrix(0,ncol=ncol(TLWfiltF),nrow=length(myset))
                TLWremF<-matrix(0,ncol=ncol(TLWfiltF),nrow=length(myset))
                for(i in 1:ncol(TLWfiltF)){
                    ## uncovered left of insert
                    tmpL<-which(TLWfiltF[,i]== -0.5)
                    uncL<-sum(rowSums(TLWkeptI[tmpL,,drop=FALSE])<1)
                    ## uncovered right of insert
                    tmpR<-which(TLWfiltF[,i]==0.5)
                    uncR<-sum(rowSums(TLWkeptI[tmpR,,drop=FALSE])<1)
                    if(uncL>0 & uncR>0){
                        ## decide which flank to keep
                        ##  - one contained in other for Q-node:
                        ##    keep contained one
                        ##  - else: keep smaller one
                        m<-2
                        cont<-0
                        eLp<-integer(length(tmpL))
                        eRp<-integer(length(tmpR))
                        while(cont==0 & m<=nhier){
                            if(is.element("Q",tree[myset[tmpL],2+(m-1)*2]) &
                               is.element("Q",tree[myset[tmpR],2+(m-1)*2])){
                                eL<-tree[myset[tmpL],2+(m-1)*2+1]
                                eR<-tree[myset[tmpR],2+(m-1)*2+1]
                                eL2<-eL[!is.na(eL)]
                                eR2<-eR[!is.na(eR)]
                                ## potentially pad elements, if needed
                                ## (>>>> might need more testing)
                                p<-1
                                while(sum(c(eL2,eR2)%%p < c(eL2,eR2))>0){
                                    p<-p*10
                                }
                                eL<-eL+p + eLp*p*10
                                eR<-eR+p + eRp*p*10
                                eLp<-eL
                                eRp<-eR
                                eL<-eL[!is.na(eL)]
                                eR<-eR[!is.na(eR)]
                                if(length(eL)>0 & length(eR)>0){
                                    if(min(eL)<min(eR) & max(eL)>max(eR)){
                                        ## R contained in L
                                        cont<-1
                                    }else if(min(eR)<min(eL) & max(eR)>max(eL)){
                                        ## L contained in R
                                        cont<-2
                                    }else if(max(eL)<min(eR) | min(eL)>max(eR)){
                                        ## neither contained in other
                                        cont<- -1
                                    } ## else: overlap -> test next level
                                }else{
                                    cont<- -1
                                }
                            }else{
                                cont<- -1
                            }
                            m<-m+1
                        }
                        if(cont==1){ ## R contained in L
                            ## keep right flank
                            TLWkeptF[tmpR,i]<-rep(1,length(tmpR))
                            TLWremF[tmpL,i]<-rep(1,length(tmpL))
                        }else if(cont==2){ ## L contained in R
                            ## keep left flank
                            TLWkeptF[tmpL,i]<-rep(1,length(tmpL))
                            TLWremF[tmpR,i]<-rep(1,length(tmpR))
                        }else{
                            ## keep smaller flank
                            if(length(tmpL)<length(tmpR)){
                                TLWkeptF[tmpL,i]<-rep(1,length(tmpL))
                                TLWremF[tmpR,i]<-rep(1,length(tmpR))
                            }else if(length(tmpR)<length(tmpL)){
                                TLWkeptF[tmpR,i]<-rep(1,length(tmpR))
                                TLWremF[tmpL,i]<-rep(1,length(tmpL))
                            }else{
                                ## keep both flanks tagged with 0.5
                                TLWkeptF[c(tmpR,tmpL),i]<-rep(0.5,length(c(tmpR,tmpL)))
                            }
                        }
                    }else if(uncL>0 & uncR==0){
                        ## keep right flank
                        TLWkeptF[tmpR,i]<-rep(1,length(tmpR))
                        TLWremF[tmpL,i]<-rep(1,length(tmpL))
                    }else if(uncL==0 & uncR>0){
                        ## keep left flank
                        TLWkeptF[tmpL,i]<-rep(1,length(tmpL))
                        TLWremF[tmpR,i]<-rep(1,length(tmpR))
                    }else{
                        ## keep both flanks tagged with 0.5, or remove both
                        ##  (make sure that when removing, matrix dimensions
                        ##   need to stay the same for TLWkeptF, TLWremF, TLWremI,
                        ##   and TLWfiltF, as same dimensions are assumed below)
                        TLWkeptF[c(tmpR,tmpL),i]<-rep(0.5,length(c(tmpR,tmpL)))
                    }
                } ## close loop over ncol(TLWfiltF)
            } ## close if(ncol(TLWfiltF)>0){
        } ## close if(isTRUE(ncol(TLWS)>0)){

        ## test for immediate adjacencies of flanks with 1.0 tags
        ##  (probably rare event >>>> would need some more testing)
        if(ncol(TLWkeptF)>1){
            ## get columns that have 1.0 tags
            totest<-which(apply(TLWkeptF, 2, function(x) sum(x==1)>0)==TRUE)
            if(length(totest)>1){
                toRm<-integer()
                toKeep<-integer()
                for(m in 1:(length(totest)-1)){
                    for(w in (m+1):length(totest)){
                        mgaps<-sum(diff(which(TLWkeptF[,totest[m]]==1))>1)
                        wgaps<-sum(diff(which(TLWkeptF[,totest[w]]==1))>1)
                        tmp<-pmax(TLWkeptF[,totest[m]],TLWkeptF[,totest[w]])
                        tgaps<-sum(diff(which(tmp==1))>1)
                        ## require that no gap will be removed in combined
                        ##  vector, requiring that one or both single vectors
                        ##  have zero gaps (otherwise things might cancel out)
                        if((mgaps + wgaps == tgaps) & (mgaps==0 | wgaps==0)){
                            ## test that there is no overlap of flank tags
                            if(sum(rowSums(TLWkeptF[,totest[c(m,w)]])>1)==0){
                                ## check if exactly one has removed insert tag
                                ##  (automatically takes into account that zero
                                ##   gaps are allowed)
                                tmp1<-TLWremI[,totest[m]]*TLWkeptF[,totest[w]]
                                tmp2<-TLWremI[,totest[w]]*TLWkeptF[,totest[m]]
                                if(sum(tmp1>0)==sum(TLWremI[,totest[w]]>0) &
                                   sum(tmp2)==0){
                                    toRm<-c(toRm,totest[m])
                                    toKeep<-c(toKeep,totest[w])
                                }else if(sum(tmp2>0)==sum(TLWremI[,totest[m]]>0) &
                                         sum(tmp1)==0){
                                    toRm<-c(toRm,totest[w])
                                    toKeep<-c(toKeep,totest[m])
                                }
                            }
                        }
                    }
                }
                ## make sure no flanks will be removed that are the ones
                ##  that should be kept in another pair
                toRm<-toRm[!is.element(toRm,toKeep) & !is.element(toKeep,toRm)]
                ## this double-check is needed because these are pairs
                if(length(toRm)>0){
                    for(i in unique(toRm)){
                        ngaps<-sum(diff(which(TLWkeptF[,i]!=0))>1)
                        if(ngaps>0){
                            stop("Removed flank contained a gap")
                        }
                        TLWremF[,i]<-TLWremF[,i]+TLWkeptF[,i]
                        TLWkeptF[,i]<-0
                    }
                    TLWfiltF<-TLWfiltF[,-unique(toRm),drop=FALSE]
                }
            }
        }

        ## test for immediate adjacencies of inserts (always 1.0 tags)
        ##  (probably rare event >>>> would need some more testing)
        if(ncol(TLWkeptI)>1){
            toRm<-integer()
            toKeep<-integer()
            toKeepBoth<-integer()
            for(m in 1:(ncol(TLWkeptI)-1)){
                for(w in (m+1):ncol(TLWkeptI)){
                    mgaps<-sum(diff(which(TLWkeptI[,m]==1))>1)
                    wgaps<-sum(diff(which(TLWkeptI[,w]==1))>1)
                    tmp<-pmax(TLWkeptI[,m],TLWkeptI[,w])
                    tgaps<-sum(diff(which(tmp==1))>1)
                    ## require that zero gaps in combined and both
                    ##  single vectors, and no overlap of tags
                    if((mgaps + wgaps + tgaps)==0 &
                       sum(rowSums(TLWkeptI[,c(m,w)])>1)==0){
                        ## check that there is no other insert directly
                        ##  adjacent to these two to make sure that
                        ##  no inserts will be removed that
                        ##  might need to be kept in another pair
                        insertends<-range(which(tmp==1))
                        a1<-which(TLWkeptI[insertends[1],]==0 &
                                      TLWkeptI[insertends[1]-1,]==1)
                        a2<-which(TLWkeptI[insertends[2],]==0 &
                                      TLWkeptI[insertends[2]+1,]==1)
                        if(length(c(a1,a2))==0){
                            ## get size of inserts, keep smaller insert
                            nm<-sum(TLWkeptI[,m]==1)
                            nw<-sum(TLWkeptI[,w]==1)
                            if(nm > nw){
                                toRm<-c(toRm,m)
                                toKeep<-c(toKeep,w)
                            }else if(nm < nw){
                                toRm<-c(toRm,w)
                                toKeep<-c(toKeep,m)
                            }else{
                                toKeepBoth<-c(toKeepBoth,c(m,w))
                            }
                        }
                    }
                }
            }
            toRm<-unique(toRm)
            toKeep<-unique(toKeep)
            toKeepBoth<-unique(toKeepBoth)
            if(length(unique(c(toRm,toKeep,toKeepBoth)))<
               length(c(toRm,toKeep,toKeepBoth))){
                stop("Removed/kept inserts overlap")
            }
            ## adjust weight for removed/kept inserts
            ##  (will work if column indices are numeric())
            TLWkeptI[,toRm]<-TLWkeptI[,toRm]*remWgt
            TLWkeptI[,toKeep]<-TLWkeptI[,toKeep]*(1-remWgt)
            TLWkeptI[,toKeepBoth]<-TLWkeptI[,toKeepBoth]*0.5
        }


        ## bind TLWkeptI and TLWkeptF together, adjusting for weight
        ##  of removed/kept flanks
        TLWkeptF[TLWkeptF==1]<-1-remWgt ## don't adjust 0.5 tags
        TLWremF[TLWremF==1]<-remWgt
        TLWfilt<-cbind(TLWkeptI,TLWkeptF+TLWremF)

        ## make subnode IDs
        if(ncol(TLWkeptI)+ncol(TLWfiltF)>0){
            tagmat<-matrix(NA,ncol=0,nrow=length(myset))
            ## flanks
            if(ncol(TLWfiltF)>0){
                for(i in 1:ncol(TLWfiltF)){
                    ## left of insert
                    tmpL<-which(TLWfiltF[,i]== -0.5)
                    ## right of insert
                    tmpR<-which(TLWfiltF[,i]==0.5)
                    ## make different tags
                    tag<-numeric(length(myset))
                    tag[tmpL]<-rep(1,length(tmpL))
                    tag[tmpR]<-rep(2,length(tmpR))
                    tagmat<-cbind(tagmat,tag)
                }
            }
            ## inserts
            tmpmat<-TLWkeptI
            if(remWgt>0){
                tmpmat[tmpmat<=remWgt]<-0
                tmpmat[tmpmat>=0.5]<-1
                ## unlike for flanks, setting all 0.5 to 1 is okay
                ##  as they are always in different columns
            }
            tagmat<-cbind(tagmat,tmpmat)

            subtags<-integer(length(myset))
            tmptags<-apply(tagmat,1,function(x) paste0(x,collapse="-"))
            unitags<-unique(tmptags)
            ## exclude markers that are untagged
            notag<-paste0(rep(0,ncol(tagmat)),collapse="-")
            unitags<-unitags[unitags!=notag]
            for(i in 1:length(unitags)){
                tagpos<-which(tmptags==unitags[i])
                subtags[tagpos]<-rep(i,length(tagpos))
            }
            SYNT$subnode[myset,1]<-subtags
        }

        ## adjust column dimensions and append data to SYNT
        SYNT$TLWS<-appendData(SYNT$TLWS,TLWfilt)

        ## obtain breakpoints for TLWcar
        bptW<-getBreakpntsSE(TLWfilt)
        ## returns $bptS and $bptE
        ## adjust column dimensions and append data to SYNT
        SYNT$TLWSbS<-appendData(SYNT$TLWSbS,bptW$bptS)
        SYNT$TLWSbE<-appendData(SYNT$TLWSbE,bptW$bptE)


        ## Translocations between scaffolds
        ## --------------------------------
        if(isTRUE(ncol(TLBS)>0)){

            TLcand<-as.vector(which(colSums(TLBS)!=0))

            if(length(TLcand)>0){

                GapFlank<-getGapsFlanks(TLBS[,TLcand,drop=FALSE])
                ## returns $Gaps and $Flanks

                ## keep flank only if flank is not a best hit
                ## adjust tags otherwise

                for(i in 1:length(TLcand)){

                    flankCar<-unique(tree$car[myset][TLBS[,TLcand[i]]>0])
                    if(length(flankCar)!=1){
                        stop("Flanking CAR is not unique")
                    }
                    flankCarPos<-match(flankCar,mycars)

                    ## check whether flanks are best hits
                    flBH<-scafcarbest[s,flankCarPos]>=bestHitOpt

                    if(isTRUE(flBH)){ ## to modify in TLBS, if any
                        toMod<-c(toMod,TLcand[i])
                    }
                }
            }

            ## potentially adjust tags in TLL$TLbetween
            ##  (set best hits to zero and tag fragments of CAR on other
            ##   scaffolds with 1.0 instead of 0.5 [division by 2 below])
            ## >>>> Note that this will remove tags for multiple
            ##      scaffolds if the CAR has a best hit on more than
            ##      one scaffold. Although this might appear incorrect,
            ##      it is required to avoid that huge parts of some
            ##      scaffolds are tagged; see note in 'getBestHits()'
            if(length(toMod)>0){
                for(i in 1:length(toMod)){
                    TLL$TLbetween[myset,toMod[i]]<-0
                    tmp<-which(TLL$TLbetween[,toMod[i]]>0)
                    TLL$TLbetween[tmp,toMod[i]]<-2
                }
            }
        }

        ## obtain breakpoints for TLBcar
        bptB<-getBreakpntsSE(TLL$TLbetween[myset,,drop=FALSE]/2)
        ## returns $bptS and $bptE
        ## adjust column dimensions and append data to SYNT
        SYNT$TLBSbS<-appendData(SYNT$TLBSbS,bptB$bptS)
        SYNT$TLBSbE<-appendData(SYNT$TLBSbE,bptB$bptE)
    }

    SYNT$TLBS<-TLL$TLbetween/2


    rownames(SYNT$TLWS)<-rownames(TLL$TLwithin)
    rownames(SYNT$TLWSbS)<-rownames(TLL$TLwithin)
    rownames(SYNT$TLWSbE)<-rownames(TLL$TLwithin)
    rownames(SYNT$TLBSbS)<-rownames(TLL$TLwithin)
    rownames(SYNT$TLBSbE)<-rownames(TLL$TLwithin)

    return(SYNT)
}

## ------------------------------------------------------------------


## identify best hits for scaffold - CAR pairs (set to 1)
getBestHits<-function(markers,tree,myscafs,mycars){
    ## get marker numbers for each scaff - CAR pair
    scafcarasso<-matrix(0,nrow=length(myscafs),ncol=length(mycars))
    rownames(scafcarasso)<-myscafs
    colnames(scafcarasso)<-mycars

    for(s in 1:length(myscafs)){
        scaffpos<-which(markers$scaff==myscafs[s])
        tmpcars<-sort(unique(tree$car[scaffpos]))
        if(length(tmpcars)>0){
            for(i in tmpcars){
                scafcarasso[s,mycars==i]<-sum(tree$car[scaffpos]==i)
            }
        }
    }

    scafmax<-apply(scafcarasso,1,max)
    carmax<-apply(scafcarasso,2,max)

    ## identify best car per scaffold (can be ties)
    bestcar<-matrix(0,nrow=length(myscafs),ncol=length(mycars))
    rownames(bestcar)<-myscafs
    colnames(bestcar)<-mycars

    for(s in 1:length(myscafs)){
        tmpcars<-which(scafcarasso[s,]==scafmax[s])
        bestcar[s,tmpcars]<-1
    }

    ## identify best scaffold per car (can be ties)
    bestscaf<-matrix(0,nrow=length(myscafs),ncol=length(mycars))
    rownames(bestscaf)<-myscafs
    colnames(bestscaf)<-mycars

    for(i in 1:length(mycars)){
        tmpscafs<-which(scafcarasso[,i]==carmax[i])
        bestscaf[tmpscafs,i]<-1
    }

    ## get reciprocal best hits (there can still be ties)
    scafcarbest<-bestcar*bestscaf
    ## solve ties
    scafties<-which(rowSums(scafcarbest)>1)
    carties<-which(colSums(scafcarbest)>1)

    ## scafties
    ## decide based on relative size to next largest scaf for involved cars
    if(length(scafties)>0){
        for(s in scafties){
            tmpcars<-which(scafcarasso[s,]==scafmax[s])
            poc<-scafcarasso[s,tmpcars]/
                (apply(scafcarasso[,tmpcars],2,function(x)
                    sort(x,decreasing=TRUE)[2]))
            keep<-which(poc==max(poc)) ## can still be ties
            ## set cars with the smaller ratio to zero
            scafcarbest[s,tmpcars[-keep]]<-0
        }
    }
    ## carties
    ## decide based on relative size to next largest car for involved scafs
    if(length(carties)>0){
        for(i in carties){
            tmpscafs<-which(scafcarasso[,i]==carmax[i])
            pos<-scafcarasso[tmpscafs,i]/
                (apply(scafcarasso[tmpscafs,],1,function(x)
                    sort(x,decreasing=TRUE)[2]))
            keep<-which(pos==max(pos)) ## can still be ties
            ## set scafs with the smaller ratio to zero
            scafcarbest[tmpscafs[-keep],i]<-0
            ## it could be possible that scaf to keep has been set to zero
            ##  already above (not sure if it can ever happen), check
            if(sum(scafcarbest[,i])==0){
                scafcarbest[tmpscafs[keep],i]<-1
            }
        }
    }
    ## check if still ties, and if yes, set to zero
    ##  (admitting that there is no easy solution)
    scafties<-which(rowSums(scafcarbest)>1)
    carties<-which(colSums(scafcarbest)>1)
    if(length(scafties)>0){
        for(s in scafties){
            scafcarbest[s,]<-0
        }
    }
    if(length(carties)>0){
        for(i in carties){
            scafcarbest[,i]<-0
        }
    }

    if(sum(rowSums(scafcarbest)>1)>0 | sum(colSums(scafcarbest)>1)>0){
        stop("Problem while identifying mutual best hits")
    }

    ## secondary best hits:
    ## (1) no mutual best hit, but car is the largest on a scaffold
    ##     (-> 'bestcar'; there can be ties but that's fine)
    ## (2) no mutual best hit, but car has most of its markers on a scaffold
    ##     (-> 'bestscaf'; there can be ties and need to be solved)
    ## remove mutual best hits
    for(s in 1:length(myscafs)){
        bestscaf[s,scafcarbest[s,]>0]<-0
    }
    ## check for ties
    carties<-which(colSums(bestscaf)>1)
    ## if ties, set to zero (admitting that there is no easy solution)
    if(length(carties)>0){
        for(i in carties){
            bestscaf[,i]<-0
        }
    }

    ## combine matrices, giving weights
    scafcarbestS<-pmax(scafcarbest*3,bestscaf*2,bestcar)

    ## ## make sure that no car has more than one scaffold assigned
    ## ## (it's okay if a scaffold has multiple cars assigned)
    ## ## >>>> Note that this can potentially lead to tagging of
    ## ##      a huge part of a scaffold as TLBS if a car merges
    ## ##      two scaffolds but it just happened that there were
    ## ##      small translocations of other cars at the scaffold ends
    ## ##      so that they cannot be joined; has to be addressed
    ## ##      together with the tests in 'filterCars3()'
    ## for(i in 1:length(mycars)){
    ##     if(sum(scafcarbestS[,i]>0)>1){
    ##         tmp<-which(scafcarbestS[,i]==max(scafcarbestS[,i]))
    ##         if(length(tmp)>1){
    ##             ## set to zero (admitting that there is no easy solution)
    ##             scafcarbestS[,i]<-0
    ##         }else{
    ##             scafcarbestS[-tmp,i]<-0
    ##         }
    ##     }
    ## }

    return(scafcarbestS)
}

## ------------------------------------------------------------------


## make blocks for Q-nodes
##  Arguments: block start, block end, block elements, leafelem,
##             iter (required because function calls itself),
##             tomask (pre-determined elements that will not
##               be incorporated in a block)
makeBlocks<-function(blocksL,blstart,blend,blelem,leafelem,
                     iter,tomask=logical(0),OneIter=FALSE){

    tmprank<-myRank(blelem)

    ## identify duplicates
    dupel<-unique(blelem[duplicated(blelem)])

    ## if no pre-determined masks, set tomask vector to FALSE
    if(length(tomask)==0){
        tomask<-rep(FALSE,length(blelem))
    }

    ## -----------------
    ## initialize first block
    ## start, end, start rank element, end rank element,
    ##  count of leaf markers, order a-/descending, block rank
    blocks<-matrix(c(blstart[1],blend[1],rep(tmprank[1],2),
                     sum(leafelem[(blstart[1]):(blend[1])]),9,NA),
                   nrow=1,ncol=7)
    ## making blocks in presence of duplicated elements requires
    ##  additional if-else, as duplicated elements can result
    ##  in conflicts among block content and wrong ranks
    ## additionally, if pre-determined masks exist, those elements
    ##  should always form their own single-element block
    if(length(blelem)>1){
        for(i in 2:length(blelem)){
            ## potentially update current block if true
            if(abs(diff(tmprank[(i-1):i]))==1 &
               (blocks[nrow(blocks),6]==9 |
                    blocks[nrow(blocks),6]==diff(tmprank[(i-1):i]))){
                ## -> note: testing for abs(diff)==1 can result in
                ##    duplicated elements being added to different blocks
                ##    (but not the same block because of test for
                ##    direction); duplicated elements will never be adjacent
                if((is.element(blelem[i-1],dupel) &
                        (blocks[nrow(blocks),6]!=9 |
                             is.element(blelem[i],dupel))) |
                   tomask[i-1]==TRUE | tomask[i]==TRUE){
                    ## make new block if previous element was duplicated,
                    ##  but wasn't the first element of the block or
                    ##  the current element is also duplicated
                    ## -> duplicated elements are only allowed at either
                    ##    end of block, and are not allowed to be
                    ##    consecutive on the same block
                    ## additionally, if pre-determined masks exist,
                    ##  those elements will form their own block

                    ## as for iter == 1
                    ##  (this distinction is probably not necessary)
                    blocks<-rbind(blocks,c(blstart[i],blend[i],rep(tmprank[i],2),sum(leafelem[(blstart[i]):(blend[i])]),9,NA))
                }else{
                    ## update current block
                    blocks[nrow(blocks),2]<-blend[i]
                    blocks[nrow(blocks),4]<-tmprank[i]
                    blocks[nrow(blocks),5]<-blocks[nrow(blocks),5]+
                        sum(leafelem[(blstart[i]):(blend[i])])
                    blocks[nrow(blocks),6]<-diff(tmprank[(i-1):i]) ## +/- 1
                }
            }else{ ## make new block
                blocks<-rbind(blocks,c(blstart[i],blend[i],rep(tmprank[i],2),
                                       sum(leafelem[(blstart[i]):(blend[i])]),9,NA))
            }
        }
    }
    ## make block rank (this should be save when no duplicated
    ##  elements exist, but need adjustment otherwise - should
    ##  be covered now)
    blocks[,7]<-myRank(rowMeans(blocks[,3:4,drop=FALSE]))

    if(iter>1){
        if(nrow(blocksL[[iter-1]])==nrow(blocks)){
            ## no further simplification
            return(blocksL)
        }
    }
    blocksL[[iter]]<-blocks

    if(nrow(blocks)==1){
        ## no further simplification
        return(blocksL)
    }

    if(isTRUE(OneIter)){
        return(blocksL)
    }


    ## transfer pre-determined masks to next iteration
    if(nrow(blocks)<length(tomask)){
        nexttomask<-rep(FALSE,nrow(blocks))
        if(sum(tomask)>0){
            tmp<-which(tomask==TRUE)
            tomaskElem<-blstart[tmp]
            if(sum(tomaskElem!=blend[tmp])>0){
                stop("Masked blocks cannot contain several elements")
            }
            for(i in 1:length(tomaskElem)){
                k<-which(blocks[,1]==tomaskElem[i] & blocks[,2]==tomaskElem[i])
                if(length(k)!=1){
                    stop("Failed to transfer masks to next block iteration")
                }
                nexttomask[k]<-TRUE
            }
        }
        tomask<-nexttomask
    }


    makeBlocks(blocksL,blocksL[[iter]][,1],blocksL[[iter]][,2],
               blocksL[[iter]][,7],leafelem,iter+1,tomask)


}

## ------------------------------------------------------------------


## mask duplicated elements for each level of blocksL;
##  alternatively, use pre-determined masks
setMasks<-function(blocksL,allelem,dupli,tomask=logical(0)){

    if(length(tomask)==0){
        mask<-blocksL[[1]][,6]==9 &
            is.element(allelem[blocksL[[1]][,1]],dupli)
    }else{
        nexttomask<-rep(FALSE,nrow(blocksL[[1]]))
        if(sum(tomask)>0){
            tmp<-which(tomask==TRUE)
            for(i in tmp){
                k<-which(blocksL[[1]][,1]==i & blocksL[[1]][,2]==i)
                if(length(k)!=1){
                    stop("Failed to transfer pre-masks to masks")
                }
                nexttomask[k]<-TRUE
            }
        }
        mask<-nexttomask
    }

    maskL<-vector("list",length(blocksL))
    maskL[[1]]<-mask
    if(length(maskL)>1){
        ## this looks whether masked blocks still remain as
        ##  separate blocks when combining blocks
        for(z in 2:length(maskL)){
            maskL[[z]]<-rep(FALSE,nrow(blocksL[[z]]))
        }
        if(sum(mask)>0){
            ## determine whether any masked elements are still
            ##  not bound into larger blocks
            for(w in which(mask==TRUE)){
                maskelempos<-blocksL[[1]][w,1:2]
                for(z in 2:length(maskL)){
                    tmp<-which(blocksL[[z]][,1]==maskelempos[1] &
                                   blocksL[[z]][,2]==maskelempos[2])
                    if(length(tmp)>0){
                        maskL[[z]][tmp]<-TRUE
                    }
                }
            }
        }
    }
    return(maskL)
}

## ------------------------------------------------------------------


## assign orientation to Q-nodes
##  +1 ascending; -1 descending; 9 no direction
assignOri3<-function(blocksL,maskL,allelem,elemrows,
                     leafelem,leaves,testorientation){

    tmpblockori<-NA

    if(nrow(blocksL[[length(blocksL)]])==1){
        ## final combination of blocks exist

        ## tmpblockori=9 can happen if single elements
        ##  forms one final block;
        ##  this has been incorporated in downstream analysis

        ## test if level below top has just
        ##  two elements and if yes, potentiall modify orientation
        ##  dependent on orientation of these elements
        ##  (excluding masked elements)
        tmpblockori<-assignOriTwoElem(blocksL,maskL,elemrows,
                                      leaves,testorientation)
        ## should be 1, -1, or 9 with existence of final block

        if(is.na(tmpblockori)){
            stop("Something went wrong with assignment of block orientation")
        }

    }else if(nrow(blocksL[[length(blocksL)]])>1){
        ## final combination of blocks wasn't possible,
        ##  (can happen with four or more blocks, e.g. 2.4.1.3)
        ##  decide for order of last simplification based on largest
        ##  block(s) unless first and last element conform to order

        z<-length(blocksL)
        ## but first, try if removal of masked elements is sufficient
        ## (>>>> hasn't been exhaustively tested)
        unmaskedbl<-which(maskL[[z]]==FALSE)
        if(length(unmaskedbl)>0 & length(unmaskedbl)<length(maskL[[z]])){
            ## some, but not all are masked

            ## same two-element tests included as above;
            ##  make new blocks even if there is only
            ##  a single unmasked block for simplicity
            ## get elements within unmasked blocks
            unmaskedelem<-integer()
            for(k in unmaskedbl){
                ## get columns in elemrows
                unmaskedelem<-c(unmaskedelem,
                                (blocksL[[z]][k,1]):(blocksL[[z]][k,2]))
            }
            ## make blocks
            unmaskedblocksL<-vector("list",0)
            unmaskedblocksL<-makeBlocks(unmaskedblocksL,
                                        1:length(unmaskedelem),
                                        1:length(unmaskedelem),
                                        allelem[unmaskedelem],
                                        leafelem[unmaskedelem],1)
            if(nrow(unmaskedblocksL[[length(unmaskedblocksL)]])==1){
                ## final simplification exists
                tmpelem<-allelem[unmaskedelem]
                tmpdupli<-unique(tmpelem[duplicated(tmpelem)])
                unmaskedmaskL<-setMasks(unmaskedblocksL,tmpelem,tmpdupli)
                ## get correct subset for leaves and testorientation
                ## expand block over elements to markers
                unmaskedmarkers<-integer()
                for(k in unmaskedbl){
                    ## get columns in elemrows
                    tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
                    ## get position of markers
                    for(j in tmp){
                        unmaskedmarkers<-c(unmaskedmarkers,
                                           (elemrows[1,j]):(elemrows[2,j]))
                    }
                }
                ## >>>> get correct from-to positions for elements
                ##      requires below change to unmaskedelemrows
                ##      from elemrows[,unmaskedelem,drop=FALSE]
                unmaskedelemrows<-matrix(1,nrow=2,ncol=length(unmaskedelem))
                cntr<-0
                for(i in 1:length(unmaskedelem)){
                    cntr<-cntr+1
                    unmaskedelemrows[1,i]<-cntr
                    cntr<-cntr+(elemrows[2,unmaskedelem[i]]-elemrows[1,unmaskedelem[i]])
                    unmaskedelemrows[2,i]<-cntr
                }

                tmpblockori<-assignOriTwoElem(unmaskedblocksL,unmaskedmaskL,
                                              unmaskedelemrows,
                                              leaves[unmaskedmarkers],
                                              testorientation[unmaskedmarkers])
            }
        }

        ## removal of masked elements was sufficient
        if(!is.na(tmpblockori)){
            return(tmpblockori)
        }

        ## else, check if first and last (unmasked) element conform
        ##  to either ascending or descending order
        if(length(unmaskedbl)>1){
            tmpelem<-blocksL[[z]][unmaskedbl,7]
            if(tmpelem[1]==min(tmpelem) &
               tail(tmpelem,n=1L)==max(tmpelem)){
                tmpblockori<-1
                return(tmpblockori)
            }else if(tmpelem[1]==max(tmpelem) &
                     tail(tmpelem,n=1L)==min(tmpelem)){
                tmpblockori<- -1
                return(tmpblockori)
            }
        }

        ## else, only use largest block(s) that are not masked
        ## (use number of markers in block)
        blsize<-integer(nrow(blocksL[[z]]))
        ## expand block over elements to markers
        for(k in 1:nrow(blocksL[[z]])){
            ## get columns in elemrows
            tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
            ## get number of markers
            for(j in tmp){
                blsize[k]<-blsize[k]+diff(elemrows[,j])+1
            }
        }
        blnelem<-blocksL[[z]][,2]-blocksL[[z]][,1]+1

        ## masked elements can negatively affect order determination and
        ##  should be excluded, if possible
        ## exclude blocks that have many elements but formed by just one node
        candbl<-intersect(which(blnelem>1),unmaskedbl)
        if(length(candbl)>0){
            maxblsize<-max(blsize[candbl])
            ## use max #markers when at least 3 (or 2) elements in block
            largebl<-which(blnelem>2 & blsize==maxblsize & maskL[[z]]==FALSE)
            if(length(largebl)==0){
                largebl<-which(blnelem>1 & blsize==maxblsize & maskL[[z]]==FALSE)
            }
        }else{
            ## don't exclude masked elements
            candbl<-which(blnelem>1)
            if(length(candbl)>0){
                maxblsize<-max(blsize[candbl])
                largebl<-which(blnelem>2 & blsize==maxblsize)
                if(length(largebl)==0){
                    largebl<-which(blnelem>1 & blsize==maxblsize)
                }
            }else{
                maxblsize<-max(blsize)
                largebl<-which(blsize==maxblsize)
            }
        }

        ## same two-element tests included as above;
        ##  make new blocks even if there is only
        ##  a single largest block for simplicity
        ## get elements within largest block(s)
        largeelem<-integer()
        for(k in largebl){
            ## get columns in elemrows
            largeelem<-c(largeelem,
                            (blocksL[[z]][k,1]):(blocksL[[z]][k,2]))
        }
        largeblocksL<-vector("list",0)
        largeblocksL<-makeBlocks(largeblocksL,
                                 1:length(largeelem),
                                 1:length(largeelem),
                                 allelem[largeelem],
                                 leafelem[largeelem],1)
        if(nrow(largeblocksL[[length(largeblocksL)]])>1){
            ## still no final simplification, give up here
            ##  and assign order based on first and last element
            if(largeblocksL[[length(largeblocksL)]][1,7]<=
               largeblocksL[[length(largeblocksL)]][length(largebl),7]){
                tmpblockori<-1
            }else{
                tmpblockori<- -1
            }
        }else{
            ## final simplification exists
            tmpelem<-allelem[largeelem]
            tmpdupli<-unique(tmpelem[duplicated(tmpelem)])
            largemaskL<-setMasks(largeblocksL,tmpelem,tmpdupli)
            ## get correct subset for leaves and testorientation
            ## expand block over elements to markers
            largemarkers<-integer()
            for(k in largebl){
                ## get columns in elemrows
                tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
                ## get position of markers
                for(j in tmp){
                    largemarkers<-c(largemarkers,
                                       (elemrows[1,j]):(elemrows[2,j]))
                }
            }
            ## >>>> get correct from-to positions for elements
            ##      requires below change to largeelemrows
            ##      from elemrows[,largeelem,drop=FALSE]
            largeelemrows<-matrix(1,nrow=2,ncol=length(largeelem))
            cntr<-0
            for(i in 1:length(largeelem)){
                cntr<-cntr+1
                largeelemrows[1,i]<-cntr
                cntr<-cntr+(elemrows[2,largeelem[i]]-elemrows[1,largeelem[i]])
                largeelemrows[2,i]<-cntr
            }

            tmpblockori<-assignOriTwoElem(largeblocksL,largemaskL,
                                          largeelemrows,
                                          leaves[largemarkers],
                                          testorientation[largemarkers],
                                          useMask=FALSE)
        }


        if(tmpblockori==9){
            ## give up here and assign order based on first and last element
            if(blocksL[[z]][1,7]<=blocksL[[z]][nrow(blocksL[[z]]),7]){
                tmpblockori<-1
            }else{
                tmpblockori<- -1
            }

        }

        if(is.na(tmpblockori) | tmpblockori==9){
            stop("Something went wrong with assignment of block orientation")
        }

    }else{
        stop("There should be at least one block")
    }

    return(tmpblockori)

}

## ------------------------------------------------------------------


## modified order decision in case second-last block contains
##  only two elements (or there is only one level with one block);
##  orientation might be modified dependent on orientation of these
##  two elements, and also dependent on number of markers per element
assignOriTwoElem<-function(blocksL,maskL,elemrows,leaves,
                           testorientation,useMask=TRUE){

    if(useMask==FALSE){
        mymaskL<-maskL
        for(k in 1:length(mymaskL)){
            mymaskL[[k]]<-rep(FALSE,length(mymaskL[[k]]))
        }
    }else if(useMask==TRUE){
        mymaskL<-maskL
    }else{
        stop("useMask has to be TRUE or FALSE")
    }

    ## adjust elemrows to start at 1 and to have no gaps
    ##  (should not be the case if function is called correctly)
    myelemrows<-elemrows
    ## se<-1
    ## for(k in 1:ncol(myelemrows)){
    ##     ee<-se+diff(myelemrows[,k])
    ##     myelemrows[,k]<-c(se,ee)
    ##     se<-ee+1
    ## }
    if(sum(elemrows[2,]-elemrows[1,])+ncol(elemrows)!=length(leaves) |
       length(leaves)!=length(testorientation) |
       max(elemrows)!=length(leaves) | min(elemrows)!=1){
        stop("Indices in elemrows do not match number of elements")
    }


    modblockori<-NA

    bllev<-length(blocksL)
    nfinalblocks<-nrow(blocksL[[bllev]][!mymaskL[[bllev]],,drop=FALSE])

    ## test if highest multi-element block level, excluding masks,
    ##  has exactly two elements
    while(nfinalblocks==1 & bllev>1){
        bllev<-bllev-1
        nfinalblocks<-nrow(blocksL[[bllev]][!mymaskL[[bllev]],,drop=FALSE])
        ## store first encountered orientation of final block(s)
        ##  (potentially overwritten below)
        if(is.na(modblockori) | modblockori==9){
            modblockori<-blocksL[[bllev+1]][!mymaskL[[bllev+1]],,drop=FALSE][1,6]
        }
    }

    if(nfinalblocks==2){

        ## assign order of the two final blocks (with masked elements,
        ##  single final block might exist below highest bllev)
        modblockori<-blocksL[[bllev+1]][!mymaskL[[bllev+1]],,drop=FALSE][1,6]

        idx1<-which(!mymaskL[[bllev]])[1]
        idx2<-which(!mymaskL[[bllev]])[2]
        res1<-findOri(blocksL,bllev,idx1,myelemrows,testorientation,leaves)
        res2<-findOri(blocksL,bllev,idx2,myelemrows,testorientation,leaves)
        ## returns $ord (orientation) and $bllev (where ori!=9 was found)

        if(!is.na(modblockori) & modblockori!=9){
            ## ## count number of elements in each block
            ##e1<-(blocksL[[bllev]][idx1,1]):(blocksL[[bllev]][idx1,2])
            ##ne1<-0
            ##for(j in e1){
            ##    ne1<-ne1+myelemrows[2,j]-myelemrows[1,j]+1
            ##}
            ##e2<-(blocksL[[bllev]][idx2,1]):(blocksL[[bllev]][idx2,2])
            ##ne2<-0
            ##for(j in e2){
            ##    ne2<-ne2+myelemrows[2,j]-myelemrows[1,j]+1
            ##}
            ## potentially switch order
            ## (also switch if smaller block has no orientation)
            ## could be made more stringent with 2*ne1 >= ne2
            if(modblockori==1){
                if(res1$ord== -1 & res2$ord== -1){
                    modblockori<- -1
                }##else if((res1$ord== -1 & ne1 > ne2) |
                ##        (res2$ord== -1 & ne2 > ne1)){
                ##   ##modblockori<- -1
                ##   modblockori<-1 ## (don't switch)
                ##}
            }else if(modblockori== -1){
                if(res1$ord==1 & res2$ord==1){
                    modblockori<-1
                }##else if((res1$ord==1 & ne1 > ne2) |
                ##        (res2$ord==1 & ne2 > ne1)){
                ##   ##modblockori<-1
                ##   modblockori<- -1 ## (don't switch)
                ##}
            }
        }
    }else if(nfinalblocks==1){ ## bllev==1 because of while loop above

        ## assign order of single final block (with masked elements,
        ##  single final block might exist below highest bllev)
        modblockori<-blocksL[[bllev]][!mymaskL[[bllev]],,drop=FALSE][1,6]

        ## check number of elements of single final block and if two,
        ##  determine orientation of these elements if they are leaves

        ## get correct subset of leaves and testorientation
        tmp<-(blocksL[[bllev]][!mymaskL[[bllev]],1]):
            (blocksL[[bllev]][!mymaskL[[bllev]],2])
        mysub<-integer()
        for(j in tmp){
            mysub<-c(mysub,myelemrows[1,j]:myelemrows[2,j])
        }
        myleaves<-leaves[mysub]
        mytestorientation<-testorientation[mysub]

        if(length(myleaves)==2 & sum(myleaves)==2){
            ## (with only two leaves, none can be larger)
            if(!is.na(modblockori) & modblockori!=9){
                ord1<-mytestorientation[1]
                ord2<-mytestorientation[2]
                ## potentially switch order
                if(modblockori==1 & !is.na(ord1) & ord1== -1 &
                   !is.na(ord2) & ord2== -1){
                    modblockori<- -1
                }else if(modblockori== -1 & !is.na(ord1) & ord1==1 &
                         !is.na(ord2) & ord2==1){
                    modblockori<-1
                }
            }
        }
    }

    return(modblockori)
    ## can be NA, 1, -1, or 9 (for a single-element block)
}

## ------------------------------------------------------------------


## check ascending order for rearrangements
checkAdjAscend<-function(blocks,mask,nelem,usePreMask=FALSE,
                         returnAdj=FALSE){

    adjA<-matrix(0,nrow=nrow(blocks),ncol=2)

    tpElA<-matrix(NA,nrow=nelem,ncol=0)

    ## tag wrong adjacencies (with masked duplicated elements)
    ##  (with duplicated elements masked, block ranks
    ##   should never be equal for different blocks)
    if(sum(!mask)>1){
        blockrank<-myRank(rowMeans(blocks[!mask,3:4]))
        if(usePreMask==FALSE & length(blockrank)!=length(unique(blockrank))){
            stop("Ranking of blocks is ambiguous")
        }
        blockid<-(1:nrow(blocks))[!mask]
        ## tests ends
        if(blockrank[1]!=min(blockrank)){
            adjA[blockid[1],1]<-1
        }
        if(tail(blockrank,n=1L)!=max(blockrank)){
            adjA[tail(blockid,n=1L),2]<-1
        }
        ## test rest
        if(usePreMask==FALSE){
            for(i in 2:length(blockrank)){
                if(blockrank[i]!=blockrank[i-1]+1){
                    adjA[blockid[i-1],2]<-1
                    adjA[blockid[i],1]<-1
                }
            }
        }else{ ## identical elements might be successive
            for(i in 2:length(blockrank)){
                if(blockrank[i]!=blockrank[i-1]+1 &
                   blockrank[i]!=blockrank[i-1]){
                    adjA[blockid[i-1],2]<-1
                    adjA[blockid[i],1]<-1
                }
            }
        }
        if(sum(adjA)>1){
            ## adjust for correct ends (-> excess 1's)
            tmpstart<-which(adjA[,1]==1)
            tmpend<-which(adjA[,2]==1)
            if(tmpend[1]<tmpstart[1]){
                adjA[tmpend[1],2]<-0
                tmpend<-tmpend[-1] ## remove foregoing end
            }
            if(length(tmpstart)>length(tmpend)){
                adjA[tail(tmpstart,n=1L),1]<-0 ## remove trailing start
                tmpstart<-tmpstart[-length(tmpstart)]
            }
            ## check that no end comes before start
            ##  and no previous end comes before next start
            if(sum(tmpend-tmpstart<0)>0){
                stop("Start adjacencies incorrect")
            }
            if(sum(c(0,tmpend)-c(tmpstart,nrow(adjA)+1)>=0)>0){
                stop("End adjacencies incorrect")
            }
        }
    }
    ## tag wrong adjacencies for masked elements only
    if(usePreMask==FALSE & sum(mask)>0){
        blockrank<-blocks[,7]
        blockid<-(1:nrow(blocks))[mask]
        for(i in 1:length(blockid)){
            if(i==1 & blockid[i]==1){
                if(blockrank[blockid[i]]!=min(blockrank)){
                    adjA[blockid[i],1]<-1
                }
                ## test following block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]+1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]+1]-1){
                    adjA[blockid[i],2]<-1
                }else{
                    adjA[blockid[i]+1,1]<-0
                }
            }else if(i==length(blockid) &
                     blockid[i]==nrow(blocks)){
                if(blockrank[blockid[i]]!=max(blockrank)){
                    adjA[blockid[i],2]<-1
                }
                ## test preceding block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]-1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]-1]+1){
                    adjA[blockid[i],1]<-1
                }else{
                    adjA[blockid[i]-1,2]<-0
                }
            }else{
                ## test preceding block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]-1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]-1]+1){
                    adjA[blockid[i],1]<-1
                }else{
                    adjA[blockid[i]-1,2]<-0
                }
                ## test following block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]+1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]+1]-1){
                    adjA[blockid[i],2]<-1
                }else{
                    adjA[blockid[i]+1,1]<-0
                }
            }
        }
        if(sum(adjA)>1){
            ## adjust for correct ends (-> excess 1's)
            tmpstart<-which(adjA[,1]==1)
            tmpend<-which(adjA[,2]==1)
            if(tmpend[1]<tmpstart[1]){
                adjA[tmpend[1],2]<-0
                tmpend<-tmpend[-1] ## remove foregoing end
            }
            if(length(tmpstart)>length(tmpend)){
                adjA[tail(tmpstart,n=1L),1]<-0 ## remove trailing start
                tmpstart<-tmpstart[-length(tmpstart)]
            }
            ## adjust that no end comes before start
            ##  and no previous end comes before next start
            tmpstart<-c(tmpstart,nrow(adjA)+1)
            laststart<-0
            lastend<-0
            while(length(tmpstart)>0 & length(tmpend)>0){
                if(laststart>=nrow(adjA) | lastend>=nrow(adjA)){
                    break
                }
                if(tmpend[1]-tmpstart[1]<0){
                    ## tag start
                    laststart<-lastend+1
                    tmpstart<-c(laststart,tmpstart)
                    adjA[laststart,1]<-1
                    next
                }else if(lastend-tmpstart[1]>=0){
                    ## tag end
                    tmpend<-c(lastend,tmpend)
                    lastend<-tmpstart[1]-1
                    adjA[lastend,2]<-1
                    next
                }else{
                    lastend<-tmpend[1]
                    tmpend<-tmpend[-1]
                    laststart<-tmpstart[1]
                    tmpstart<-tmpstart[-1]
                }
            }
        }
    }else if(usePreMask==TRUE & sum(mask)>0){
        ## tag wrong adjacencies for pre-masked elements only
        blockid<-(1:nrow(blocks))[mask]
        for(i in 1:length(blockid)){
            adjA[blockid[i],1]<-1
            adjA[blockid[i],2]<-1
        }
    }
    if(returnAdj==FALSE){
        ## go through adjacencies and tag elements in-between
        ## >>> this might make more tags than necessary, but also avoid
        ##     potentially complicated decision which block has been
        ##     transposed as they will often be multiple options
        ## >>> could be interesting to store info about the number of
        ##     breaks in ordering though (i.e., # columns in TLWC)
        if(sum(adjA)>1){
            tmpstart<-which(adjA[,1]==1)
            tmpend<-which(adjA[,2]==1)
            if(length(tmpstart)>length(tmpend)){
                stop("Adjacency boundaries incorrect")
            }
            if(sum(tmpend-tmpstart<0)>0){
                stop("Adjacency boundaries incorrect")
            }
            for(i in 1:length(tmpstart)){
                tpEl<-rep(0,nelem)
                ## get positions of elements
                tmp1<-tmpstart[i]
                tmp<-(blocks[tmp1,1]):(blocks[tmp1,2])
                while(tmp1<tmpend[i]){
                    tmp1<-tmp1+1
                    tmp<-c(tmp,(blocks[tmp1,1]):(blocks[tmp1,2]))
                }
                tpEl[tmp]<-rep(1,length(tmp))
                tpElA<-cbind(tpElA,tpEl)
            }
        }
        return(tpElA)
    }else{
        return(adjA)
    }
}

## ------------------------------------------------------------------


## check descending order for rearrangements
checkAdjDescend<-function(blocks,mask,nelem,usePreMask=FALSE,
                          returnAdj=FALSE){

    adjD<-matrix(0,nrow=nrow(blocks),ncol=2)

    tpElD<-matrix(NA,nrow=nelem,ncol=0)

    ## tag wrong adjacencies (with masked duplicated elements)
    ##  (with duplicated elements masked, block ranks
    ##   should never be equal for different blocks)
    if(sum(!mask)>1){
        blockrank<-myRank(rowMeans(blocks[!mask,3:4]))
        if(usePreMask==FALSE & length(blockrank)!=length(unique(blockrank))){
            stop("Ranking of blocks is ambiguous")
        }
        blockid<-(1:nrow(blocks))[!mask]
        ## tests ends
        if(blockrank[1]!=max(blockrank)){
            adjD[blockid[1],1]<-1
        }
        if(tail(blockrank,n=1L)!=min(blockrank)){
            adjD[tail(blockid,n=1L),2]<-1
        }
        ## test rest
        if(usePreMask==FALSE){
            for(i in 2:length(blockrank)){
                if(blockrank[i]!=blockrank[i-1]-1){
                    adjD[blockid[i-1],2]<-1
                    adjD[blockid[i],1]<-1
                }
            }
        }else{ ## identical elements might be successive
            for(i in 2:length(blockrank)){
                if(blockrank[i]!=blockrank[i-1]-1 &
                   blockrank[i]!=blockrank[i-1]){
                    adjD[blockid[i-1],2]<-1
                    adjD[blockid[i],1]<-1
                }
            }
        }
        if(sum(adjD)>1){
            ## adjust for correct ends (-> excess 1's)
            tmpstart<-which(adjD[,1]==1)
            tmpend<-which(adjD[,2]==1)
            if(tmpend[1]<tmpstart[1]){
                adjD[tmpend[1],2]<-0
                tmpend<-tmpend[-1] ## remove foregoing end
            }
            if(length(tmpstart)>length(tmpend)){
                adjD[tail(tmpstart,n=1L),1]<-0 ## remove trailing start
                tmpstart<-tmpstart[-length(tmpstart)]
            }
            ## check that no end comes before start
            ##  and no previous end comes before next start
            if(sum(tmpend-tmpstart<0)>0){
                stop("Start adjacencies incorrect")
            }
            if(sum(c(0,tmpend)-c(tmpstart,nrow(adjD)+1)>=0)>0){
                stop("End adjacencies incorrect")
            }
        }
    }
    ## tag wrong adjacencies for duplicated elements only
    if(usePreMask==FALSE & sum(mask)>0){
        blockrank<-blocks[,7]
        blockid<-(1:nrow(blocks))[mask]
        for(i in 1:length(blockid)){
            if(i==1 & blockid[i]==1){
                if(blockrank[blockid[i]]!=max(blockrank)){
                    adjD[blockid[i],1]<-1
                }
                ## test following block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]+1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]+1]+1){
                    adjD[blockid[i],2]<-1
                }else{
                    adjD[blockid[i]+1,1]<-0
                }
            }else if(i==length(blockid) &
                     blockid[i]==nrow(blocks)){
                if(blockrank[blockid[i]]!=min(blockrank)){
                    adjD[blockid[i],2]<-1
                }
                ## test preceding block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]-1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]-1]-1){
                    adjD[blockid[i],1]<-1
                }else{
                    adjD[blockid[i]-1,2]<-0
                }
            }else{
                ## test preceding block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]-1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]-1]-1){
                    adjD[blockid[i],1]<-1
                }else{
                    adjD[blockid[i]-1,2]<-0
                }
                ## test following block
                if(blockrank[blockid[i]]!=blockrank[blockid[i]+1] &
                   blockrank[blockid[i]]!=blockrank[blockid[i]+1]+1){
                    adjD[blockid[i],2]<-1
                }else{
                    adjD[blockid[i]+1,1]<-0
                }
            }
        }
        if(sum(adjD)>1){
            ## adjust for correct ends (-> excess 1's)
            tmpstart<-which(adjD[,1]==1)
            tmpend<-which(adjD[,2]==1)
            if(tmpend[1]<tmpstart[1]){
                adjD[tmpend[1],2]<-0
                tmpend<-tmpend[-1] ## remove foregoing end
            }
            if(length(tmpstart)>length(tmpend)){
                adjD[tail(tmpstart,n=1L),1]<-0 ## remove trailing start
                tmpstart<-tmpstart[-length(tmpstart)]
            }
            ## adjust that no end comes before start
            ##  and no previous end comes before next start
            tmpstart<-c(tmpstart,nrow(adjD)+1)
            laststart<-0
            lastend<-0
            while(length(tmpstart)>0 & length(tmpend)>0){
                if(laststart>=nrow(adjD) | lastend>=nrow(adjD)){
                    break
                }
                if(tmpend[1]-tmpstart[1]<0){
                    ## tag start
                    laststart<-lastend+1
                    tmpstart<-c(laststart,tmpstart)
                    adjD[laststart,1]<-1
                    next
                }else if(lastend-tmpstart[1]>=0){
                    ## tag end
                    tmpend<-c(lastend,tmpend)
                    lastend<-tmpstart[1]-1
                    adjD[lastend,2]<-1
                    next
                }else{
                    lastend<-tmpend[1]
                    tmpend<-tmpend[-1]
                    laststart<-tmpstart[1]
                    tmpstart<-tmpstart[-1]
                }
            }
        }
    }else if(usePreMask==TRUE & sum(mask)>0){
        ## tag wrong adjacencies for pre-masked elements only
        blockid<-(1:nrow(blocks))[mask]
        for(i in 1:length(blockid)){
            adjD[blockid[i],1]<-1
            adjD[blockid[i],2]<-1
        }
    }
    if(returnAdj==FALSE){
        ## go through adjacencies and tag elements in-between
        ## >>> this might make more tags than necessary, but also avoid
        ##     potentially complicated decision which block has been
        ##     transposed as they will often be multiple options
        ## >>> could be interesting to store info about the number of
        ##     breaks in ordering though (i.e., # columns in TLWC)
        if(sum(adjD)>1){
            tmpstart<-which(adjD[,1]==1)
            tmpend<-which(adjD[,2]==1)
            if(length(tmpstart)>length(tmpend)){
                stop("Adjacency boundaries incorrect")
            }
            if(sum(tmpend-tmpstart<0)>0){
                stop("Adjacency boundaries incorrect")
            }
            for(i in 1:length(tmpstart)){
                tpEl<-rep(0,nelem)
                ## get positions of elements
                tmp1<-tmpstart[i]
                tmp<-(blocks[tmp1,1]):(blocks[tmp1,2])
                while(tmp1<tmpend[i]){
                    tmp1<-tmp1+1
                    tmp<-c(tmp,(blocks[tmp1,1]):(blocks[tmp1,2]))
                }
                tpEl[tmp]<-rep(1,length(tmp))
                tpElD<-cbind(tpElD,tpEl)
            }
        }
        return(tpElD)
    }else{
        return(adjD)
    }
}

## ------------------------------------------------------------------


## check orientation of blocks, for positive higher-order orientation
checkOriAscend<-function(blocksL,blocklev,allelem,dupli,
                         elemrows,leaves,testorientation,
                         blockoriL,subbl,splitblockid,
                         remWgt=0.05){

    tpElA<-matrix(NA,ncol=0,nrow=length(allelem))
    ivElA<-matrix(NA,ncol=0,nrow=length(allelem))

    ## some blocks are descending
    if(sum(blockoriL[[blocklev]][subbl]== -1)>0){
        for(k in subbl){
            inv<-9
            if(blockoriL[[blocklev]][k]== -1){
                ## expand block for further testing
                ## get elemrows
                tmp<-(blocksL[[blocklev]][k,1]):(blocksL[[blocklev]][k,2])
                ## if lower-level blocks exist, get corresponding rows
                if(blocklev>1){
                    tmp1<-which(blocksL[[blocklev-1]][,1]>=
                                    blocksL[[blocklev]][k,1] &
                                        blocksL[[blocklev-1]][,2]<=
                                            blocksL[[blocklev]][k,2])
                    if(length(tmp1)>2){
                        ## at least three elements
                        ## -> inversion
                        inv<-1
                    }else if(length(tmp1)==2){
                        ## two elements -> further tests
                        e1<-blocksL[[blocklev-1]][tmp1[1],1]:
                            blocksL[[blocklev-1]][tmp1[1],2]
                        e2<-blocksL[[blocklev-1]][tmp1[2],1]:
                            blocksL[[blocklev-1]][tmp1[2],2]

                        if((length(e1)==1 &
                                sum(is.element(allelem[e1],dupli))>0) |
                           (length(e2)==1 &
                                sum(is.element(allelem[e2],dupli))>0)){
                            ## block contains duplicated element
                            ## one level down the hierarchy
                            ## (not further bound to other elements)
                            ## -> don't call it inversion
                            ## >>>> it might be an inversion; this
                            ##      could potentially be tested later
                            inv<-0
                        }else if(sum(blocksL[[blocklev-1]][tmp1,5]!=0)==0){
                            ## all elements are nodes
                            ## >>>> might be a transposition too
                            inv<-1
                        }else{
                            ## at least one leaf:
                            ## -> take orientation into account

                            ## always transposition unless both
                            ##  elements have negative orientation;
                            ##  with blocklev==2, all elements are leaves
                            ##  and majority of markers for each element have
                            ##  negative orientation (all markers will need
                            ##  to have orientation info available)

                            if(blocklev==2){
                                ## get positions of markers
                                tmpe1<-numeric()
                                tmpe2<-numeric()
                                for(j in e1){
                                    tmpe1<-c(tmpe1,(elemrows[1,j]):(elemrows[2,j]))
                                }
                                for(j in e2){
                                    tmpe2<-c(tmpe2,(elemrows[1,j]):(elemrows[2,j]))
                                }
                                ## if all are leafmarkers:
                                if(sum(leaves[c(tmpe1,tmpe2)]==0)==0){
                                    ## orientation is negative for majority
                                    ##  of markers for each element
                                    if((sum(!is.na(testorientation[tmpe1]) & testorientation[tmpe1]== -1)>length(tmpe1)/2) & (sum(!is.na(testorientation[tmpe2]) & testorientation[tmpe2]== -1)>length(tmpe2)/2)){
                                        ## -> inversion
                                        inv<-1
                                    }else{
                                        ## -> transposition
                                        inv<-0
                                    }
                                }else{
                                    ## -> transposition
                                    inv<-0
                                }
                            }else{ ## blocklev>2
                                ## get orientation:
                                res1<-findOri(blocksL,blocklev-1,tmp1[1],
                                              elemrows,testorientation,leaves)
                                res2<-findOri(blocksL,blocklev-1,tmp1[2],
                                              elemrows,testorientation,leaves)
                                ## orientation is negative for both
                                ##  (possible if some blocks were excluded
                                ##   in 'checkAdjComplexRearr()')
                                if(res1$ord== -1 & res2$ord== -1){
                                    inv<-1
                                }else{
                                    inv<-0
                                }
                            }
                        } ## close tests with blocks containing leaves
                    } ## close testing for two elements
                }else{ ## blocklev == 1

                    ## tmp has elementrows stored

                    if(diff(blocksL[[blocklev]][k,1:2])+1>2){
                        ## at least three elements
                        ## -> inversion
                        inv<-1
                    }else if(diff(blocksL[[blocklev]][k,1:2])+1==2){
                        ## only two elements -> further tests
                        e1<-blocksL[[blocklev]][k,1]
                        e2<-blocksL[[blocklev]][k,2]

                        if(sum(is.element(allelem[tmp],dupli))>0){
                            ## block contains duplicated element
                            ## -> don't call it inversion
                            ## >>>> it might be an inversion; this
                            ##      could potentially be tested later
                            inv<-0
                        }else if(blocksL[[blocklev]][k,5]==0){
                            ## all elements are nodes
                            ## >>>> might be a transposition too
                            inv<-1
                        }else{
                            ## at least one leaf

                            ## always transposition unless all are leaves
                            ##  and all have negative orientation

                            ## get positions of markers
                            tmp2<-numeric()
                            for(j in tmp){
                                tmp2<-c(tmp2,(elemrows[1,j]):(elemrows[2,j]))
                            }
                            ## if all are leafmarkers, check that orientation
                            ##  is negative for all markers (| if no orientation
                            ##  info for leaves available, tag block
                            ##  as inversion)
                            if(sum(leaves[tmp2]==0)==0){
                                if((sum(!is.na(testorientation[tmp2]) & testorientation[tmp2]== -1)==length(tmp2)) | (sum(is.na(testorientation[tmp2]))==length(tmp2))){
                                    ## -> inversion
                                    inv<-1
                                }else{
                                    ## -> transposition
                                    inv<-0
                                }
                            }else{
                                ## -> transposition
                                inv<-0
                            }
                        }
                    }
                } ## close blocklev == 1


                if(inv==0){ ## add (additional) rearrangement
                    ## append column to tmptp for each
                    tp1<-rep(0,length(allelem))
                    tp2<-rep(0,length(allelem))
                    ## get number of markers for each part
                    ne1<-0
                    ne2<-0
                    for(j in e1){
                        ne1<-ne1+elemrows[2,j]-elemrows[1,j]+1
                    }
                    for(j in e2){
                        ne2<-ne2+elemrows[2,j]-elemrows[1,j]+1
                    }
                    if(ne1<ne2){ ## left part smaller
                        ## tag positions of affected block row
                        tp1[e1]<-rep(1-remWgt,length(e1))
                        tp2[e2]<-rep(remWgt,length(e2))
                    }else if(ne1>ne2){ ## right part smaller
                        ## tag positions of affected block row
                        tp1[e1]<-rep(remWgt,length(e1))
                        tp2[e2]<-rep(1-remWgt,length(e2))
                    }else{ ## tie
                        ## -> tag the part that has the wrong orientation,
                        ##      if any, otherwise give 0.5 tags to both
                        if(blocklev>1){
                            ## get ori of two parts from lower level
                            ori1<-blocksL[[blocklev-1]][tmp1[1],6]
                            ori2<-blocksL[[blocklev-1]][tmp1[2],6]
                            ## descend to lower level, if needed
                            y<-blocklev-2
                            if(y==0){
                                pos1<-tmp1[1]
                            }
                            while(ori1==9 & y>=1){
                                pos1<-which(blocksL[[y]][,1]==blocksL[[blocklev-1]][tmp1[1],1] & blocksL[[y]][,2]==blocksL[[blocklev-1]][tmp1[1],2])
                                if(length(pos1)!=1){
                                    stop("problem while descending block levels")
                                }
                                ori1<-blocksL[[y]][pos1,6]
                                y<-y-1
                            }
                            y<-blocklev-2
                            if(y==0){
                                pos2<-tmp1[2]
                            }
                            while(ori2==9 & y>=1){
                                pos2<-which(blocksL[[y]][,1]==blocksL[[blocklev-1]][tmp1[2],1] & blocksL[[y]][,2]==blocksL[[blocklev-1]][tmp1[2],2])
                                if(length(pos2)!=1){
                                    stop("problem while descending block levels")
                                }
                                ori2<-blocksL[[y]][pos2,6]
                                y<-y-1
                            }
                            ## if ori is still '9', check if it is a leaf
                            ##  (then get orientation)
                            if(ori1==9){
                                j<-blocksL[[1]][pos1,1]
                                ## (in fact, j should be equal e1 then)
                                tmpe1<-(elemrows[1,j]):(elemrows[2,j])
                                if(length(tmpe1)==1 &
                                   sum(leaves[tmpe1]==0)==0 &
                                   sum(is.na(testorientation[tmpe1]))==0){
                                    ## single marker and a leaf and
                                    ##  orientation info available
                                    ori1<-testorientation[tmpe1] ## 1 or -1
                                }
                            }
                            if(ori2==9){
                                ## get positions of markers
                                j<-blocksL[[1]][pos2,1]
                                ## (in fact, j should be equal e2 then)
                                tmpe2<-(elemrows[1,j]):(elemrows[2,j])
                                if(length(tmpe2)==1 &
                                   sum(leaves[tmpe2]==0)==0 &
                                   sum(is.na(testorientation[tmpe2]))==0){
                                    ## single marker and a leaf and
                                    ##  orientation info available
                                    ori2<-testorientation[tmpe2] ## 1 or -1
                                }
                            }

                        }else{ ## blocklev == 1
                            ## only two elements, e1 and e2
                            ## check if they are leaves
                            ##  (then get orientation)
                            tmpe1<-(elemrows[1,e1]):(elemrows[2,e1])
                            if(length(tmpe1)==1 &
                               sum(leaves[tmpe1]==0)==0 &
                               sum(is.na(testorientation[tmpe1]))==0){
                                ## single marker and a leaf and
                                ##  orientation info available
                                ori1<-testorientation[tmpe1] ## 1 or -1
                            }else{
                                ori1<-9
                            }
                            tmpe2<-(elemrows[1,e2]):(elemrows[2,e2])
                            if(length(tmpe2)==1 &
                               sum(leaves[tmpe2]==0)==0 &
                               sum(is.na(testorientation[tmpe2]))==0){
                                ## single marker and a leaf and
                                ##  orientation info available
                                ori2<-testorientation[tmpe2] ## 1 or -1
                            }else{
                                ori2<-9
                            }
                        }
                        if(ori1== -1 & ori2!= -1){
                            ## tag positions of affected block row
                            tp1[e1]<-rep(1-remWgt,length(e1))
                            tp2[e2]<-rep(remWgt,length(e2))
                        }else if(ori1!= -1 & ori2== -1){
                            ## tag positions of affected block row
                            tp1[e1]<-rep(remWgt,length(e1))
                            tp2[e2]<-rep(1-remWgt,length(e2))
                        }else{
                            ## tag positions of affected block row
                            tp1[e1]<-rep(0.5,length(e1))
                            tp2[e2]<-rep(0.5,length(e2))
                        }
                    }
                    ## make splitblockid (whether element was untagged or not)
                    splitblockid[e1]<-paste0(splitblockid[e1],
                                             rep(".1",length(e1)))
                    splitblockid[e2]<-paste0(splitblockid[e2],
                                             rep(".2",length(e2)))

                    ## bind tags together
                    tpElA<-cbind(tpElA,tp1,tp2)
                    ## -> also need to switch orientation of
                    ##    blockrow if classified as transposition
                    ##    to correctly determine inversions on
                    ##    lower level of hierarchy
                    blockoriL[[blocklev]][k]<-1

                }else if(inv==1){ ## inversion
                    ## append column to tmpiv
                    iv<-rep(0,length(allelem))
                    ## tag positions of affected block row
                    iv[tmp]<-rep(1,length(tmp))
                    ivElA<-cbind(ivElA,iv)
                }else{
                    stop("determining inversions failed")
                }

            } ## close test for blocks that are '-1'
        } ## close loop over block rows
    } ## close some blocks are descending

    return(list(tpElA=tpElA,ivElA=ivElA,blockoriL=blockoriL,
                splitblockid=splitblockid))

}

## ------------------------------------------------------------------


## check orientation of blocks, for negative higher-order orientation
checkOriDescend<-function(blocksL,blocklev,allelem,dupli,
                          elemrows,leaves,testorientation,
                          blockoriL,subbl,splitblockid,
                          remWgt=0.05){

    tpElD<-matrix(NA,ncol=0,nrow=length(allelem))
    ivElD<-matrix(NA,ncol=0,nrow=length(allelem))

    ## some blocks are ascending
    if(sum(blockoriL[[blocklev]][subbl]==1)>0){
        for(k in subbl){
            inv<-9
            if(blockoriL[[blocklev]][k]==1){
                ## expand block for further testing
                ## get elemrows
                tmp<-(blocksL[[blocklev]][k,1]):(blocksL[[blocklev]][k,2])
                ## if lower-level blocks exist, get corresponding rows
                if(blocklev>1){
                    tmp1<-which(blocksL[[blocklev-1]][,1]>=
                                    blocksL[[blocklev]][k,1] &
                                        blocksL[[blocklev-1]][,2]<=
                                            blocksL[[blocklev]][k,2])
                    if(length(tmp1)>2){
                        ## at least three elements
                        ## -> inversion
                        inv<-1
                    }else if(length(tmp1)==2){
                        ## two elements -> further tests
                        e1<-blocksL[[blocklev-1]][tmp1[1],1]:
                            blocksL[[blocklev-1]][tmp1[1],2]
                        e2<-blocksL[[blocklev-1]][tmp1[2],1]:
                            blocksL[[blocklev-1]][tmp1[2],2]

                        if((length(e1)==1 &
                                sum(is.element(allelem[e1],dupli))>0) |
                           (length(e2)==1 &
                                sum(is.element(allelem[e2],dupli))>0)){
                            ## block contains duplicated element
                            ## one level down the hierarchy
                            ## (not further bound to other elements)
                            ## -> don't call it inversion
                            ## >>>> it might be an inversion; this
                            ##      could potentially be tested later
                            inv<-0
                        }else if(sum(blocksL[[blocklev-1]][tmp1,5]!=0)==0){
                            ## all elements are nodes
                            ## >>>> might be a transposition too
                            inv<-1
                        }else{
                            ## at least one leaf:
                            ## -> take orientation into account

                            ## always transposition unless both
                            ##  elements have positive orientation;
                            ##  with blocklev==2, all elements are leaves
                            ##  and majority of markers for each element have
                            ##  positive orientation (all markers will need
                            ##  to have orientation info available)

                            if(blocklev==2){
                                ## get positions of markers
                                tmpe1<-numeric()
                                tmpe2<-numeric()
                                for(j in e1){
                                    tmpe1<-c(tmpe1,(elemrows[1,j]):(elemrows[2,j]))
                                }
                                for(j in e2){
                                    tmpe2<-c(tmpe2,(elemrows[1,j]):(elemrows[2,j]))
                                }
                                ## if all are leafmarkers:
                                if(sum(leaves[c(tmpe1,tmpe2)]==0)==0){
                                    ## orientation is positive for majority
                                    ##  of markers for each element
                                    if((sum(!is.na(testorientation[tmpe1]) & testorientation[tmpe1]==1)>length(tmpe1)/2) & (sum(!is.na(testorientation[tmpe2]) & testorientation[tmpe2]==1)>length(tmpe2)/2)){
                                        ## -> inversion
                                        inv<-1
                                    }else{
                                        ## -> transposition
                                        inv<-0
                                    }
                                }else{
                                    ## -> transposition
                                    inv<-0
                                }
                            }else{ ## blocklev>2
                                ## get orientation:
                                res1<-findOri(blocksL,blocklev-1,tmp1[1],
                                              elemrows,testorientation,leaves)
                                res2<-findOri(blocksL,blocklev-1,tmp1[2],
                                              elemrows,testorientation,leaves)
                                ## orientation is positive for both
                                ##  (possible if some blocks were excluded
                                ##   in 'checkAdjComplexRearr()')
                                if(res1$ord==1 & res2$ord==1){
                                    inv<-1
                                }else{
                                    inv<-0
                                }
                            }
                        } ## close tests with blocks containing leaves
                    } ## close testing for two elements
                }else{ ## blocklev == 1

                    ## tmp has elementrows stored

                    if(diff(blocksL[[blocklev]][k,1:2])+1>2){
                        ## at least three elements
                        ## -> inversion
                        inv<-1
                    }else if(diff(blocksL[[blocklev]][k,1:2])+1==2){
                        ## only two elements -> further tests
                        e1<-blocksL[[blocklev]][k,1]
                        e2<-blocksL[[blocklev]][k,2]

                        if(sum(is.element(allelem[tmp],dupli))>0){
                            ## block contains duplicated element
                            ## -> don't call it inversion
                            ## >>>> it might be an inversion; this
                            ##      could potentially be tested later
                            inv<-0
                        }else if(blocksL[[blocklev]][k,5]==0){
                            ## all elements are nodes
                            ## >>>> might be a transposition too
                            inv<-1
                        }else{
                            ## at least one leaf

                            ## always transposition unless all are leaves
                            ##  and all have positive orientation

                            ## get positions of markers
                            tmp2<-numeric()
                            for(j in tmp){
                                tmp2<-c(tmp2,(elemrows[1,j]):(elemrows[2,j]))
                            }
                            ## if all are leafmarkers, check that orientation
                            ##  is positive for all markers (| if no orientation
                            ##  info for leaves available, tag block
                            ##  as inversion)
                            if(sum(leaves[tmp2]==0)==0){
                                if((sum(!is.na(testorientation[tmp2]) & testorientation[tmp2]==1)==length(tmp2)) | (sum(is.na(testorientation[tmp2]))==length(tmp2))){
                                    ## -> inversion
                                    inv<-1
                                }else{
                                    ## -> transposition
                                    inv<-0
                                }
                            }else{
                                ## -> transposition
                                inv<-0
                            }
                        }
                    }
                } ## close blocklev == 1


                if(inv==0){ ## add (additional) rearrangement
                    ## append column to tmptp for each
                    tp1<-rep(0,length(allelem))
                    tp2<-rep(0,length(allelem))
                    ## get number of markers for each part
                    ne1<-0
                    ne2<-0
                    for(j in e1){
                        ne1<-ne1+elemrows[2,j]-elemrows[1,j]+1
                    }
                    for(j in e2){
                        ne2<-ne2+elemrows[2,j]-elemrows[1,j]+1
                    }
                    if(ne1<ne2){ ## left part smaller
                        ## tag positions of affected block row
                        tp1[e1]<-rep(1-remWgt,length(e1))
                        tp2[e2]<-rep(remWgt,length(e2))
                    }else if(ne1>ne2){ ## right part smaller
                        ## tag positions of affected block row
                        tp1[e1]<-rep(remWgt,length(e1))
                        tp2[e2]<-rep(1-remWgt,length(e2))
                    }else{ ## tie
                        ## -> tag the part that has the wrong orientation,
                        ##      if any, otherwise give 0.5 tags to both
                        if(blocklev>1){
                            ## get ori of two parts from lower level
                            ori1<-blocksL[[blocklev-1]][tmp1[1],6]
                            ori2<-blocksL[[blocklev-1]][tmp1[2],6]
                            ## descend to lower level, if needed
                            y<-blocklev-2
                            if(y==0){
                                pos1<-tmp1[1]
                            }
                            while(ori1==9 & y>=1){
                                pos1<-which(blocksL[[y]][,1]==blocksL[[blocklev-1]][tmp1[1],1] & blocksL[[y]][,2]==blocksL[[blocklev-1]][tmp1[1],2])
                                if(length(pos1)!=1){
                                    stop("problem while descending block levels")
                                }
                                ori1<-blocksL[[y]][pos1,6]
                                y<-y-1
                            }
                            y<-blocklev-2
                            if(y==0){
                                pos2<-tmp1[2]
                            }
                            while(ori2==9 & y>=1){
                                pos2<-which(blocksL[[y]][,1]==blocksL[[blocklev-1]][tmp1[2],1] & blocksL[[y]][,2]==blocksL[[blocklev-1]][tmp1[2],2])
                                if(length(pos2)!=1){
                                    stop("problem while descending block levels")
                                }
                                ori2<-blocksL[[y]][pos2,6]
                                y<-y-1
                            }
                            ## if ori is still '9', check if it is a leaf
                            ##  (then get orientation)
                            if(ori1==9){
                                j<-blocksL[[1]][pos1,1]
                                ## (in fact, j should be equal e1 then)
                                tmpe1<-(elemrows[1,j]):(elemrows[2,j])
                                if(length(tmpe1)==1 &
                                   sum(leaves[tmpe1]==0)==0 &
                                   sum(is.na(testorientation[tmpe1]))==0){
                                    ## single marker and a leaf and
                                    ##  orientation info available
                                    ori1<-testorientation[tmpe1] ## 1 or -1
                                }
                            }
                            if(ori2==9){
                                ## get positions of markers
                                j<-blocksL[[1]][pos2,1]
                                ## (in fact, j should be equal e2 then)
                                tmpe2<-(elemrows[1,j]):(elemrows[2,j])
                                if(length(tmpe2)==1 &
                                   sum(leaves[tmpe2]==0)==0 &
                                   sum(is.na(testorientation[tmpe2]))==0){
                                    ## single marker and a leaf and
                                    ##  orientation info available
                                    ori2<-testorientation[tmpe2] ## 1 or -1
                                }
                            }

                        }else{ ## blocklev == 1
                            ## only two elements, e1 and e2
                            ## check if they are leaves
                            ##  (then get orientation)
                            tmpe1<-(elemrows[1,e1]):(elemrows[2,e1])
                            if(length(tmpe1)==1 &
                               sum(leaves[tmpe1]==0)==0 &
                               sum(is.na(testorientation[tmpe1]))==0){
                                ## single marker and a leaf and
                                ##  orientation info available
                                ori1<-testorientation[tmpe1] ## 1 or -1
                            }else{
                                ori1<-9
                            }
                            tmpe2<-(elemrows[1,e2]):(elemrows[2,e2])
                            if(length(tmpe2)==1 &
                               sum(leaves[tmpe2]==0)==0 &
                               sum(is.na(testorientation[tmpe2]))==0){
                                ## single marker and a leaf and
                                ##  orientation info available
                                ori2<-testorientation[tmpe2] ## 1 or -1
                            }else{
                                ori2<-9
                            }
                        }
                        if(ori1==1 & ori2!=1){
                            ## tag positions of affected block row
                            tp1[e1]<-rep(1-remWgt,length(e1))
                            tp2[e2]<-rep(remWgt,length(e2))
                        }else if(ori1!=1 & ori2==1){
                            ## tag positions of affected block row
                            tp1[e1]<-rep(remWgt,length(e1))
                            tp2[e2]<-rep(1-remWgt,length(e2))
                        }else{
                            ## tag positions of affected block row
                            tp1[e1]<-rep(0.5,length(e1))
                            tp2[e2]<-rep(0.5,length(e2))
                        }
                    }
                    ## make splitblockid (whether element was untagged or not)
                    splitblockid[e1]<-paste0(splitblockid[e1],
                                             rep(".1",length(e1)))
                    splitblockid[e2]<-paste0(splitblockid[e2],
                                             rep(".2",length(e2)))

                    ## bind tags together
                    tpElD<-cbind(tpElD,tp1,tp2)
                    ## -> also need to switch orientation of
                    ##    blockrow if classified as transposition
                    ##    to correctly determine inversions on
                    ##    lower level of hierarchy
                    blockoriL[[blocklev]][k]<- -1

                }else if(inv==1){ ## inversion
                    ## append column to tmpiv
                    iv<-rep(0,length(allelem))
                    ## tag positions of affected block row
                    iv[tmp]<-rep(1,length(tmp))
                    ivElD<-cbind(ivElD,iv)
                }else{
                    stop("determining inversions failed")
                }

            } ## close test for blocks that are '+1'
        } ## close loop over block rows
    } ## close some blocks are descending

    return(list(tpElD=tpElD,ivElD=ivElD,blockoriL=blockoriL,
                splitblockid=splitblockid))

}

## ------------------------------------------------------------------


## for Q-nodes and in the presence of TPelem, need to determine
##  first which parts of TPelem are likely to be rearranged, given
##  arrangment further down the hierarchy
## >>>> the checks included here are a rough approximation
##      and likely do not cover all special cases
## >>>> in some cases both the flanks and inserts might receive
##      tags -> could for example be improved by taking relative
##      sizes into account, as done for CARs
## >>>> makePreMasks function currently unused, as not working
##      perfectly in all cases
makePreMasks<-function(tmptree,allelem,elemrows,TPelem,n,nhier){

    ## get vector of duplicated elements (TPelem=1)
    dupli<-numeric()
    ## in addition, bind flanks of same TP element together
    TPflanks<-matrix(0,ncol=0,nrow=length(allelem))

    if(nrow(TPelem)>0){
        dupli<-unique(allelem[apply(TPelem,2,function(x) is.element(1,x))])
        for(k in 1:length(dupli)){
            ## bind flanks of same element together
            tmpflank<-numeric(length(allelem))
            tmpflank[which(allelem==dupli[k])]<-1
            TPflanks<-cbind(TPflanks,tmpflank)
        }
    }

    rankelem<-myRank(allelem)

    adjL<-vector("list",length(dupli))
    adjA<-matrix(NA,nrow=length(allelem),ncol=2)
    adjD<-matrix(NA,nrow=length(allelem),ncol=2)
    maskA<-matrix(0,nrow=length(allelem),ncol=length(dupli))
    maskD<-matrix(0,nrow=length(allelem),ncol=length(dupli))


    ## for each dupli
    for(d in 1:length(dupli)){

        flankpos<-sort(which(TPflanks[,d]==1))


        ## for ascending order -------------

        ## check if flank adjacencies are correct to outsides and inserts
        adjAofi<-adjA
        adjAofi[flankpos,]<-0

        for(f in 1:length(flankpos)){
            if(flankpos[f]==1){ ## flank at front end
                if(rankelem[flankpos[f]]!=min(rankelem)){
                    adjAofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]+1]-1){
                    adjAofi[flankpos[f],2]<-1
                }
            }else if(flankpos[f]==length(rankelem)){ ## flank at rear end
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]-1]+1){
                    adjAofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=max(rankelem)){
                    adjAofi[flankpos[f],2]<-1
                }
            }else{ ## flank in middle
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]-1]+1){
                    adjAofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]+1]-1){
                    adjAofi[flankpos[f],2]<-1
                }
            }
        }

        ## check if adjacencies of outsides and inserts are correct
        ##  (excluding flanks, or taking them into account)
        adjAoi<-adjA

        for(f in 1:length(flankpos)){
            if(f==1 & flankpos[f]==1){ ## flank at front end
                if(rankelem[flankpos[f]+1]==min(rankelem) |
                   rankelem[flankpos[f]+1]==rankelem[flankpos[f]]+1){
                    adjAoi[flankpos[f]+1,1]<-0
                }else{
                    adjAoi[flankpos[f]+1,1]<-1
                }
            }else if(f==length(flankpos) &
                     flankpos[f]==length(rankelem)){ ## flank at rear end
                if(rankelem[flankpos[f]-1]==max(rankelem) |
                   rankelem[flankpos[f]-1]==rankelem[flankpos[f]]-1){
                    adjAoi[flankpos[f]-1,2]<-0
                }else{
                    adjAoi[flankpos[f]-1,2]<-1
                }
            }else{ ## flank in middle
                if(rankelem[flankpos[f]-1]==rankelem[flankpos[f]+1]){
                    ## need additional checking at lower-level node
                    adjAoi[flankpos[f]-1,2]<-9
                    adjAoi[flankpos[f]+1,1]<-9
                }else if(rankelem[flankpos[f]-1]==rankelem[flankpos[f]+1]-1 |
                         (rankelem[flankpos[f]-1]==rankelem[flankpos[f]]-1 &
                              rankelem[flankpos[f]+1]==rankelem[flankpos[f]]+1)){
                    adjAoi[flankpos[f]-1,2]<-0
                    adjAoi[flankpos[f]+1,1]<-0
                }else{
                    adjAoi[flankpos[f]-1,2]<-1
                    adjAoi[flankpos[f]+1,1]<-1
                }
            }
        }

        ## proceed to lower-level node
        ##  (>>> take into account that there might be no lower-level)
        ##  (>>> for nodes, order could be ascending or descending;
        ##       if elements are leaves, check orientation >>> not done!)

        ## check if adjacencies of outsides and inserts are correct
        ##  (excluding flanks) from unsolved cases above ('9') in 'adjAoi'
        ## >>> ignore elements that are NA already
        ## >>> this is identical to what's done for 'adjDoi'
        m<-n+1
        while(sum(!is.na(adjAoi) & adjAoi==9)>0 & m<=nhier){
            ## get positions of markers
            tmpstart<-which(!is.na(adjAoi[,2]) & adjAoi[,2]==9)
            tmpend<-which(!is.na(adjAoi[,1]) & adjAoi[,1]==9)
            if(length(tmpstart)!=length(tmpend)){
                stop("Testing unsolved adjacencies requires equal number of start and end positions")
            }
            for(f in 1:length(tmpstart)){
                tmp1<-(elemrows[1,tmpstart[f]]):(elemrows[2,tmpstart[f]])
                tmp2<-(elemrows[1,tmpend[f]]):(elemrows[2,tmpend[f]])

                oielem<-tmptree[c(tmp1,tmp2),2+(m-1)*2+1]
                ## all elements of outsides and inserts
                ## determine elements that are NA already
                ##  (NA elements would not have been tagged with a 9)
                toexcl<-sort(unique(which(is.na(oielem) | oielem=="NA")))
                oielem<-myRank(oielem) ## puts NAs last
                ## determine position of boundary within oielem
                oilast<-1+diff(elemrows[,tmpstart[f]])

                ## if Q-node, check order of elements
                mynode<-unique(tmptree[c(tail(tmp1,n=1L),head(tmp2,n=1L)),
                                       2+(m-1)*2])
                if(length(mynode)!=1 | !is.element(mynode,c("Q","P"))){
                    stop("Incorrect node assignment")
                }
                if(mynode=="Q"){
                    ## check if last element in tmp1 has no conflict
                    ##  with first element in tmp2
                    if(oielem[oilast]==oielem[oilast+1]-1 |
                       oielem[oilast]==oielem[oilast+1]+1){
                        ## also still need to check for duplicated elements
                        adjAoi[tmpstart[f],2]<-0
                        adjAoi[tmpend[f],1]<-0
                    }else if(oielem[oilast]!=oielem[oilast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node
                        adjAoi[tmpstart[f],2]<-1
                        adjAoi[tmpend[f],1]<-1
                    }
                }else{ ## P-node
                    if(oielem[oilast]!=oielem[oilast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node;
                        ## also still need to check for duplicated elements
                        adjAoi[tmpstart[f],2]<-0
                        adjAoi[tmpend[f],1]<-0
                    }
                }
                ## check for duplicated elements
                ##  (might overwrite a '0' from above)
                ## make vector with all preceding hierarchies
                tmphiercol<-seq((2+n*2)-1,(2+(m-1)*2)-1,2)
                tmpprev<-apply(as.matrix(tmptree[c(tmp1,tmp2),tmphiercol]),1,
                               function(x) paste0(x,collapse="-"))
                oielem.c<-paste(tmpprev,oielem,sep="-")
                if(length(toexcl)>0){
                    unioi.c<-unique(oielem.c[-toexcl])
                }else{
                    unioi.c<-unique(oielem.c)
                }
                for(i in unioi.c){
                    tmp3<-which(oielem.c==unioi.c[i])
                    if(length(tmp3)>1 &
                       sum(diff(tmp3)!=1)>0){ ## there are duplicated elements
                        ## get their boundaries
                        oibnd<-numeric(nrow(tmptree))
                        oibnd[c(tmp1,tmp2)[tmp3]]<-1

                        ## test whether they align with flank breakpoints
                        ##  but same dupli doesn't cross flank breakpoints
                        ##  (only consider breakpoints, not genes in-between)
                        if((oibnd[elemrows[2,tmpstart[f]]]==1 &
                                oibnd[elemrows[1,tmpend[f]]]==0)|
                           (oibnd[elemrows[1,tmpend[f]]]==1 &
                                oibnd[elemrows[2,tmpstart[f]]]==0)){
                            adjAoi[tmpstart[f],2]<-1
                            adjAoi[tmpend[f],1]<-1
                        }
                    }
                }
            }
            m<-m+1
        }



        ## for descending order -------------

        ## check if flank adjacencies are correct to outsides and inserts
        adjDofi<-adjD
        adjDofi[flankpos,]<-0

        for(f in 1:length(flankpos)){
            if(flankpos[f]==1){ ## flank at front end
                if(rankelem[flankpos[f]]!=max(rankelem)){
                    adjDofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]+1]+1){
                    adjDofi[flankpos[f],2]<-1
                }
            }else if(flankpos[f]==length(rankelem)){ ## flank at rear end
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]-1]-1){
                    adjDofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=min(rankelem)){
                    adjDofi[flankpos[f],2]<-1
                }
            }else{ ## flank in middle
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]-1]-1){
                    adjDofi[flankpos[f],1]<-1
                }
                if(rankelem[flankpos[f]]!=rankelem[flankpos[f]+1]+1){
                    adjDofi[flankpos[f],2]<-1
                }
            }
        }

        ## check if adjacencies of outsides and inserts are correct
        ##  (excluding flanks, or taking them into account)
        adjDoi<-adjD

        for(f in 1:length(flankpos)){
            if(f==1 & flankpos[f]==1){ ## flank at front end
                if(rankelem[flankpos[f]+1]==max(rankelem) |
                   rankelem[flankpos[f]+1]==rankelem[flankpos[f]]-1){
                    adjDoi[flankpos[f]+1,1]<-0
                }else{
                    adjDoi[flankpos[f]+1,1]<-1
                }
            }else if(f==length(flankpos) &
                     flankpos[f]==length(rankelem)){ ## flank at rear end
                if(rankelem[flankpos[f]-1]==min(rankelem) |
                   rankelem[flankpos[f]-1]==rankelem[flankpos[f]]+1){
                    adjDoi[flankpos[f]-1,2]<-0
                }else{
                    adjDoi[flankpos[f]-1,2]<-1
                }
            }else{ ## flank in middle
                if(rankelem[flankpos[f]-1]==rankelem[flankpos[f]+1]){
                    ## need additional checking at lower-level node
                    adjDoi[flankpos[f]-1,2]<-9
                    adjDoi[flankpos[f]+1,1]<-9
                }else if(rankelem[flankpos[f]-1]==rankelem[flankpos[f]+1]+1 |
                         (rankelem[flankpos[f]-1]==rankelem[flankpos[f]]+1 &
                              rankelem[flankpos[f]+1]==rankelem[flankpos[f]]-1)){
                    adjDoi[flankpos[f]-1,2]<-0
                    adjDoi[flankpos[f]+1,1]<-0
                }else{
                    adjDoi[flankpos[f]-1,2]<-1
                    adjDoi[flankpos[f]+1,1]<-1
                }
            }
        }

        ## proceed to lower-level node
        ##  (>>> take into account that there might be no lower-level)
        ##  (>>> for nodes, order could be ascending or descending;
        ##       if elements are leaves, check orientation >>> not done!)

        ## check if adjacencies of outsides and inserts are correct
        ##  (excluding flanks) from unsolved cases above ('9') in 'adjDoi'
        ## >>> ignore elements that are NA already
        ## >>> this is identical to what's done for 'adjAoi'
        m<-n+1
        while(sum(!is.na(adjDoi) & adjDoi==9)>0 & m<=nhier){
            ## get positions of markers
            tmpstart<-which(!is.na(adjDoi[,2]) & adjDoi[,2]==9)
            tmpend<-which(!is.na(adjDoi[,1]) & adjDoi[,1]==9)
            if(length(tmpstart)!=length(tmpend)){
                stop("Testing unsolved adjacencies requires equal number of start and end positions")
            }
            for(f in 1:length(tmpstart)){
                tmp1<-(elemrows[1,tmpstart[f]]):(elemrows[2,tmpstart[f]])
                tmp2<-(elemrows[1,tmpend[f]]):(elemrows[2,tmpend[f]])

                oielem<-tmptree[c(tmp1,tmp2),2+(m-1)*2+1]
                ## all elements of outsides and inserts
                ## determine elements that are NA already
                ##  (NA elements would not have been tagged with a 9)
                toexcl<-sort(unique(which(is.na(oielem) | oielem=="NA")))
                oielem<-myRank(oielem) ## puts NAs last
                ## determine position of boundary within oielem
                oilast<-1+diff(elemrows[,tmpstart[f]])

                ## if Q-node, check order of elements
                mynode<-unique(tmptree[c(tail(tmp1,n=1L),head(tmp2,n=1L)),
                                       2+(m-1)*2])
                if(length(mynode)!=1 | !is.element(mynode,c("Q","P"))){
                    stop("Incorrect node assignment")
                }
                if(mynode=="Q"){
                    ## check if last element in tmp1 has no conflict
                    ##  with first element in tmp2
                    if(oielem[oilast]==oielem[oilast+1]-1 |
                       oielem[oilast]==oielem[oilast+1]+1){
                        ## also still need to check for duplicated elements
                        adjDoi[tmpstart[f],2]<-0
                        adjDoi[tmpend[f],1]<-0
                    }else if(oielem[oilast]!=oielem[oilast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node
                        adjDoi[tmpstart[f],2]<-1
                        adjDoi[tmpend[f],1]<-1
                    }
                }else{ ## P-node
                    if(oielem[oilast]!=oielem[oilast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node;
                        ## also still need to check for duplicated elements
                        adjDoi[tmpstart[f],2]<-0
                        adjDoi[tmpend[f],1]<-0
                    }
                }
                ## check for duplicated elements
                ##  (might overwrite a '0' from above)
                ## make vector with all preceding hierarchies
                tmphiercol<-seq((2+n*2)-1,(2+(m-1)*2)-1,2)
                tmpprev<-apply(as.matrix(tmptree[c(tmp1,tmp2),tmphiercol]),1,
                               function(x) paste0(x,collapse="-"))
                oielem.c<-paste(tmpprev,oielem,sep="-")
                if(length(toexcl)>0){
                    unioi.c<-unique(oielem.c[-toexcl])
                }else{
                    unioi.c<-unique(oielem.c)
                }
                for(i in unioi.c){
                    tmp3<-which(oielem.c==unioi.c[i])
                    if(length(tmp3)>1 &
                       sum(diff(tmp3)!=1)>0){ ## there are duplicated elements
                        ## get their boundaries
                        oibnd<-numeric(nrow(tmptree))
                        oibnd[c(tmp1,tmp2)[tmp3]]<-1

                        ## test whether they align with flank breakpoints
                        ##  but same dupli doesn't cross flank breakpoints
                        ##  (only consider breakpoints, not genes in-between)
                        if((oibnd[elemrows[2,tmpstart[f]]]==1 &
                                oibnd[elemrows[1,tmpend[f]]]==0)|
                           (oibnd[elemrows[1,tmpend[f]]]==1 &
                                oibnd[elemrows[2,tmpstart[f]]]==0)){
                            adjDoi[tmpstart[f],2]<-1
                            adjDoi[tmpend[f],1]<-1
                        }
                    }
                }
            }
            m<-m+1
        }



        ## for either order -------------

        ## proceed to lower-level node
        ##  (>>> take into account that there might be no lower-level)
        ##  (>>> for nodes, order could be ascending or descending;
        ##       if elements are leaves, check orientation >>> not done!)

        ## check if flank adjacencies are correct when excluding inserts
        adjf1<-adjA
        adjf1[flankpos,]<-0

        ## get positions of markers
        tmp<-numeric()
        for(f in flankpos){
            tmp<-c(tmp,(elemrows[1,f]):(elemrows[2,f]))
        }
        ## check for duplicated elements
        flankelem<-tmptree[tmp,3+n*2] ## all elements of flank
        flankelem<-myRank(flankelem)
        unifl<-unique(flankelem)
        for(i in unifl){
            tmp2<-which(flankelem==unifl[i])
            if(length(tmp2)>1 &
               sum(diff(tmp2)!=1)>0){ ## there are duplicated elements
                ## get their boundaries
                flbnd<-numeric(nrow(tmptree))
                flbnd[tmp[tmp2]]<-1
                ## test whether they align with flank breakpoints
                ##  but same dupli doesn't cross flank breakpoints
                ##  (only consider breakpoints, not genes in-between)
                for(f in 1:length(flankpos)){
                    if(flankpos[f]==1){ ## start of tree
                        if(flbnd[1]==1){
                            adjf1[1,1]<-1
                        }
                        if(flbnd[elemrows[2,1]]==1 &
                           flbnd[elemrows[1,flankpos[f+1]]]==0){
                            adjf1[1,2]<-1
                        }
                    }else if(flankpos[f]==length(allelem)){ ## end of tree
                        if(flbnd[elemrows[2,flankpos[f-1]]]==0 &
                           flbnd[elemrows[1,flankpos[f]]]==1){
                            adjf1[flankpos[f],1]<-1
                        }
                        if(flbnd[length(flbnd)]==1){
                            adjf1[flankpos[f],2]<-1
                        }
                    }else{ ## middle
                        if(f==1){
                            if(flbnd[elemrows[1,flankpos[f]]]==1){
                                adjf1[flankpos[f],1]<-1
                            }
                            if(flbnd[elemrows[2,flankpos[f]]]==1 &
                               flbnd[elemrows[1,flankpos[f+1]]]==0){
                                adjf1[flankpos[f],2]<-1
                            }
                        }else if(f==length(flankpos)){
                            if(flbnd[elemrows[2,flankpos[f-1]]]==0 &
                               flbnd[elemrows[1,flankpos[f]]]==1){
                                adjf1[flankpos[f],1]<-1
                            }
                            if(flbnd[elemrows[2,flankpos[f]]]==1){
                                adjf1[flankpos[f],2]<-1
                            }
                        }else{
                            if(flbnd[elemrows[2,flankpos[f-1]]]==0 &
                               flbnd[elemrows[1,flankpos[f]]]==1){
                                adjf1[flankpos[f],1]<-1
                            }
                            if(flbnd[elemrows[2,flankpos[f]]]==1 &
                               flbnd[elemrows[1,flankpos[f+1]]]==0){
                                adjf1[flankpos[f],2]<-1
                            }
                        }
                    }
                }
            }
        }

        ## if Q-node, check order of flank elements
        ##  (allow both ascending and descending order)
        adjf2<-adjA

        ## determine start/end positions of flanks within flankelem
        flstart<-1
        flend<-diff(elemrows[,flankpos[1]])+1
        for(f in 2:length(flankpos)){
            flstart<-c(flstart,flend[length(flend)]+1)
            flend<-c(flend,diff(elemrows[,flankpos[f]])+
                         flstart[length(flstart)])
        }

        mynode<-unique(tmptree[tmp,2+n*2])
        if(length(mynode)!=1 | !is.element(mynode,c("Q","P"))){
            stop("Incorrect node assignment")
        }

        if(mynode=="Q"){
            for(f in 2:length(flankpos)){
                if(flankelem[flend[f-1]]==flankelem[flstart[f]]){
                    ## need additional checking at lower-level node
                    adjf2[flankpos[f-1],2]<-9
                    adjf2[flankpos[f],1]<-9
                }else if(flankelem[flend[f-1]]==flankelem[flstart[f]]-1 |
                         flankelem[flend[f-1]]==flankelem[flstart[f]]+1){
                    adjf2[flankpos[f-1],2]<-0
                    adjf2[flankpos[f],1]<-0
                }else{
                    adjf2[flankpos[f-1],2]<-1
                    adjf2[flankpos[f],1]<-1
                }
            }
        }else{ ## P-node
            for(f in 2:length(flankpos)){
                if(flankelem[flend[f-1]]==flankelem[flstart[f]]){
                    ## need additional checking at lower-level node
                    adjf2[flankpos[f-1],2]<-9
                    adjf2[flankpos[f],1]<-9
                }else{
                    adjf2[flankpos[f-1],2]<-0
                    adjf2[flankpos[f],1]<-0
                }
            }
        }


        ## check if flank adjacencies are correct when excluding inserts
        ## from unsolved cases above ('9') in 'adjf2'
        ## >>> ignore elements that are NA already
        m<-n+2
        while(sum(!is.na(adjf2) & adjf2==9)>0 & m<=nhier){
            ## get positions of markers
            tmpstart<-which(!is.na(adjf2[,2]) & adjf2[,2]==9)
            tmpend<-which(!is.na(adjf2[,1]) & adjf2[,1]==9)
            if(length(tmpstart)!=length(tmpend)){
                stop("Testing unsolved adjacencies requires equal number of start and end positions")
            }
            for(f in 1:length(tmpstart)){
                tmp1<-(elemrows[1,tmpstart[f]]):(elemrows[2,tmpstart[f]])
                tmp2<-(elemrows[1,tmpend[f]]):(elemrows[2,tmpend[f]])

                flelem<-tmptree[c(tmp1,tmp2),2+(m-1)*2+1]
                ## all elements of the two flanks under consideration
                ## determine elements that are NA already
                ##  (NA elements would not have been tagged with a 9)
                toexcl<-sort(unique(which(is.na(flelem) | flelem=="NA")))
                flelem<-myRank(flelem) ## puts NAs last
                ## determine position of boundary within flelem
                fllast<-1+diff(elemrows[,tmpstart[f]])

                ## if Q-node, check order of elements
                mynode<-unique(tmptree[c(tail(tmp1,n=1L),head(tmp2,n=1L)),
                                       2+(m-1)*2])
                if(length(mynode)!=1 | !is.element(mynode,c("Q","P"))){
                    stop("Incorrect node assignment")
                }
                if(mynode=="Q"){
                    ## check if last element in tmp1 has no conflict
                    ##  with first element in tmp2
                    if(flelem[fllast]==flelem[fllast+1]-1 |
                       flelem[fllast]==flelem[fllast+1]+1){
                        ## also still need to check for duplicated elements
                        adjf2[tmpstart[f],2]<-0
                        adjf2[tmpend[f],1]<-0
                    }else if(flelem[fllast]!=flelem[fllast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node
                        adjf2[tmpstart[f],2]<-1
                        adjf2[tmpend[f],1]<-1
                    }
                }else{ ## P-node
                    if(flelem[fllast]!=flelem[fllast+1]){
                        ## otherwise, need additional checking at
                        ##  next lower-level node;
                        ## also still need to check for duplicated elements
                        adjf2[tmpstart[f],2]<-0
                        adjf2[tmpend[f],1]<-0
                    }
                }
                ## check for duplicated elements
                ##  (might overwrite a '0' from above)
                ## make vector with all preceding hierarchies
                tmphiercol<-seq((2+(n+1)*2)-1,(2+(m-1)*2)-1,2)
                tmpprev<-apply(as.matrix(tmptree[c(tmp1,tmp2),tmphiercol]),1,
                               function(x) paste0(x,collapse="-"))
                flelem.c<-paste(tmpprev,flelem,sep="-")
                if(length(toexcl)>0){
                    unifl.c<-unique(flelem.c[-toexcl])
                }else{
                    unifl.c<-unique(flelem.c)
                }
                for(i in unifl.c){
                    tmp3<-which(flelem.c==unifl.c[i])
                    if(length(tmp3)>1 &
                       sum(diff(tmp3)!=1)>0){ ## there are duplicated elements
                        ## get their boundaries
                        flbnd<-numeric(nrow(tmptree))
                        flbnd[c(tmp1,tmp2)[tmp3]]<-1

                        ## test whether they align with flank breakpoints
                        ##  but same dupli doesn't cross flank breakpoints
                        ##  (only consider breakpoints, not genes in-between)
                        if((flbnd[elemrows[2,tmpstart[f]]]==1 &
                                flbnd[elemrows[1,tmpend[f]]]==0)|
                           (flbnd[elemrows[1,tmpend[f]]]==1 &
                                flbnd[elemrows[2,tmpstart[f]]]==0)){
                            adjf2[tmpstart[f],2]<-1
                            adjf2[tmpend[f],1]<-1
                        }
                    }
                }
            }
            m<-m+1
        }

        adjL[[d]]<-list(adjAofi=adjAofi,adjDofi=adjDofi,
                        adjAoi=adjAoi,adjDoi=adjDoi,
                        adjf1=adjf1,adjf2=adjf2)
    }






    ## for each dupli, decide whether to mask flanks or inserts

    for(d in 1:length(dupli)){

        flankpos<-sort(which(TPflanks[,d]==1))


        ## for ascending order -------------

        fifA<-numeric(2*length(flankpos)-1)

        for(f in 1:length(flankpos)){
            if(f==1 & flankpos[f]==1){ ## flank at front end
                if(adjL[[d]]$adjAofi[flankpos[f],1]==1){
                    fifA[1]<-1
                }
                if(adjL[[d]]$adjAofi[flankpos[f],2]==1){
                    if(adjL[[d]]$adjAoi[flankpos[f]+1,1]==0){
                        fifA[1]<-1
                    }else{
                        fifA[2]<-1
                    }
                }
                if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                   adjL[[d]]$adjf2[flankpos[f],2]==0){
                    fifA[2]<-1
                }
            }else if(f==length(flankpos) &
                     flankpos[f]==length(rankelem)){ ## flank at rear end
                if(adjL[[d]]$adjAofi[flankpos[f],2]==1){
                    fifA[length(fifA)]<-1
                }
                if(adjL[[d]]$adjAofi[flankpos[f],1]==1){
                    if(adjL[[d]]$adjAoi[flankpos[f]-1,2]==0){
                        fifA[length(fifA)]<-1
                    }else{
                        fifA[length(fifA)-1]<-1
                    }
                }
                if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                   adjL[[d]]$adjf2[flankpos[f],1]==0){
                    fifA[length(fifA)-1]<-1
                }
            }else{ ## flank in middle
                if(f==1){
                    if(adjL[[d]]$adjAofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]-1,2]==0){
                        fifA[1]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]+1,1]==0){
                        fifA[1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==1){
                        fifA[1]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjf2[flankpos[f],2]==1){
                        fifA[1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                       adjL[[d]]$adjf2[flankpos[f],2]==0){
                        fifA[2]<-1
                    }
                }else if(f==length(flankpos)){
                    if(adjL[[d]]$adjAofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]+1,1]==0){
                        fifA[length(fifA)]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]-1,2]==0){
                        fifA[length(fifA)]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==1){
                        fifA[length(fifA)]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjf2[flankpos[f],1]==1){
                        fifA[length(fifA)]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                       adjL[[d]]$adjf2[flankpos[f],1]==0){
                        fifA[length(fifA)-1]<-1
                    }
                }else{
                    if(adjL[[d]]$adjAofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]-1,2]==0){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjAoi[flankpos[f]+1,1]==0){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==1){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==1){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjf2[flankpos[f],1]==1){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjAofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjf2[flankpos[f],2]==1){
                        fifA[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                       adjL[[d]]$adjf2[flankpos[f],1]==0){
                        fifA[f*2-2]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                       adjL[[d]]$adjf2[flankpos[f],2]==0){
                        fifA[f*2]<-1
                    }
                }
            }
        }


        ## for descending order -------------

        fifD<-numeric(2*length(flankpos)-1)

        for(f in 1:length(flankpos)){
            if(f==1 & flankpos[f]==1){ ## flank at front end
                if(adjL[[d]]$adjDofi[flankpos[f],1]==1){
                    fifD[1]<-1
                }
                if(adjL[[d]]$adjDofi[flankpos[f],2]==1){
                    if(adjL[[d]]$adjDoi[flankpos[f]+1,1]==0){
                        fifD[1]<-1
                    }else{
                        fifD[2]<-1
                    }
                }
                if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                   adjL[[d]]$adjf2[flankpos[f],2]==0){
                    fifD[2]<-1
                }
            }else if(f==length(flankpos) &
                     flankpos[f]==length(rankelem)){ ## flank at rear end
                if(adjL[[d]]$adjDofi[flankpos[f],2]==1){
                    fifD[length(fifD)]<-1
                }
                if(adjL[[d]]$adjDofi[flankpos[f],1]==1){
                    if(adjL[[d]]$adjDoi[flankpos[f]-1,2]==0){
                        fifD[length(fifD)]<-1
                    }else{
                        fifD[length(fifD)-1]<-1
                    }
                }
                if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                   adjL[[d]]$adjf2[flankpos[f],1]==0){
                    fifD[length(fifD)-1]<-1
                }
            }else{ ## flank in middle
                if(f==1){
                    if(adjL[[d]]$adjDofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]-1,2]==0){
                        fifD[1]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]+1,1]==0){
                        fifD[1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==1){
                        fifD[1]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjf2[flankpos[f],2]==1){
                        fifD[1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                       adjL[[d]]$adjf2[flankpos[f],2]==0){
                        fifD[2]<-1
                    }
                }else if(f==length(flankpos)){
                    if(adjL[[d]]$adjDofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]+1,1]==0){
                        fifD[length(fifD)]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]-1,2]==0){
                        fifD[length(fifD)]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==1){
                        fifD[length(fifD)]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjf2[flankpos[f],1]==1){
                        fifD[length(fifD)]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                       adjL[[d]]$adjf2[flankpos[f],1]==0){
                        fifD[length(fifD)-1]<-1
                    }
                }else{
                    if(adjL[[d]]$adjDofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]-1,2]==0){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjDoi[flankpos[f]+1,1]==0){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==1){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==1){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],1]==1 &
                       adjL[[d]]$adjf2[flankpos[f],2]==1){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjDofi[flankpos[f],2]==1 &
                       adjL[[d]]$adjf2[flankpos[f],1]==1){
                        fifD[f*2-1]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],1]==0 &
                       adjL[[d]]$adjf2[flankpos[f],1]==0){
                        fifD[f*2-2]<-1
                    }
                    if(adjL[[d]]$adjf1[flankpos[f],2]==0 &
                       adjL[[d]]$adjf2[flankpos[f],2]==0){
                        fifD[f*2]<-1
                    }
                }
            }
        }

        ## store elements to mask
        for(i in 1:length(fifA)){
            if(fifA[i]==0){
                next
            }
            if(i%%2==1){ ## flank, (i+1)/2 gives flankpos
                maskA[flankpos[(i+1)/2],d]<-1
            }else{ ## insert, between flanks i/2 and i/2+1
                maskA[(flankpos[i/2]+1):(flankpos[i/2+1]-1),d]<-1
            }
        }
        for(i in 1:length(fifD)){
            if(fifD[i]==0){
                next
            }
            if(i%%2==1){ ## flank, (i+1)/2 gives flankpos
                maskD[flankpos[(i+1)/2],d]<-1
            }else{ ## insert, between flanks i/2 and i/2+1
                maskD[(flankpos[i/2]+1):(flankpos[i/2+1]-1),d]<-1
            }
        }

        ## for testing
        if(sum(fifA)==0 | sum(fifD)==0){
            stop("duplicated elements exist but neither flanks nor inserts were pre-masked")
        }

    }


    A<-rowSums(maskA)!=0
    D<-rowSums(maskD)!=0

    return(list(A=A,D=D))
}
## decide which order to keep later in tagTP2
## >>>> function will need a bit more testing - so far tested with
##      TOY24_nested_markers_test1.txt,
##      TOY24_markers_test1.txt, and
##      TOY24_markers_test2.txt (here, makes more tags than might
##        be necessary because each dupli is checked separately,
##        independent on tags for other dupli's)

## ------------------------------------------------------------------


## get matrices that have marker position of breakpoints starts and ends
##  from matrix with tags (matrix has to be a subset of one scaffold,
##  and all gaps in tags indicate breakpoints)
##  tags for breakpoints will get same value as original tags
getBreakpntsSE<-function(submat){

    if(!is.matrix(submat)){
        stop("Require matrix as input")
    }
    bptS<-matrix(0,nrow=nrow(submat),ncol=0)
    bptE<-bptS
    if(ncol(submat)>0){
        for(i in 1:ncol(submat)){
            bstart<-which(diff(c(0,submat[,i]))!=0)
            bend<-which(diff(c(submat[,i],0))!=0)

            ## if(length(bstart)>0){
            ##     bptS<-cbind(bptS,numeric(nrow(submat)))
            ##     bptS[bstart,ncol(bptS)]<-submat[bstart,i]
            ## }
            ## if(length(bend)>0){
            ##     bptE<-cbind(bptE,numeric(nrow(submat)))
            ##     bptE[bend,ncol(bptE)]<-submat[bend,i]
            ## }
            ## >>> changed here to make sure that dimension of matrix with
            ##     breakpoints corresponds to matrix with tag values
            bptS<-cbind(bptS,numeric(nrow(submat)))
            if(length(bstart)>0){
                bptS[bstart,i]<-submat[bstart,i]
            }
            bptE<-cbind(bptE,numeric(nrow(submat)))
            if(length(bend)>0){
                bptE[bend,i]<-submat[bend,i]
            }
        }
    }
    return(list(bptS=bptS,bptE=bptE))
}

## ------------------------------------------------------------------


## append rows to SYNT (adjust column dimensions if necessary)
appendData<-function(oldMatrix,matrixAppendix){

    if(!is.matrix(oldMatrix) | !is.matrix(matrixAppendix)){
        stop("Require matrix as input")
    }
    if(ncol(oldMatrix)!=ncol(matrixAppendix)){
        if(ncol(oldMatrix)>ncol(matrixAppendix)){
            matrixAppendix<-cbind(matrixAppendix,
                           matrix(0,nrow=nrow(matrixAppendix),
                                  ncol=ncol(oldMatrix)-ncol(matrixAppendix)))
        }else{
            oldMatrix<-cbind(oldMatrix,
                                matrix(0,nrow=nrow(oldMatrix),
                                       ncol=ncol(matrixAppendix)-ncol(oldMatrix)))
        }
    }
    newMatrix<-rbind(oldMatrix,matrixAppendix)
    return(newMatrix)
}

## ------------------------------------------------------------------


## if all block elements received a tp tag, adjust so that
##  largest block won't get a tag
## blocks and mask are at tested blocklevel, blocks1 is at
##  blocklevel==1
adjustTPtags<-function(tpmat,blocks,mask,blocks1,elemrows,remWgt){

    if(!is.matrix(tpmat)){
        stop("Require matrix as input")
    }
    if(sum(rowSums(tpmat)>=1)!=nrow(tpmat)){
        return(tpmat)
    }else{
        untag<-integer()
        ## get block sizes
        blsize<-integer(nrow(blocks))
        for(k in 1:length(blsize)){
            ## expand block over elements to markers
            ## get columns in elemrows
            tmp<-(blocks[k,1]):(blocks[k,2])
            ## get number of markers
            for(j in tmp){
                blsize[k]<-blsize[k]+diff(elemrows[,j])+1
            }
        }
        untag<-which(blsize==max(blsize))
        if(length(untag)>1){
            if(sum(!mask[untag])==1){
                ## only one of the largest blocks is unmasked
                untag<-untag[!mask[untag]]
            }else{
                ## get number of subblocks
                nsubbl<-integer(length(blsize))
                for(k in 1:length(blsize)){
                    nsubbl[k]<-sum(blocks1[,1]>=blocks[k,1] &
                                       blocks1[,2]<=blocks[k,2])
                }
                ## get the one(s) with the fewest subblocks
                ## >>>> ! excluding blocks that were removed in
                ##        checkAdjComplexRearr() if called from there
                nsubbl<-nsubbl[untag]
                untag<-untag[which(nsubbl==min(nsubbl))]
            }
        }
        if(length(untag)==1){
            ## get row and column in tpmat
            tagrow<-(blocks[untag,1]):(blocks[untag,2])
            tagcln<-which(colSums(tpmat[tagrow,,drop=FALSE])==length(tagrow))
            if(length(tagcln)==1){
                tagpos<-tpmat[,tagcln]==1
                tpmat[tagpos,tagcln]<-remWgt
                ## note: this might remove tags beyond the
                ##  boundaries of the largest block (which should
                ##  be reasonable as everything tagged within a
                ##  column should have correct adjacencies within it)
                if(ncol(tpmat)==2){
                    tagpos<-tpmat[,-tagcln]==1
                    tpmat[tagpos,-tagcln]<-1-remWgt
                }
            }
        }else if(length(untag)==2 & ncol(tpmat)==2){
            ## two blocks, both tagged and of equal size
            tpmat[tpmat[,1]==1,1]<-0.5
            tpmat[tpmat[,2]==1,2]<-0.5
            if(sum(rowSums(tpmat)==0.5)!=nrow(tpmat)){
                stop("Tags for transposed blocks are not as they should be")
            }
        }
        ## note: if length(untag)!=1, don't adjust tags,
        ##  unless only two adjacency sets
        return(tpmat)
    }
}

## ------------------------------------------------------------------


## if no final block can be made, exclude smallest blocks
##   and test if issue can be solved; return blocks to keep
keepBlocks<-function(blocksL,maskL,elemrows,allelem,leafelem,testlim=100){

    ## block columns:
    ##   start, end, start rank element, end rank element,
    ##   count of leaf markers, order a-/descending, block rank

    z<-length(blocksL)

    if(nrow(blocksL[[z]])==1){ ## nothing to do
        return(TRUE)
    }


    ## number of markers in block
    blsize<-integer(nrow(blocksL[[z]]))
    ## expand block over elements to markers
    for(k in 1:nrow(blocksL[[z]])){
        ## get columns in elemrows
        tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
        ## get number of markers
        for(j in tmp){
            blsize[k]<-blsize[k]+diff(elemrows[,j])+1
        }
    }
    blsizeAdj<-blsize

    duplibl<-which(maskL[[z]]==TRUE)
    if(length(duplibl)>0 & length(duplibl)<length(maskL[[z]])){
        ## some, but not all are masked (i.e., duplicated elements
        ##   that are not bound into larger blocks)
        ## set blsize for duplibl to zero for first round
        ##   of testing below
        blsizeAdj[duplibl]<-0
    }

    ## exclude small block(s) and re-cluster remaining blocks
    ##   until final simplification can be achieved;
    ## however do not enter loop if <=3 blocks remain, as
    ##   they always can be clustered
    minblsize<-min(blsizeAdj)
    tokeep<-blsizeAdj>minblsize
    final<-FALSE
    while(final==FALSE & sum(tokeep)>3){
        ## make new block with kept rows of old block
        tmpblocks<-blocksL[[z]][tokeep,]
        ## get elements within largest blocks, and adjust
        ##   positions in tmpblocks
        largeelem<-integer()
        tmpstart<-1
        for(k in 1:nrow(tmpblocks)){
            ## get columns in elemrows
            largeelem<-c(largeelem,
                         (tmpblocks[k,1]):(tmpblocks[k,2]))
            tmpend<-(tmpblocks[k,2]-tmpblocks[k,1])+tmpstart
            tmpblocks[k,1:2]<-c(tmpstart,tmpend)
            tmpstart<-tmpend+1
        }
        tmpblocks[,7]<-myRank(rowMeans(tmpblocks[,3:4,drop=FALSE]))
        ## continue clustering with reduced set of block elements
        largeblocksL<-list(tmpblocks)
        largeblocksL<-makeBlocks(largeblocksL,
                                 tmpblocks[,1],
                                 tmpblocks[,2],
                                 tmpblocks[,7],
                                 leafelem[largeelem],
                                 iter=2)
        if(nrow(largeblocksL[[length(largeblocksL)]])==1){
            ## final simplification
            final<-TRUE
        }else{
            ## exclude next smallest block(s)
            minblsize<-min(blsizeAdj[tokeep])
            tokeep<-blsizeAdj>minblsize
        }
    }
    ## for each excluded block, put it back and test
    ##   if final simplification still possible
    putback<-integer()
    if(sum(tokeep==FALSE)>1){
        cand<-which(tokeep==FALSE)
        for(cn in cand){
            totest<-tokeep
            totest[cn]<-TRUE
            ## make new block with kept rows of old block
            tmpblocks<-blocksL[[z]][totest,]
            ## get elements within remaining blocks, and adjust
            ##   positions in tmpblocks
            largeelem<-integer()
            tmpstart<-1
            for(k in 1:nrow(tmpblocks)){
                ## get columns in elemrows
                largeelem<-c(largeelem,
                             (tmpblocks[k,1]):(tmpblocks[k,2]))
                tmpend<-(tmpblocks[k,2]-tmpblocks[k,1])+tmpstart
                tmpblocks[k,1:2]<-c(tmpstart,tmpend)
                tmpstart<-tmpend+1
            }
            tmpblocks[,7]<-myRank(rowMeans(tmpblocks[,3:4,drop=FALSE]))
            ## continue clustering with reduced set of block elements
            largeblocksL<-list(tmpblocks)
            largeblocksL<-makeBlocks(largeblocksL,
                                     tmpblocks[,1],
                                     tmpblocks[,2],
                                     tmpblocks[,7],
                                     leafelem[largeelem],
                                     iter=2)
            if(nrow(largeblocksL[[length(largeblocksL)]])==1){
                ## final simplification
                putback<-c(putback,cn)
            }
        }
    }
    if(length(putback)==0){ ## none can be put back
        ##    -> just stick with set in tokeep
        tokeep<-tokeep
    }else if(length(putback)==1){  ## one can be put back
        ##    -> put it back
        tokeep[putback]<-TRUE
    }else{ ## some or all can be put back
        ##    -> test pairs / randomly select one
        ##       (better: keep testing combinations, and
        ##        also check effects on final rearrangement
        ##        numbers for different sets of removed blocks)
        tuple<-2
        putback1<-as.matrix(t(putback))
        putback2<-matrix(NA,nrow=tuple,ncol=0)
        foundcand<-FALSE
        while(foundcand==FALSE & tuple<=sum(!tokeep)){
            ## make temporary matrix with all tests
            cn1cn2<-putback2
            for(i in 1:ncol(putback1)){
                cn1<-putback1[,i]
                cn2set<-which(tokeep==FALSE)
                cn2set<-cn2set[!is.element(cn2set,cn1)]
                for(cn2 in cn2set){
                    cn1cn2<-cbind(cn1cn2,c(cn1,cn2))
                }
            }
            ## randomly select subset if too many test sets
            if(ncol(cn1cn2)>testlim){
                redset<-sample(1:ncol(cn1cn2),testlim)
                ## (standard sample is save)
                cn1cn2<-cn1cn2[,redset]
            }
            for(i in 1:ncol(cn1cn2)){
                totest<-tokeep
                totest[cn1cn2[,i]]<-TRUE

                ## make new block with kept rows of old block
                tmpblocks<-blocksL[[z]][totest,]
                ## get elements within remaining blocks, and adjust
                ##   positions in tmpblocks
                largeelem<-integer()
                tmpstart<-1
                for(k in 1:nrow(tmpblocks)){
                    ## get columns in elemrows
                    largeelem<-c(largeelem,
                                 (tmpblocks[k,1]):(tmpblocks[k,2]))
                    tmpend<-(tmpblocks[k,2]-tmpblocks[k,1])+tmpstart
                    tmpblocks[k,1:2]<-c(tmpstart,tmpend)
                    tmpstart<-tmpend+1
                }
                tmpblocks[,7]<-myRank(rowMeans(tmpblocks[,3:4,drop=FALSE]))
                ## continue clustering with reduced set of block elements
                largeblocksL<-list(tmpblocks)
                largeblocksL<-makeBlocks(largeblocksL,
                                         tmpblocks[,1],
                                         tmpblocks[,2],
                                         tmpblocks[,7],
                                         leafelem[largeelem],
                                         iter=2)
                if(nrow(largeblocksL[[length(largeblocksL)]])==1){
                    ## final simplification
                    putback2<-cbind(putback2,sort(cn1cn2[,i]))
                }
            }
            putback2<-putback2[,!duplicated(putback2,MARGIN=2),drop=FALSE]
            if(ncol(putback2)==0){
                putback2<-putback1
                tuple<-tuple-1
                foundcand<-TRUE
            }else if(ncol(putback2)==1){
                foundcand<-TRUE
            }else{
                putback1<-putback2
                tuple<-tuple+1
                putback2<-matrix(NA,nrow=tuple,ncol=0)
            }
        }
        if(ncol(putback2)==1){
            ## one set can be put back
            tokeep[putback2[,1]]<-TRUE
        }else if(ncol(putback2)>1){
            ## randomly sample one pair from largest putback2
            candsz<-apply(putback2,2,function(x) sum(blsizeAdj[x]))
            candps<-(1:ncol(putback2))[candsz==max(candsz)]
            rsm<-candps[sample.int(length(candps),1)]
            ## standard sample is NOT save if length of set can be 1
            ## (which should not be possible)
            tokeep[putback2[,rsm]]<-TRUE
        }
    }
    ## when <3 blocks kept, keep largest /
    ##   randomly sample three blocks
    if(sum(tokeep==TRUE)<3){
        maxblsize<-max(blsizeAdj)
        tokeep<-blsizeAdj==maxblsize
        while(sum(tokeep==TRUE)<3){
            ## get next largest block(s)
            maxblsize<-max(blsizeAdj[!tokeep])
            tokeep<-blsizeAdj>=maxblsize
        }
        if(sum(tokeep==TRUE)>3){
            ## got too many, remove some randomly
            cand<-which(blsizeAdj==maxblsize)
            sampsize<-sum(tokeep==TRUE)-3
            getout<-cand[sample.int(length(cand),sampsize)]
            ## standard sample is NOT save to use here
            ##   because length(cand) can be 1
            tokeep[getout]<-FALSE
        }
    }
    return(tokeep)
}

## ------------------------------------------------------------------


## if no final block can be made, exclude smallest blocks
##   and test for wrong adjacencies when re-clustering blocks
##   so that final combination is possible; tag each original
##   block in-between incorrect adjacencies separately, as well
##   as all excluded blocks
checkAdjComplexRearr<-function(blocksL,maskL,allelem,elemrows,
                               leafelem,leaves,testorientation,
                               tokeep,remWgt,ori=NULL){

    if(sum(tokeep==TRUE)==0 | sum(tokeep==FALSE)==0){
        stop("tokeep requires a mix of TRUE and FALSE values")
    }
    if(!is.null(ori)){
        if(!is.element(ori,c(1,-1))){
            stop("Unexpected orientation")
        }
    }

    tmptp<-matrix(NA,nrow=length(allelem),ncol=0)
    tmpiv<-matrix(NA,nrow=length(allelem),ncol=0)
    invelem<-rep(NA,length(allelem))
    splitblockid<-rep("",length(allelem))

    z<-length(blocksL)

    ## exclude blocks from original blocksL and maskL for
    ##  temporary use (i.e., for assignOri3 below)
    keptblocksL<-vector("list",z)
    keptmaskL<-vector("list",z)
    keptblocksL[[z]]<-blocksL[[z]][tokeep,,drop=FALSE]
    keptmaskL[[z]]<-maskL[[z]][tokeep]
    idx1<-blocksL[[z]][!tokeep,1,drop=FALSE]
    idx2<-blocksL[[z]][!tokeep,2,drop=FALSE]
    bllev<-z-1
    while(bllev>=1){
        toexcl<-integer()
        for(x in 1:length(idx1)){
            toexcl<-c(toexcl,which(blocksL[[bllev]][,1]>=idx1[x] &
                                       blocksL[[bllev]][,2]<=idx2[x]))
        }
        keptblocksL[[bllev]]<-blocksL[[bllev]][-toexcl,,drop=FALSE]
        keptmaskL[[bllev]]<-maskL[[bllev]][-toexcl]
        bllev<-bllev-1
    }
    ## adjust positions in keptblocksL (important)
    bllev<-z
    while(bllev>=1){
        tmpstart<-1
        for(k in 1:nrow(keptblocksL[[bllev]])){
            tmpend<-(keptblocksL[[bllev]][k,2]-
                         keptblocksL[[bllev]][k,1])+tmpstart
            keptblocksL[[bllev]][k,1:2]<-c(tmpstart,tmpend)
            tmpstart<-tmpend+1
        }
        bllev<-bllev-1
    }


    ## make new block with kept rows of old block
    tmpblocks<-blocksL[[z]][tokeep,,drop=FALSE]
    ## get elements within remaining blocks, and adjust
    ##   positions in tmpblocks (important)
    unmaskedelem<-integer()
    tmpstart<-1
    for(k in 1:nrow(tmpblocks)){
        ## get columns in elemrows
        unmaskedelem<-c(unmaskedelem,
                        (tmpblocks[k,1]):(tmpblocks[k,2]))
        tmpend<-(tmpblocks[k,2]-tmpblocks[k,1])+tmpstart
        tmpblocks[k,1:2]<-c(tmpstart,tmpend)
        tmpstart<-tmpend+1
    }
    tmpblocks[,7]<-myRank(rowMeans(tmpblocks[,3:4,drop=FALSE]))
    ## continue clustering with reduced set of block elements
    unmaskedblocksL<-list(tmpblocks)
    unmaskedblocksL<-makeBlocks(unmaskedblocksL,
                                tmpblocks[,1],
                                tmpblocks[,2],
                                tmpblocks[,7],
                                leafelem[unmaskedelem],
                                iter=2)
    if(nrow(unmaskedblocksL[[length(unmaskedblocksL)]])!=1){
        stop("Something went wrong when excluding blocks")
    }


    ## find duplicated elements not yet bound into larger blocks
    tmpelem<-allelem[unmaskedelem]
    tmpdupli<-unique(tmpelem[duplicated(tmpelem)])
    ## masks used in assignOri3, checkAdj*, adjustTPtags
    ##unmaskedmaskL<-setMasks(unmaskedblocksL,tmpelem,tmpdupli,
    ##                        maskL[[z]][tokeep])
    unmaskedmaskL<-setMasks(unmaskedblocksL,tmpelem,tmpdupli)
    ## get correct subset for leaves and testorientation
    ## expand block over elements to markers
    unmaskedmarkers<-integer()
    for(j in unmaskedelem){
        ## get position of markers
        unmaskedmarkers<-c(unmaskedmarkers,
                           (elemrows[1,j]):(elemrows[2,j]))
    }
    ## get correct from-to positions for elements
    unmaskedelemrows<-matrix(1,nrow=2,ncol=length(unmaskedelem))
    cntr<-0
    for(i in 1:length(unmaskedelem)){
        cntr<-cntr+1
        unmaskedelemrows[1,i]<-cntr
        cntr<-cntr+(elemrows[2,unmaskedelem[i]]-elemrows[1,unmaskedelem[i]])
        unmaskedelemrows[2,i]<-cntr
    }

    ## bind keptblocksL and unmaskedblocksL (same for maskL)
    allblocksL<-keptblocksL
    allmaskL<-keptmaskL
    allblocksL[[z]]<-unmaskedblocksL[[1]]
    allmaskL[[z]]<-unmaskedmaskL[[1]]
    for(bl in 2:length(unmaskedblocksL)){
        allblocksL[[z+(bl-1)]]<-unmaskedblocksL[[bl]]
        allmaskL[[z+(bl-1)]]<-unmaskedmaskL[[bl]]
    }

    ## identify top-level orientation
    if(is.null(ori)){
        unmaskedori<-assignOri3(allblocksL,allmaskL,
                                allelem[unmaskedelem],
                                unmaskedelemrows,
                                leafelem[unmaskedelem],
                                leaves[unmaskedmarkers],
                                testorientation[unmaskedmarkers])
    }else{
        unmaskedori<-ori
    }
    ## check adjacencies
    allevel<-length(allblocksL)-1
    if(unmaskedori==1){
        unmaskedtp<-checkAdjAscend(allblocksL[[allevel]],
                                   allmaskL[[allevel]],
                                   length(allelem[unmaskedelem]))
    }else if(unmaskedori== -1){
        unmaskedtp<-checkAdjDescend(allblocksL[[allevel]],
                                    allmaskL[[allevel]],
                                    length(allelem[unmaskedelem]))
    }else{
        stop("Node has unexpected orientation")
    }

    if(ncol(unmaskedtp)>0){
        ## if all block elements received a tp tag, adjust so that
        ##  largest block won't get a tag
        if(sum(rowSums(unmaskedtp)>=1)==nrow(unmaskedtp)){
            unmaskedtp<-adjustTPtags(unmaskedtp,
                                     allblocksL[[allevel]],
                                     allmaskL[[allevel]],
                                     allblocksL[[1]],
                                     unmaskedelemrows,
                                     remWgt)
        }
    }

    ## storage for block element orientation
    ##  (can be overwritten in functions below)
    allblockoriL<-vector("list",length(allblocksL))
    for(az in 1:length(allblocksL)){
        allblockoriL[[az]]<-allblocksL[[az]][,6]
    }


    ## store splitblockid (if TP in 'checkOriAscend'/'checkOriDescend')
    unmaskedsplitblockid<-rep("",length(unmaskedelem))

    ## storage for inversions
    unmaskediv<-matrix(NA,nrow=length(unmaskedelem),ncol=0)

    ## for using allblocksL below it is important that
    ##   columns 1,2,5,6 are correct (i.e., match to subset
    ##   of allelem[unmaskedelem])


    ## in checkOri* functions, special treatment for
    ##   blocklevel==2: however here checks will only be
    ##   valid if alllevel==1 corresponds to
    ##   blocklevel==1; if this is not the case, checkOri*
    ##   needs to work as when the blocklevel would be >2

    ## === check orientation of blocks (top level) ===
    if(unmaskedori==1){
        xxx<-checkOriAscend(allblocksL,allevel,
                            allelem[unmaskedelem],tmpdupli,
                            unmaskedelemrows,
                            leaves[unmaskedmarkers],
                            testorientation[unmaskedmarkers],
                            allblockoriL,
                            subbl=1:nrow(allblocksL[[allevel]]),
                            unmaskedsplitblockid,remWgt=remWgt)
        ## returns tpElA, ivElA, (modified) blockoriL, splitblockid
        unmaskedtp<-cbind(unmaskedtp,xxx$tpElA)
        unmaskediv<-cbind(unmaskediv,xxx$ivElA)
        allblockoriL<-xxx$blockoriL
        unmaskedsplitblockid<-xxx$splitblockid
        ## keep track of expected orientation of leaves
        unmaskedinvelem<-rep(0,length(unmaskedelem))
        if(ncol(xxx$ivElA)>0){
            for(i in 1:ncol(xxx$ivElA)){
                newinv<-which(xxx$ivElA[,i]==1 & unmaskedinvelem==0)
                backinv<-which(xxx$ivElA[,i]==1 & unmaskedinvelem==1)
                if(length(newinv)>0){
                    unmaskedinvelem[newinv]<-rep(1,length(newinv))
                }
                if(length(backinv)>0){
                    unmaskedinvelem[backinv]<-rep(0,length(backinv))
                }
            }
        }
    }else if(unmaskedori == -1){
        xxx<-checkOriDescend(allblocksL,allevel,
                             allelem[unmaskedelem],tmpdupli,
                             unmaskedelemrows,
                             leaves[unmaskedmarkers],
                             testorientation[unmaskedmarkers],
                             allblockoriL,
                             subbl=1:nrow(allblocksL[[allevel]]),
                             unmaskedsplitblockid,remWgt=remWgt)
        ## returns tpElD, ivElD, (modified) blockoriL,splitblockid
        unmaskedtp<-cbind(unmaskedtp,xxx$tpElD)
        unmaskediv<-cbind(unmaskediv,xxx$ivElD)
        allblockoriL<-xxx$blockoriL
        unmaskedsplitblockid<-xxx$splitblockid
        ## keep track of expected orientation of leaves
        unmaskedinvelem<-rep(1,length(unmaskedelem))
        if(ncol(xxx$ivElD)>0){
            for(i in 1:ncol(xxx$ivElD)){
                newinv<-which(xxx$ivElD[,i]==1 & unmaskedinvelem==0)
                backinv<-which(xxx$ivElD[,i]==1 & unmaskedinvelem==1)
                if(length(newinv)>0){
                    unmaskedinvelem[newinv]<-rep(1,length(newinv))
                }
                if(length(backinv)>0){
                    unmaskedinvelem[backinv]<-rep(0,length(backinv))
                }
            }
        }
    }else{ ## unmaskedori==9 (should never be the case here)
        unmaskedinvelem<-rep(NA,length(unmaskedelem))
    }
    ## ===

    allevel<-allevel-1

    ## === check orientation of blocks (remaining levels) ===
    ## loop through remaining levels, separately for each higher-level
    ##  block (adjacencies will always be correct by definition)
    ## don't do the last level as then the tests
    ##  should switch to blocksL[[blocklevel]]

    if(allevel>z){
        for(az in allevel:(z+1)){
            for(k in 1:(nrow(allblocksL[[az+1]]))){
                ## get elements to consider
                tmp<-which(allblocksL[[az]][,1]>=
                               allblocksL[[az+1]][k,1] &
                                   allblocksL[[az]][,2]<=
                                       allblocksL[[az+1]][k,2])
                passedOri<-9
                if(allblockoriL[[az+1]][k]==9){
                    ## pass orientation from higher level
                    y<-az+2
                    while(y<=length(allblocksL) & passedOri==9){
                        if(nrow(allblocksL[[y]])==1){
                            ## avoid taking original orientation
                            ##  but take assigned ori instead
                            break
                        }
                        tmp2<-which(allblocksL[[y]][,1]<=
                                        allblocksL[[az+1]][k,1] &
                                            allblocksL[[y]][,2]>=
                                                allblocksL[[az+1]][k,2])
                        passedOri<-allblockoriL[[y]][tmp2]
                        y<-y+1
                    }
                    ## if still 9, use ori (if ori would be 9 too then
                    ##  this loop would never have been entered)
                    if(passedOri==9){
                        passedOri<-unmaskedori
                    }
                    allblockoriL[[az+1]][k]<-passedOri
                }
                ## get the ones with opposite of expected orientation
                if(allblockoriL[[az+1]][k]==1){
                    xxx<-checkOriAscend(allblocksL,az,
                                        allelem[unmaskedelem],tmpdupli,
                                        unmaskedelemrows,
                                        leaves[unmaskedmarkers],
                                        testorientation[unmaskedmarkers],
                                        allblockoriL,
                                        subbl=tmp,
                                        unmaskedsplitblockid,remWgt=remWgt)
                    ## returns tpElA, ivElA, (modified) blockoriL,
                    ##  splitblockid
                    unmaskedtp<-cbind(unmaskedtp,xxx$tpElA)
                    unmaskediv<-cbind(unmaskediv,xxx$ivElA)
                    allblockoriL<-xxx$blockoriL
                    unmaskedsplitblockid<-xxx$splitblockid
                    ## keep track of expected orientation of leaves
                    if(ncol(xxx$ivElA)>0){
                        for(i in 1:ncol(xxx$ivElA)){
                            newinv<-which(xxx$ivElA[,i]==1 &
                                              unmaskedinvelem==0)
                            backinv<-which(xxx$ivElA[,i]==1 &
                                               unmaskedinvelem==1)
                            if(length(newinv)>0){
                                unmaskedinvelem[newinv]<-rep(1,length(newinv))
                            }
                            if(length(backinv)>0){
                                unmaskedinvelem[backinv]<-rep(0,length(backinv))
                            }
                        }
                    }
                }else if(allblockoriL[[az+1]][k]== -1){
                    xxx<-checkOriDescend(allblocksL,az,
                                         allelem[unmaskedelem],tmpdupli,
                                         unmaskedelemrows,
                                         leaves[unmaskedmarkers],
                                         testorientation[unmaskedmarkers],
                                         allblockoriL,
                                         subbl=tmp,
                                         unmaskedsplitblockid,remWgt=remWgt)
                    ## returns tpElD, ivElD, (modified) blockoriL,
                    ##  splitblockid
                    unmaskedtp<-cbind(unmaskedtp,xxx$tpElD)
                    unmaskediv<-cbind(unmaskediv,xxx$ivElD)
                    allblockoriL<-xxx$blockoriL
                    unmaskedsplitblockid<-xxx$splitblockid
                    ## keep track of expected orientation of leaves
                    if(ncol(xxx$ivElD)>0){
                        for(i in 1:ncol(xxx$ivElD)){
                            newinv<-which(xxx$ivElD[,i]==1 &
                                              unmaskedinvelem==0)
                            backinv<-which(xxx$ivElD[,i]==1 &
                                               unmaskedinvelem==1)
                            if(length(newinv)>0){
                                unmaskedinvelem[newinv]<-rep(1,length(newinv))
                            }
                            if(length(backinv)>0){
                                unmaskedinvelem[backinv]<-rep(0,length(backinv))
                            }
                        }
                    }
                }
            } ## close loop over rows
        } ## close loop over allevels
    } ## close allevel>z

    ## allblocksL[[z]] is not tested here


    ## ===== transfer identified ori to full blocksL =====

    ## identify expected blockorientation for each row in
    ##   allblocksL[[z]] (and then blocksL[[z]])
    for(k in 1:(nrow(allblocksL[[z+1]]))){
        ## get elements to consider
        tmp<-which(allblocksL[[z]][,1]>=
                       allblocksL[[z+1]][k,1] &
                           allblocksL[[z]][,2]<=
                               allblocksL[[z+1]][k,2])
        passedOri<-9
        if(allblockoriL[[z+1]][k]==9){
            ## pass orientation from higher level
            y<-z+2
            while(y<=length(allblocksL) & passedOri==9){
                if(nrow(allblocksL[[y]])==1){
                    ## avoid taking original orientation
                    ##  but take assigned ori instead
                    break
                }
                tmp2<-which(allblocksL[[y]][,1]<=
                                allblocksL[[z+1]][k,1] &
                                    allblocksL[[y]][,2]>=
                                        allblocksL[[z+1]][k,2])
                passedOri<-allblockoriL[[y]][tmp2]
                y<-y+1
            }
            ## if still 9, use ori (if ori would be 9 too then
            ##  this loop would never have been entered)
            if(passedOri==9){
                passedOri<-unmaskedori
            }
            allblockoriL[[z+1]][k]<-passedOri
        }
        ## transfer ori to first-level blocks (this will be the expected ori)
        allblockoriL[[z]][tmp]<-allblockoriL[[z+1]][k]
    }

    ## transfer expected ori to elements in blocksL[[z]]
    expectedOri<-rep(9,length(tokeep))
    expectedOri[tokeep]<-allblockoriL[[z]]
    ## as simplification, assign unmaskedori to excluded elements
    ##   (>>>> better would be to find the first clustering in
    ##         allblocksL where they would integrate and
    ##         take that ori, but would be a bit of work)
    expectedOri[!tokeep]<-unmaskedori

    ## ===== transfer identified rearrs to full blocksL =====

    ## expand unmaskedtp to all elements (will result in gaps
    ##   of tags when excluded block is crossed)
    if(ncol(unmaskedtp)>0){
        for(k in 1:ncol(unmaskedtp)){
            tpEl<-rep(0,length(allelem))
            tpEl[unmaskedelem]<-unmaskedtp[,k]
            tmptp<-cbind(tmptp,tpEl)
        }
    }

    ## add tags for masked blocks
    maskedbl<-sort(which(tokeep==FALSE))
    ## make tags for all elements in these blocks
    for(k in maskedbl){
        tpEl<-rep(0,length(allelem))
        tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
        tpEl[tmp]<-rep(1,length(tmp))
        tmptp<-cbind(tmptp,tpEl)
    }

    ## expand unmaskediv to all elements (will result in gaps
    ##   of tags when excluded block is crossed)
    if(ncol(unmaskediv)>0){
        for(k in 1:ncol(unmaskediv)){
            ivEl<-rep(0,length(allelem))
            ivEl[unmaskedelem]<-unmaskediv[,k]
            tmpiv<-cbind(tmpiv,ivEl)
        }
    }

    ## expand unmaskedinvelem to all elements
    invelem[unmaskedelem]<-unmaskedinvelem
    ## add values for excluded elements depending on their
    ##   expected ori
    for(k in maskedbl){
        tmp<-(blocksL[[z]][k,1]):(blocksL[[z]][k,2])
        if(expectedOri[k]==1){
            invelem[tmp]<-0
        }else if(expectedOri[k]== -1){
            invelem[tmp]<-1
        }else{
            invelem[tmp]<-NA
        }
    }

    ## expand unmaskedsplitblockid to all elements
    splitblockid[unmaskedelem]<-unmaskedsplitblockid

    ## =====

    return(list(tmptp=tmptp,tmpiv=tmpiv,ori=unmaskedori,
                invelem=invelem,splitblockid=splitblockid,
                expectedOri=expectedOri))
}

## ------------------------------------------------------------------


## find orientation of block element by stepping from current
##  blocklevel all the way down to leaves, if necessary
##  (orientation assignment might still not possible, i.e.,
##   when no leaf or no orientation info available)
findOri<-function(blocksL,bllev,idx,elemrows,
                  testorientation,leaves){
    ## bllev: bocklevel to be tested
    ## idx: row with element to be tested at bllev

    if(bllev<1){
        stop("bllev needs to be at least 1")
    }
    if(length(idx)!=1){
        stop("require unique index to find element orientation")
    }

    ord<-blocksL[[bllev]][idx,6]

    while(ord==9 & bllev>0 & length(idx)==1){
        elms<-blocksL[[bllev]][idx,1:2]
        bllev<-bllev-1
        ## get index for next lower level
        if(bllev>0){
            idx<-which(blocksL[[bllev]][,1]==elms[1] &
                           blocksL[[bllev]][,2]==elms[2])
            ord<-blocksL[[bllev]][idx,6]
        }else if(bllev==0){
            ## test if leaf and if yes, get its orientation
            if(length(idx)==1){
                tmp<-(blocksL[[1]][idx,1]):(blocksL[[1]][idx,2])
                tmp2<-(elemrows[1,tmp[1]]):(elemrows[2,tail(tmp,n=1L)])
                if(length(tmp2)==1 & sum(!is.na(testorientation[tmp2]) &
                                             leaves[tmp2]==1)==1){
                    ord<-testorientation[tmp2]
                }
            }
        }
    }

    if(!is.element(ord,c(1, -1,9))){
        stop("something went wrong when trying to find orientation of element")
    }
    return(list(ord=ord,bllev=bllev))
}

## ------------------------------------------------------------------
