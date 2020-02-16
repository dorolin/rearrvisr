## ------------------------------------------------------------------------
## make summary table for synteny blocks resulting from the alignment
##  of the focal genome to the ancestral genome reconstruction
## ------------------------------------------------------------------------

#' Summarize Synteny Blocks
#'
#' For each synteny block, summarize rearrangements and information on the
#' alignment between the focal genome and the compared genome
#'
#' @param SYNT A list of matrices that store data on different classes of
#'   rearrangements and additional information. \code{SYNT} must have been
#'   generated with the \code{\link{computeRearrs}} function (optionally
#'   filtered with the \code{\link{filterRearrs}} function).
#' @param focalgenome Data frame representing the focal genome, containing the
#'   mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
#'   \code{$end}, and \code{$strand}, and optional further columns. Markers need
#'   to be ordered by their map position.
#' @param compgenome Data frame representing the compared genome (e.g., an
#'   ancestral genome reconstruction, or an extant genome), with the first three
#'   columns \code{$marker}, \code{$orientation}, and \code{$car}, followed by
#'   columns alternating node type and node element. Markers need to be ordered
#'   by their node elements. \code{compgenome} must be the same data frame that
#'   was used to generate the list \code{SYNT} with the
#'   \code{\link{computeRearrs}} function.
#' @param ordfocal Character vector with the IDs of the focal genome segments
#'   that will be summarized. Have to match (a subset of) IDs in
#'   \code{focalgenome$scaff}.
#'
#' @details
#'
#'   \code{focalgenome} must contain the column \code{$marker}, a vector of
#'   either characters or integers with unique ortholog IDs that can be matched
#'   to the values in the rownames of \code{SYNT} and the \code{$marker} column
#'   of \code{compgenome}. Values can be \code{NA} for markers that have no
#'   ortholog. \code{$scaff} must be a character vector giving the name of the
#'   focal genome segment (i.e., chromosome or scaffold) of origin of each
#'   marker. \code{$start} and \code{$end} must be numeric vectors giving the
#'   location of each marker on its focal genome segment. \code{$strand} must be
#'   a vector of \code{"+"} and \code{"-"} characters giving the reading
#'   direction of each marker. Additional columns are ignored and may store
#'   custom information, such as marker names. Markers need to be ordered by
#'   their map position within each focal genome segment, for example by running
#'   the \code{\link{orderGenomeMap}} function. \code{focalgenome} may contain
#'   additional rows that were absent when running the
#'   \code{\link{computeRearrs}} function. However, all markers present in
#'   \code{SYNT} need to be contained in \code{focalgenome}, with the subset of
#'   shared markers being in the same order.
#'
#' @section References:
#'
#'   Booth, K.S. & Lueker, G.S. (1976). Testing for the consecutive ones
#'   property, interval graphs, and graph planarity using \emph{PQ}-Tree
#'   algorithms. \emph{Journal of Computer and System Sciences}, \strong{13},
#'   335--379. doi:
#'   \href{https://doi.org/10.1016/S0022-0000(76)80045-1}{10.1016/S0022-0000(76)80045-1}.
#'
#'   Chauve, C. & Tannier, E. (2008). A methodological framework for the
#'   reconstruction of contiguous regions of ancestral genomes and its
#'   application to mammalian genomes. \emph{PLOS Computational Biology},
#'   \strong{4}, e1000234. doi:
#'   \href{https://doi.org/10.1371/journal.pcbi.1000234}{10.1371/journal.pcbi.1000234}.
#'
#' @return A list of lists that summarizes the alignment between the focal
#'   genome and each \emph{PQ-tree}, and records whether synteny blocks are part
#'   of different classes of rearrangements. The top-level list elements are
#'   focal genome segments, and the lower-level list elements contain
#'   information on synteny blocks and rearrangements for each focal genome
#'   segment. For details on \emph{PQ-trees} see the description of the
#'   \code{"compgenome"} class in the Details section of the
#'   \code{\link{checkInfile}} function, Booth & Lueker 1976, Chauve & Tannier
#'   2008, or the package vignette.
#'
#'   The names of the top-level list elements correspond to the strings in
#'   \code{ordfocal}. Each list element is itself a list containing the data
#'   frame \code{$blocks} and five numeric matrices \code{$NM1}, \code{$NM2},
#'   \code{$SM}, \code{$IV}, and \code{$IVsm}, described below. In all six
#'   list elements, each synteny block is represented by a row. Note that
#'   separate blocks are also generated when the hierarchical structure of the
#'   underlying \emph{PQ-tree} changes, therefore not all independent rows are
#'   caused by a rearrangement.
#'
#'   \code{$blocks} contains information on the alignment and structure of each
#'   \emph{PQ-tree}. The columns \code{$blocks$start} and \code{$blocks$end}
#'   give the start and end positions of the synteny block in \code{SYNT}
#'   (positions start at \code{1} separately for each focal genome segment).
#'   \code{$blocks$markerS} and \code{$blocks$markerE} give the marker IDs of
#'   the first and last marker per block. \code{$blocks$car} gives the ID of the
#'   CAR. Nine columns per hierarchy level describe the structure of each
#'   \emph{PQ-tree} and its alignment to the focal genome. Hierarchy levels of
#'   the \emph{PQ-trees} are indicated by suffixes \code{{1, 2, ...}}.
#'   \code{$blocks$type} gives the node type. \code{$blocks$elemS} and
#'   \code{$blocks$elemE} give the first and last ID of the node elements per
#'   block. They correspond to the IDs in the odd columns of \code{compgenome}
#'   (note that some IDs within blocks or in-between might be missing when
#'   markers in the compared genome are absent from the focal genome).
#'   \code{$blocks$node} indicates whether the block contains \emph{PQ-tree}
#'   nodes (value is \code{1}) or only leaf elements (value is \code{0}). The
#'   columns \code{$blocks$nodeori}, \code{$blocks$subnode},
#'   \code{$blocks$blockid}, \code{$blocks$blockori}, and \code{$blocks$premask}
#'   summarize for each block the values in the list elements of \code{SYNT}
#'   with the corresponding names (described in the Value section in the
#'   documentation of the \code{\link{computeRearrs}} function). The column
#'   \code{$blocks$nodeori1}, for example, summarizes for each block the values
#'   in the second column (i.e., the first node level) of \code{SYNT$nodeori}.
#'
#'   The numeric matrices \code{$NM1}, \code{$NM2}, \code{$SM}, \code{$IV},
#'   and \code{$IVsm} indicate whether blocks are part of different classes of
#'   rearrangements. \code{$NM1} stores \code{T}rans\code{L}ocations between
#'   CARs \code{B}etween focal \code{S}egments; \code{$NM2} stores
#'   \code{T}rans\code{L}ocations between CARs \code{W}ithin focal
#'   \code{S}egments; \code{$SM} stores \code{T}rans\code{L}ocations within
#'   CARs \code{W}ithin focal \code{S}egments; \code{$IV} and \code{$IVsm} store
#'   \code{I}n\code{V}ersions within CARs within focal segments. In \code{$IV},
#'   blocks that are part of a multi-marker inversion are tagged with \code{1},
#'   while in \code{$IVsm}, integers \code{>0} indicate the positions of
#'   single-marker inversions (i.e., markers with switched orientation) within
#'   their blocks. Each rearrangement is represented by a separate column, and
#'   blocks that are part of a rearrangement have a tag value of \code{>0}. Note
#'   that some columns in \code{$NM2} or \code{$SM} may be duplicated due to
#'   the functioning of the underlying algorithm in \code{\link{computeRearrs}};
#'   although corresponding to the same rearrangement, these duplicated columns
#'   are nevertheless included for completeness. By default these columns will
#'   not be visualized with the \code{\link{genomeRearrPlot}} function. If no
#'   rearrangements were detected for a certain class, the matrix has zero
#'   columns. See the package vignette or the Value section in the documentation
#'   of the \code{\link{computeRearrs}} function for details on the meaning of
#'   different tag values in these matrices. Note that if \code{SYNT} has been
#'   filtered with the \code{\link{filterRearrs}} function, only the above
#'   matrices will be affected, while the information in \code{$blocks} will
#'   remain unchanged.
#'
#'   The returned data can be visualized with the \code{\link{genomeRearrPlot}}
#'   function.
#'
#' @seealso \code{\link{checkInfile}}, \code{\link{computeRearrs}},
#'   \code{\link{filterRearrs}}, \code{\link{genomeRearrPlot}}.
#'
#' @examples
#' SYNT <- computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' BLOCKS <- summarizeBlocks(SYNT, TOY24_focalgenome, TOY24_compgenome,
#'                           c("1","2","3"))
#'
#' \dontrun{
#'
#' ## show summary for first focal genome segment
#' BLOCKS[[1]]
#' }
#'
#' @export
#' @importFrom utils head tail

summarizeBlocks<-function(SYNT,focalgenome,compgenome,ordfocal){

    ## checks
    ## -------------------------------------------

    ## ERRORS
    ## check SYNT
    checkInfile(SYNT, myclass="SYNT")
    ## check focalgenome
    checkInfile(focalgenome, myclass="focalgenome", checkorder = TRUE)
    ## check compgenome
    checkInfile(compgenome, myclass="compgenome", checkorder = FALSE)
    ## class of marker in tree and marker file must be the same
    ## (either character or integer)
    if(class(focalgenome$marker) != class(compgenome$marker)){
        stop("class for '$marker' in focalgenome and compgenome need to be identical")
    }
    ## check settings
    if(length(ordfocal) < 1){
        stop("require at least one focal genome segment in ordfocal")
    }
    if(class(ordfocal) != "character"){
        stop("class for ordfocal needs to be 'character'")
    }
    if(anyDuplicated(ordfocal) > 0){
        stop("some focal segments in ordfocal are duplicated")
    }
    if(length(intersect(ordfocal,focalgenome$scaff)) != length(ordfocal)){
        stop(paste("focal segment", setdiff(ordfocal,intersect(ordfocal,focalgenome$scaff))[1], "is absent in focalgenome"))
    }


    ## initial processing
    ## -------------------------------------------

    ## subset focalgenome to those retained in SYNT
    markerpos<-match(as.character(focalgenome$marker),rownames(SYNT$SM))
    markerpos<-which(!is.na(markerpos))
    if(length(unique(markerpos)) != nrow(SYNT$SM)){
        stop("some markers in SYNT are absent in focalgenome")
    }
    markers<-focalgenome[markerpos,,drop=FALSE]
    if(sum(rownames(SYNT$SM) == as.character(markers$marker)) != nrow(markers)){
        stop("SYNT is not ordered to the subset of shared markers\n    in focalgenome. Rerun 'computeRearrs' and retry.")
    }
    if(length(intersect(ordfocal,markers$scaff)) != length(ordfocal)){
        stop(paste("focal segment", setdiff(ordfocal,intersect(ordfocal,markers$scaff))[1], "is absent in the subset of\n    shared markers in SYNT"))
    }

    ## subset compgenome to those retained in SYNT
    treepos<-match(markers$marker,compgenome$marker)
    treepos<-treepos[!is.na(treepos)]
    if(length(unique(treepos)) != nrow(SYNT$SM)){
        stop("some markers in SYNT are absent in compgenome")
    }
    tree<-compgenome[treepos,,drop=FALSE]
    if(sum(rownames(SYNT$SM) == as.character(tree$marker)) != nrow(tree)){
        stop("could not subset compgenome to the set of shared markers\n    in SYNT. Rerun 'computeRearrs' and retry.")
    }

    if(sum(sum(markers$marker != tree$marker)) > 0){
        stop("markers in compgenome and focalgenome could not correctly be subset")
    }

    nhier<-ncol(SYNT$nodeori)

    BLOCKS<-vector("list",length(ordfocal))
    names(BLOCKS)<-ordfocal


    ## -------------------------------------------
    ## summarize blocks and rearrangements
    ## -------------------------------------------

    for(s in 1:length(ordfocal)){

        testscaff<-ordfocal[s]
        scaffrows<-which(markers$scaff==testscaff)

        ## ----------------
        ## summarize blocks
        ## ----------------

        ## get nhier for current scaffold
        nhier2<-0
        for(n in 2:nhier){
            treecol<-2+(n-1)*2 ## column in tree with node type
            if(sum(is.element(c("Q","P"),tree[scaffrows,treecol]))>0){
                nhier2<-n
            }
        }

        blocks<-matrix(NA,nrow=0,ncol=5+9*(nhier2-1))

        tmpblockid<-matrix(NA,nrow=length(scaffrows),ncol=0)
        for(n in 2:nhier){
            treecol<-1+(n-1)*2 ## column in tree with car/node id
            tmpelem1<-tree[scaffrows,treecol]
            ## set to NA those that are NA at next level (leaves)
            tmpelem2<-tree[scaffrows,treecol+2]
            tmpelem1[is.na(tmpelem2)]<-NA
            tmpblockid<-cbind(tmpblockid,tmpelem1,SYNT$blockid[scaffrows,n])
        }
        colnames(tmpblockid)<-paste0("[,",1:ncol(tmpblockid),"]")

        ## get start and end positions in markers
        lastrow<-rep(NA,ncol(tmpblockid))
        curr<-0
        for(i in 1:length(scaffrows)){
            if(identical(lastrow,tmpblockid[i,])==TRUE){
                blocks[curr,2]<-i ## end
            }else{ ## new block
                blocks<-rbind(blocks,rep(NA,ncol(blocks)))
                curr<-curr+1
                blocks[curr,1]<-i ## start
                blocks[curr,2]<-i ## end
                lastrow<-tmpblockid[i,]
            }
        }

        for(k in 1:nrow(blocks)){
            ## get information from tree and SYNT
            tmprows<-scaffrows[(blocks[k,1]):(blocks[k,2])]
            blocks[k,5]<-tree$car[tmprows[1]]
            for(n in 2:nhier2){
                treecol<-2+(n-1)*2 ## column in tree with node type
                ## the following only if there are still elements
                if(sum(is.element(c("Q","P"),tree[tmprows,treecol]))>0){
                    ## start of columns for that hierarchy level
                    blockcol<-6+(n-2)*9
                    ## node type
                    blocks[k,blockcol]<-tree[tmprows[1],treecol]
                    ## start and end element for Q-nodes
                    if(is.element("Q",blocks[k,blockcol])==TRUE){
                        blocks[k,blockcol+1]<-tree[head(tmprows,n=1L),
                                                   treecol+1]
                        blocks[k,blockcol+2]<-tree[tail(tmprows,n=1L),
                                                   treecol+1]
                    }
                    ## contains nodes
                    blocks[k,blockcol+3]<-sum(is.element(c("Q","P"),
                                                         tree[tmprows,
                                                                 treecol+2]))
                    ## nodeori
                    blocks[k,blockcol+4]<-unique(SYNT$nodeori[tmprows,n])
                    ## subnode
                    blocks[k,blockcol+5]<-unique(SYNT$subnode[tmprows,n-1])
                    ## blockid
                    blocks[k,blockcol+6]<-unique(SYNT$blockid[tmprows,n])
                    ## blockori
                    blocks[k,blockcol+7]<-unique(SYNT$blockori[tmprows,n])
                    ## premask
                    blocks[k,blockcol+8]<-unique(SYNT$premask[tmprows,n])
                }
            }
        }

        colnames(blocks)<-c("start","end","markerS","markerE","car",paste0(rep(c("type","elemS","elemE","node","nodeori","subnode","blockid","blockori","premask"),nhier2-1),rep(1:(nhier2-1),each=9)))

        blocks<-as.data.frame(blocks,stringsAsFactors=FALSE)
        ## -> all columns are characters, only change first five
        blocks$start<-as.numeric(blocks$start)
        blocks$end<-as.numeric(blocks$end)
        class(blocks$markerS)<-class(focalgenome$marker)
        class(blocks$markerE)<-class(focalgenome$marker)
        blocks$car<-as.numeric(blocks$car)

        ## get start and end marker IDs from markers
        blocks$markerS<-markers$marker[scaffrows][blocks$start]
        blocks$markerE<-markers$marker[scaffrows][blocks$end]


        ## ------------------------
        ## summarize rearrangements
        ## ------------------------

        ## NM1
        ## -------
        Re<-SYNT$NM1[scaffrows,]
        tmp<-removeZeros(Re,scaffrows,
                         rownames(SYNT$NM1)[scaffrows],testscaff,1)
        blockNM1<-matrix(NA,ncol=ncol(tmp),nrow=nrow(blocks))

        if(ncol(tmp)>0){
            for(k in 1:nrow(blocks)){
                for(l in 1:ncol(tmp)){
                    val<-unique(tmp[(blocks[k,1]):(blocks[k,2]),l])
                    if(length(val)==1){
                        blockNM1[k,l]<-val
                    }else{
                        stop("Unexpected number of tags for TLcarB")
                    }
                }
            }
        }

        ## NM2
        ## -------
        Re<-SYNT$NM2[scaffrows,]
        tmp<-removeZeros(Re,scaffrows,
                         rownames(SYNT$NM2)[scaffrows],testscaff,1)
        blockNM2<-matrix(NA,ncol=ncol(tmp),nrow=nrow(blocks))

        if(ncol(tmp)>0){
            for(k in 1:nrow(blocks)){
                for(l in 1:ncol(tmp)){
                    val<-unique(tmp[(blocks[k,1]):(blocks[k,2]),l])
                    if(length(val)==1){
                        blockNM2[k,l]<-val
                    }else{
                        stop("Unexpected number of tags for TLcarW")
                    }
                }
            }
        }


        ## SM
        ## ---
        Re<-SYNT$SM[scaffrows,]
        tmp<-removeZeros(Re,scaffrows,
                         rownames(SYNT$SM)[scaffrows],testscaff,1)
        blockSM<-matrix(NA,ncol=ncol(tmp),nrow=nrow(blocks))

        if(ncol(tmp)>0){
            for(k in 1:nrow(blocks)){
                for(l in 1:ncol(tmp)){
                    val<-unique(tmp[(blocks[k,1]):(blocks[k,2]),l])
                    if(length(val)==1){
                        blockSM[k,l]<-val
                    }else{
                        stop("Unexpected number of tags for TP")
                    }
                }
            }
        }


        ## IV
        ## ---
        Re<-SYNT$IV[scaffrows,]
        ## the below assumes that inversion tags are always 1,
        ##  or 0 otherwise, no other values allowed
        if(sum(!is.element(unique(as.vector(Re)),c(0,1)))>0){
            stop("inversion tags in 'SYNT$IV' have to be 0 or 1")
        }
        tmp<-removeZeros(Re,scaffrows,
                         rownames(SYNT$IV)[scaffrows],testscaff,1)
        ## block inversions
        blockIV<-matrix(0,ncol=ncol(tmp),nrow=nrow(blocks))
        ## single-marker inversions
        blockIVsm<-matrix(0,ncol=ncol(tmp),nrow=nrow(blocks))

        if(ncol(tmp)>0){
            smi<-which(colSums(tmp)==1)
            for(k in 1:nrow(blocks)){
                for(l in 1:ncol(tmp)){
                    val<-unique(tmp[(blocks[k,1]):(blocks[k,2]),l])
                    if(length(val)==1){
                        if(!is.element(l,smi)){
                            blockIV[k,l]<-val ## multi-marker inversion
                        }else if(is.element(l,smi) & val==1){
                            blockIVsm[k,l]<-1 ## single-marker inversion
                        }
                    }else{
                        if(sum(tmp[(blocks[k,1]):(blocks[k,2]),l])==1){
                            pos<-which(tmp[(blocks[k,1]):
                                               (blocks[k,2]),l]==1)
                            blockIVsm[k,l]<-pos
                            ## save position of single-marker inversion
                            ##  within block
                        }else{
                            stop("Unexpected number of tags for IV")
                        }
                    }
                }
            }
        }
        blockIV<-removeZeros(blockIV,1:nrow(blockIV),
                              paste0("[",1:nrow(blockIV),",]"),testscaff,1)
        blockIVsm<-removeZeros(blockIVsm,1:nrow(blockIVsm),
                                paste0("[",1:nrow(blockIVsm),",]"),testscaff,1)


        BLOCKS[[s]]<-list(blocks=blocks,NM1=blockNM1,
                          NM2=blockNM2,SM=blockSM,
                          IV=blockIV,IVsm=blockIVsm)
    }

    return(BLOCKS)
}

## ------------------------------------------------------------------------
