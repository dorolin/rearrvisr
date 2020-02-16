## ------------------------------------------------------------------------
## make summary table for number and type of rearrangement events and
##  number of breakpoints for each focal genome segment
## ------------------------------------------------------------------------

#' Summarize Rearrangements
#'
#' For each focal genome segment, summarize the number and type of rearrangement
#' events and the number of breakpoints
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
#' @param remWgt A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Needs to match the value for \code{remWgt} used in the
#'   \code{\link{computeRearrs}} function.
#' @param remThld A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Controls whether breakpoints for components of rearrangements
#'   that are less parsimonious to have changed position relative to the
#'   alternative components will be output. To output all breakpoints,
#'   \code{remThld} needs to be smaller than \code{remWgt} used in the
#'   \code{\link{computeRearrs}} function.
#'
#' @details
#'
#'   Only rearrangements that have components tagged with values larger than
#'   \code{remWgt} will be counted. For proper functioning, \code{remWgt} should
#'   correspond to the value that has been used to generate \code{SYNT}. The
#'   number of nonsyntenic moves is computed as the maximum of class I and class
#'   II nonsyntenic moves per focal genome segment.
#'
#'   The number of breakpoints is computed based on the
#'   \code{\link{getBreakpoints}} function. To include the breakpoint of origin
#'   for nonsyntenic and syntenic moves in the estimate, \code{remThld} needs to
#'   be set to zero (which is the default). Note that this may nevertheless
#'   underestimate the number of breakpoints as the location of origin is not
#'   determined for all rearrangements. When the input is a filtered version of
#'   \code{SYNT} (i.e., filtered with the \code{\link{filterRearrs}} function),
#'   the number of breakpoints may be overestimated. This can be prevented by
#'   increasing the value of \code{remThld} to match the value of \code{remWgt}.
#'   However, this may then underestimate the number of breakpoints as some
#'   breakpoints of origin will not be counted. Breakpoints that fall on
#'   identical positions are only counted once.
#'   
#'
#' @return A matrix with the number of identified nonsyntenic moves, syntenic
#'   moves, inversions, and breakpoints in columns, for the set of focal genome
#'   segments in \code{ordfocal} in rows.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{filterRearrs}},
#'   \code{\link{getBreakpoints}}.
#'
#' @examples
#' SYNT <- computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' summarizeRearrs(SYNT, TOY24_focalgenome, TOY24_compgenome, c("1","2","3"))
#'
#' @export


summarizeRearrs<-function(SYNT, focalgenome, compgenome, ordfocal,
                          remWgt = 0.05, remThld = 0){

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
    if(length(remWgt)!=1){
        stop("remWgt needs to be a single numeric value within [0,0.5)")
    }
    if(remWgt<0 | remWgt>=0.5){
        stop("remWgt needs to be a single numeric value within [0,0.5)")
    }
    if (length(remThld) != 1) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if (remThld < 0 | remThld >= 0.5) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
    }


    ## -------------------------------------------
    ## summarize rearrangements and breakpoints
    ## -------------------------------------------

    ## compute and simplify BLOCKS
    ## ---------------------------
    myblocks<-summarizeBlocks(SYNT,focalgenome,compgenome,ordfocal)
    myblocks<-simplifyBlockTags(myblocks,remThld=remWgt)

    ## compute breakpoints
    ## -------------------
    mybrpts<-getBreakpoints(SYNT,focalgenome,ordfocal,remThld=remThld)

    ## compute stats
    ## -------------
    rearrsum<-matrix(0,nrow=length(ordfocal),ncol=4)
    colnames(rearrsum)<-c("nonsyntenic","syntenic",
                          "inversions","breakpoints")
    rownames(rearrsum)<-ordfocal


    for(s in 1:length(ordfocal)){
        ## rearrangement events
        ## --------------------
        ## TLBS should give min number per chromosome
        allBS<-sum(apply(myblocks[[s]]$TLBS,2,function(x) max(x)>remWgt))
        allWS<-sum(apply(myblocks[[s]]$TLWS,2,function(x) max(x)>remWgt))
        rearrsum[s,1]<-max(allBS,allWS)
        ## special counting for TLWC (potential underestimate for P-nodes)
        fullWC<-sum(apply(myblocks[[s]]$TLWC,2,function(x) max(x)>=1-remWgt))
        halfWC<-sum(apply(myblocks[[s]]$TLWC,2,function(x) max(x)==0.5))/2
        rearrsum[s,2]<-fullWC+ceiling(halfWC)
        ## IV and IVsm
        allIV<-sum(apply(myblocks[[s]]$IV,2,function(x) max(x)>remWgt))
        allIVsm<-sum(apply(myblocks[[s]]$IVsm,2,function(x) max(x)>remWgt))
        rearrsum[s,3]<-allIV+allIVsm
        ## breakpoints
        ## -----------
        tmp<-mybrpts[[s]]
        if(is.null(tmp)){
            rearrsum[s,4]<-0
        }else{
            rearrsum[s,4]<-nrow(tmp)
        }
    }

    ## >>>> note that for TLWS, and for TLWC for P-nodes, tags for
    ##      one event are in the same column (e.g., 0.05 and 0.95,
    ##      but also 0.5 tags as they are separated by gaps anyway);
    ##      for TLWC for Q-nodes, tags for same event can be in separate
    ##      columns (e.g., two columns per transposition pair, one with
    ##      0.05 and one with 0.95 tags, or two columns with 0.5 tags)
    ## >>>> this can lead to underestimation of TLWC's with P-nodes
    ##      with the code above (will require >1 0.5'er event)

    return(rearrsum)
}

## ------------------------------------------------------------------------
