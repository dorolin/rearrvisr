## ------------------------------------------------------------------------
## filter rearrangements and retain only those that comprise a certain
##  number of markers
## ------------------------------------------------------------------------

#' Filter Rearrangements
#'
#' Remove rearrangements that comprise less than a minimum or more than a
#' maximum number of markers
#'
#' @param SYNT A list of matrices that store data on different classes of
#'   rearrangements and additional information. \code{SYNT} must have been
#'   generated with the \code{\link{computeRearrs}} function.
#' @param focalgenome Data frame representing the focal genome, containing the
#'   mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
#'   \code{$end}, and \code{$strand}, and optional further columns. Markers need
#'   to be ordered by their map position.
#' @param filterMin A numerical vector of the form \code{c(nm1, nm2, sm,
#'   iv)} that specifies the minimum number of markers a rearrangement has to
#'   comprise to be retained. \code{nm1} is the minimum number of markers in
#'   \code{SYNT$NM1}, \code{nm2} is the minimum number of markers in
#'   \code{SYNT$NM2}, \code{sm} is the minimum number of markers in
#'   \code{SYNT$SM}, and \code{iv} is the minimum number of markers in
#'   \code{SYNT$IV}.
#' @param filterMax A numerical vector of the form \code{c(nm1, nm2, sm,
#'   iv)} that specifies the maximum number of markers a rearrangement is
#'   allowed to comprise to be retained. \code{nm1} is the maximum number of
#'   markers in \code{SYNT$NM1}, \code{nm2} is the maximum number of markers
#'   in \code{SYNT$NM2}, \code{sm} is the maximum number of markers in
#'   \code{SYNT$SM}, and \code{iv} is the maximum number of markers in
#'   \code{SYNT$IV}.
#'
#' @details
#'
#'   Parameters \code{SYNT} and \code{focalgenome} need to be
#'   specified.
#'
#'   \code{focalgenome} must contain the column \code{$marker}, a vector of
#'   either characters or integers with unique ortholog IDs that can be matched
#'   to the values in the rownames of \code{SYNT}. Values can be \code{NA} for
#'   markers that have no ortholog. \code{$scaff} must be a character vector
#'   giving the name of the focal genome segment (e.g., chromosome or scaffold)
#'   of origin of each marker. \code{$start} and \code{$end} must be numeric
#'   vectors giving the location of each marker on its focal genome segment.
#'   \code{$strand} must be a vector of \code{"+"} and \code{"-"} characters
#'   giving the reading direction of each marker. Additional columns are ignored
#'   and may store custom information, such as marker names. Markers need to be
#'   ordered by their map position within each focal genome segment, for example
#'   by running the \code{\link{orderGenomeMap}} function. \code{focalgenome}
#'   may contain additional rows that were absent when running the
#'   \code{\link{computeRearrs}} function. However, all markers present in
#'   \code{SYNT} need to be contained in \code{focalgenome}, with the subset of
#'   shared markers being in the same order.
#'
#'   Rearrangements are stored in \code{SYNT} and include the following
#'   rearrangement classes: NM1 are class I nonsyntenic moves; NM2 are class II
#'   nonsyntenic moves; SM are syntenic moves; IV are inversions.
#'
#' @return A filtered version of \code{SYNT}. An additional list element
#'   \code{$filter} is created that specifies the applied filter.
#'
#'   Note that for rearrangements that have more than one component, only the
#'   component that falls in the specified filter range is removed. This may
#'   result in an overestimation of the number of breakpoints when a filtered
#'   version of \code{SYNT} is used as input for the
#'   \code{\link{summarizeRearrs}} function.
#'   
#' @seealso \code{\link{computeRearrs}}, \code{\link{genomeImagePlot}},
#'   \code{\link{summarizeRearrs}}.
#'
#' @examples
#' SYNT <- computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' ## only retain inversions comprising at least two markers
#' SYNT_filt<-filterRearrs(SYNT, TOY24_focalgenome, filterMin = c(0, 0, 0, 2))
#'
#' @export


filterRearrs<-function(SYNT,focalgenome,filterMin=c(NA,NA,NA,NA),
                       filterMax=c(NA,NA,NA,NA)){

    concatFlanks<-FALSE
    ## FALSE: each fragment of a flank stored in the same
    ##        column is considered separately
    ## NOTE: when changing this to TRUE the test whether SYNT was
    ##       already filtered needs to take this into account
    ##       (previous or current filtering with concatFlanks==TRUE),
    ##       as two rounds of filtering will lead to wrong results


    ## checks
    ## -------------------------------------------

    ## ERRORS
    ## check SYNT
    checkInfile(SYNT, myclass="SYNT")

    ## check focalgenome
    checkInfile(focalgenome, myclass="focalgenome", checkorder = TRUE)
    ## check settings
    if(length(filterMin) != 4 | !is.vector(filterMin) |
       !(is.numeric(filterMin) | sum(is.na(filterMin))==4)){
        stop("filterMin must be a numeric vector of length four")
    }
    if(length(filterMax) != 4 | !is.vector(filterMax) |
       !(is.numeric(filterMax) | sum(is.na(filterMax))==4)){
        stop("filterMax must be a numeric vector of length four")
    }

    ## WARNINGS
    ## check that SYNT wasn't already filtered
    if(exists('filter',where=SYNT)){
        filterMin<-pmax(filterMin,as.numeric(SYNT$filter[,1]),
                        na.rm=TRUE)
        filterMax<-pmin(filterMax,as.numeric(SYNT$filter[,2]),
                        na.rm=TRUE)
        warning(paste0("SYNT was already filtered.\n    Using adjusted filters:\n    filterMin: ",paste0(filterMin,collapse=" "),"; filterMax: ",paste0(filterMax,collapse=" ")),immediate.=TRUE)
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


    ## -------------------------------------------
    ## filter rearrangements
    ## -------------------------------------------

    SYNT_filt<-SYNT

    ## NM1
    ## ----
    if(isTRUE(ncol(SYNT$NM1)>0) &
       (!is.na(filterMin[1]) | !is.na(filterMax[1]))){

        fmin<-filterMin[1]
        fmin[is.na(fmin)]<-0
        fmax<-filterMax[1]
        fmax[is.na(fmax)]<-Inf

        for(l in 1:ncol(SYNT$NM1)){
            tmp<-SYNT$NM1[,l,drop=FALSE]
            if(colSums(tmp)>0){
                GapFlank<-getGapsFlanks(tmp)
                ## returns $Gaps and $Flanks
                if(concatFlanks==TRUE){
                    flanksize<-rep(sum(GapFlank$Flanks[[1]][1,]),
                                   ncol(GapFlank$Flanks[[1]]))
                }else{
                    flanksize<-GapFlank$Flanks[[1]][1,]
                }
                ## modify values in SYNT_filt for small rearrangements
                for(j in 1:ncol(GapFlank$Flanks[[1]])){
                    if(flanksize[j]<fmin | flanksize[j]>fmax){
                        SYNT_filt$NM1[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j]),l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                        SYNT_filt$NM1bS[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j]),l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                        SYNT_filt$NM1bE[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j]),l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                    }
                }
            }
        }
    }


    myscafs<-unique(markers$scaff)

    for(s in 1:length(myscafs)){
        myset<-which(markers$scaff==myscafs[s])


        ## NM2
        ## ----

        if(isTRUE(ncol(SYNT$NM2)>0) &
           (!is.na(filterMin[2]) | !is.na(filterMax[2]))){

            fmin<-filterMin[2]
            fmin[is.na(fmin)]<-0
            fmax<-filterMax[2]
            fmax[is.na(fmax)]<-Inf

            testset<-SYNT$NM2[myset,,drop=FALSE]

            for(l in 1:ncol(testset)){
                tmp<-testset[,l,drop=FALSE]
                if(colSums(tmp)>0){
                    GapFlank<-getGapsFlanks(tmp)
                    ## returns $Gaps and $Flanks
                    if(concatFlanks==TRUE){
                        flanksize<-rep(sum(GapFlank$Flanks[[1]][1,]),
                                       ncol(GapFlank$Flanks[[1]]))
                    }else{
                        flanksize<-GapFlank$Flanks[[1]][1,]
                    }
                    ## modify values in SYNT_filt for small rearrangements
                    for(j in 1:ncol(GapFlank$Flanks[[1]])){
                        if(flanksize[j]<fmin | flanksize[j]>fmax){
                            SYNT_filt$NM2[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$NM2bS[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$NM2bE[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                        }
                    }
                }
            }
        }

        ## SM
        ## ----

        if(isTRUE(ncol(SYNT$SM)>0) &
           (!is.na(filterMin[3]) | !is.na(filterMax[3]))){

            fmin<-filterMin[3]
            fmin[is.na(fmin)]<-0
            fmax<-filterMax[3]
            fmax[is.na(fmax)]<-Inf

            testset<-SYNT$SM[myset,,drop=FALSE]

            for(l in 1:ncol(testset)){
                tmp<-testset[,l,drop=FALSE]
                if(colSums(tmp)>0){
                    GapFlank<-getGapsFlanks(tmp)
                    ## returns $Gaps and $Flanks
                    if(concatFlanks==TRUE){
                        flanksize<-rep(sum(GapFlank$Flanks[[1]][1,]),
                                       ncol(GapFlank$Flanks[[1]]))
                    }else{
                        flanksize<-GapFlank$Flanks[[1]][1,]
                    }
                    ## modify values in SYNT_filt for small rearrangements
                    for(j in 1:ncol(GapFlank$Flanks[[1]])){
                        if(flanksize[j]<fmin | flanksize[j]>fmax){
                            SYNT_filt$SM[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$SMbS[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$SMbE[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                        }
                    }
                }
            }
        }

        ## IV
        ## ----

        if(isTRUE(ncol(SYNT$IV)>0) &
           (!is.na(filterMin[4]) | !is.na(filterMax[4]))){

            fmin<-filterMin[4]
            fmin[is.na(fmin)]<-0
            fmax<-filterMax[4]
            fmax[is.na(fmax)]<-Inf

            testset<-SYNT$IV[myset,,drop=FALSE]

            for(l in 1:ncol(testset)){
                tmp<-testset[,l,drop=FALSE]
                if(colSums(tmp)>0){
                    GapFlank<-getGapsFlanks(tmp)
                    ## returns $Gaps and $Flanks
                    if(concatFlanks==TRUE){
                        flanksize<-rep(sum(GapFlank$Flanks[[1]][1,]),
                                       ncol(GapFlank$Flanks[[1]]))
                    }else{
                        flanksize<-GapFlank$Flanks[[1]][1,]
                    }
                    ## modify values in SYNT_filt for small rearrangements
                    for(j in 1:ncol(GapFlank$Flanks[[1]])){
                        if(flanksize[j]<fmin | flanksize[j]>fmax){
                            SYNT_filt$IV[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$IVbS[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                            SYNT_filt$IVbE[myset[(GapFlank$Flanks[[1]][2,j]):(GapFlank$Flanks[[1]][3,j])],l]<-rep(0,GapFlank$Flanks[[1]][1,j])
                        }
                    }
                }
            }
        }
    }

    filtspecs<-cbind(filterMin,filterMax)
    rownames(filtspecs)<-c("NM1","NM2","SM","IV")
    SYNT_filt$filter<-filtspecs

    return(SYNT_filt)
}

## ------------------------------------------------------------------------

