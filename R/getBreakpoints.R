## ------------------------------------------------------------------------
## extract breakpoint coordinates for focal genome segments
## ------------------------------------------------------------------------

#' Get Breakpoints
#'
#' Extract breakpoint coordinates for focal genome segments
#'
#' @param SYNT A list of matrices that store data on different classes of
#'   rearrangements and additional information. \code{SYNT} must have been
#'   generated with the \code{\link{computeRearrs}} function (optionally
#'   filtered with the \code{\link{filterRearrs}} function).
#' @param focalgenome Data frame representing the focal genome, containing the
#'   mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
#'   \code{$end}, and \code{$strand}, and optional further columns. Markers need
#'   to be ordered by their map position.
#' @param ordfocal Character vector with the IDs of the focal genome segments
#'   that will be summarized. Have to match (a subset of) IDs in
#'   \code{focalgenome$scaff}. The default \code{ordfocal = NULL} extracts
#'   breakpoint coordinates for all focal genome segments in
#'   \code{focalgenome$scaff}.
#' @param remThld A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Controls whether breakpoints for components of rearrangements
#'   that are less parsimonious to have changed position relative to the
#'   alternative components will be output. To output all breakpoints,
#'   \code{remThld} needs to be smaller than \code{remWgt} used in the
#'   \code{\link{computeRearrs}} function.
#'
#' @details
#'
#'   Parameters \code{SYNT} and \code{focalgenome} need to be specified.
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
#' @return A list with breakpoint coordinates for the set of focal genome
#'   segments in \code{ordfocal} that have orthologous markers in \code{SYNT}.
#'
#'   If no breakpoints exist, the list element for the focal genome segment is
#'   \code{NULL}, otherwise it is a data frame with breakpoint coordinates in
#'   rows. Columns \code{$bptmin} and \code{$bptmax} give the minimum and
#'   maximum possible coordinates of a breakpoint as the end and start positions
#'   of the two orthologs in \code{SYNT} adjacent to a rearrangement boundary,
#'   and are obtained from \code{focalgenome$end} and \code{focalgenome$start}.
#'   Column \code{$bptmid} gives the breakpoint coordinate as midpoint between
#'   \code{$bptmin} and \code{$bptmax}. Column \code{$maxtagval} gives the
#'   maximum tag value of rearrangements sharing the same breakpoint coordinate.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{filterRearrs}},
#'   \code{\link{genomeImagePlot}}.
#'
#' @examples
#' SYNT<-computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' getBreakpoints(SYNT, TOY24_focalgenome, c("1", "2", "3"))
#'
#' @export


getBreakpoints<-function(SYNT,focalgenome,ordfocal = NULL,remThld = 0){

    ## checks
    ## -------------------------------------------

    ## ERRORS
    ## check SYNT
    checkInfile(SYNT, myclass="SYNT")

    ## check focalgenome
    checkInfile(focalgenome, myclass="focalgenome", checkorder = TRUE)

    ## define ordfocal if not supplied
    if(is.null(ordfocal)){
        ordfocal<-unique(focalgenome$scaff)
    }

    ## check settings
    if(length(ordfocal) < 1){
        stop("require at least one focal segment in ordfocal")
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
    if (length(remThld) != 1) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if (remThld < 0 | remThld >= 0.5) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
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
    ## extract breakpoints for each segment
    ## -------------------------------------------

    scafset<-intersect(ordfocal,unique(markers$scaff))
    if(length(scafset)==0){
        stop("none of the genome segments in 'ordfocal' has markers contained in 'SYNT'")
    }

    brptsL<-vector("list",length(scafset))
    names(brptsL)<-scafset

    for(s in 1:length(scafset)){

        mpos<-which(markers$scaff==scafset[s])

        ## breakpoints (combined for all rearrangement classes)
        tmpS<-cbind(SYNT$IVbS[mpos,,drop=FALSE],SYNT$SMbS[mpos,,drop=FALSE],
                    SYNT$NM2bS[mpos,,drop=FALSE],SYNT$NM1bS[mpos,,drop=FALSE])
        tmpE<-cbind(SYNT$IVbE[mpos,,drop=FALSE],SYNT$SMbE[mpos,,drop=FALSE],
                    SYNT$NM2bE[mpos,,drop=FALSE],SYNT$NM1bE[mpos,,drop=FALSE])
        if(sum(tmpS)>remThld | sum(tmpE)>remThld){
            brpts<-getBreakpnts2BP(tmpS,tmpE,markers$start[mpos],
                                   markers$end[mpos],remThld)
            if(nrow(brpts)>0){
                brptsL[[s]]<-brpts
            } ## otherwise list element remains NULL
        } ## otherwise list element remains NULL
    }

    return(brptsL)
}

## ------------------------------------------------------------------------

