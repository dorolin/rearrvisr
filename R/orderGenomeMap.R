## ------------------------------------------------------------------------
## order markers in a genome map file by scaffold and map position
## ------------------------------------------------------------------------

#' Order genome map
#'
#' Order a genome map by genome segments and by the position of markers within
#' genome segments
#'
#' @param genomemap Data frame representing the genome map to be ordered,
#'   containing the mandatory columns \code{$marker}, \code{$scaff},
#'   \code{$start}, \code{$end}, and \code{$strand}, and optional further
#'   columns.
#' @param ordnames Character vector with the names of genome segments (i.e.,
#'   chromosomes or scaffolds) to which the \code{genomemap} will be sorted. The
#'   IDs in the column \code{$scaff} of \code{genomemap} will be matched to the
#'   names in \code{ordnames}, and sorted according to their appearance in
#'   \code{ordnames}.
#' @param partial Integer of value \code{0} or \code{1}. Indicates whether IDs
#'   in the column \code{$scaff} of \code{genomemap} have to match exactly
#'   (\code{partial = 0}) or partially (\code{partial = 1}) to the the names in
#'   \code{ordnames}.
#' @param ordpfx String that is prefix to the names in \code{ordnames}, allowing
#'   additional matches.
#' @param ordsfx String that is suffix to the names in \code{ordnames}, thereby
#'   restricting matches. Only relevant when \code{partial = 1}.
#' @param sortby String indicating whether genome segments that do not have a
#'   unique match to the names in \code{ordnames} will be sorted by their size
#'   (\code{sortby = "size"}) or by their name (\code{sortby = "name"}).
#'
#' @details
#'
#'   \code{genomemap} must contain the mandatory columns \code{$marker} (a
#'   character or integer vector that gives the IDs of markers), \code{$scaff}
#'   (a character vector that gives the ID of the genome segment of origin of
#'   each marker), \code{$start} and \code{$end} (numeric vectors that specify
#'   the location of each marker on its genome segment), and \code{$strand} (a
#'   vector of \code{"+"} and \code{"-"} characters that indicate the reading
#'   direction of each marker). Additional columns are ignored and may store
#'   custom information.
#'
#'   If \code{partial = 0}, only IDs in the column \code{$scaff} of
#'   \code{genomemap} will be considered that either match exactly to the names
#'   in \code{ordnames}, or that match exactly to the combined string of
#'   \code{ordpfx} and \code{ordnames}. \code{ordsfx} is ignored.
#'
#'   If \code{partial = 1}, all IDs in the column \code{$scaff} of
#'   \code{genomemap} will be considered that either start with the combined
#'   string of \code{ordnames} and \code{ordsfx}, or that start with the
#'   combined string of \code{ordpfx}, \code{ordnames}, and \code{ordsfx}.
#'
#'   If more than one ID in the column \code{$scaff} of \code{genomemap} matches
#'   (partially) to a specified genome segment name, matching genome segments
#'   will be sorted by their size (if \code{sortby = "size"}) or by their name
#'   (if \code{sortby = "name"}). All genome segments without any match will
#'   similarly be sorted by their size or name (as specified by \code{sortby})
#'   and appended to the end of the ordered genome map.
#'
#' @return A data frame containing an ordered version of \code{genomemap}. The
#'   IDs in column \code{$scaff} in \code{genomemap} are sorted according to
#'   their appearance in \code{ordnames}. Markers within each genome segment are
#'   sorted by their map position (i.e., the midpoint between positions given by
#'   the columns \code{$start} and \code{$end} in \code{genomemap}).
#'
#' @examples
#' \dontrun{
#'
#' ## specify genome segment names that should appear at the top of
#' ## the genome map, and sort remaining genome segments by their size:
#' SIM_ord1<-orderGenomeMap(SIM_markers, ordnames = c("2", "3", "X"),
#'                          ordpfx = "chr", partial = 1, sortby = "size")
#' head(unique(SIM_ord1$scaff), n = 20L)
#'
#' ## sort all genome segments by name:
#' SIM_ord2<-orderGenomeMap(SIM_markers, ordnames = "all", sortby = "name")
#' head(unique(SIM_ord2$scaff), n = 20L)
#' ## ordnames = "all" is used as a non-matching dummy name
#'
#' ## only sort map positions, keeping original order of genome segments:
#' SIM_ord3<-orderGenomeMap(SIM_markers, ordnames = unique(SIM_markers$scaff))
#' head(unique(SIM_ord3$scaff), n = 20L)
#' }
#'
#' @export

orderGenomeMap<-function(genomemap,ordnames,partial=0,
                         ordpfx="",ordsfx="",sortby="size"){

    ## check genomemap
    checkInfile(genomemap, myclass="focalgenome", checkorder = FALSE)
    ## check settings
    if(length(ordnames) < 1){
        stop("require at least one focal scaffold in ordnames")
    }
    if(class(ordnames) != "character"){
        stop("class for ordnames needs to be 'character'")
    }
    if(anyDuplicated(ordnames) > 0){
        stop("some focal scaffolds in ordnames are duplicated")
    }
    if(length(partial) != 1){
        stop("partial must be of length 1")
    }
    if(!is.element(partial,c(0,1))){
        stop(paste("partial need to be 0 or 1", partial))
    }
    if(length(ordpfx) != 1){
        stop("ordpfx must be of length 1")
    }
    if(class(ordpfx) != "character"){
        stop("class for ordpfx needs to be 'character'")
    }
    if(length(ordsfx) != 1){
        stop("ordsfx must be of length 1")
    }
    if(class(ordsfx) != "character"){
        stop("class for ordsfx needs to be 'character'")
    }
    if(length(sortby) != 1){
        stop("sortby must be of length 1")
    }
    if(!is.element(sortby,c("name","size"))){
        stop("sortby must be one of \"name\" or \"size\"")
    }

    scaford<-makeScafOrd(genomemap,ordnames,partial,ordpfx,ordsfx,sortby)

    orderedmap<-genomemap[order(scaford,genomemap$scaff,
                                   genomemap$start+
                                       (genomemap$end-genomemap$start)/2),]

    ## ## add order ID to genomemap
    ## genomemap$order<-1:nrow(genomemap)

    return(orderedmap)
}

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## order markers by specified scaffolds
## ------------------------------------------------------------------------

makeScafOrd<-function(genomemap,ordnames,partial=0,ordpfx="",ordsfx="",
                      sortby="size"){

    scaford<-numeric(nrow(genomemap))
    for(i in 1:length(ordnames)){
        if(partial==0){
            tmp1<-which(grepl(paste0("^",ordnames[i],"$"),genomemap$scaff,
                              perl=TRUE)==TRUE)
            tmp2<-which(grepl(paste0("^",ordpfx,ordnames[i],"$"),
                              genomemap$scaff,perl=TRUE)==TRUE)
            ## order exact and ordpfx+exact matches
            if(sortby=="size"){
                ## order by scaffold size
                if(length(tmp1)>=length(tmp2)){
                    scaford[tmp1]<-i
                    scaford[tmp2]<-i+0.5
                }else{
                    scaford[tmp1]<-i+0.5
                    scaford[tmp2]<-i
                }
            }else if(sortby=="name"){
                ## order by name
                nord<-order(c(ordnames[i],paste0(ordpfx,ordnames[i])))
                if(nord[1]<nord[2]){
                    scaford[tmp1]<-i
                    scaford[tmp2]<-i+0.5
                }else{
                    scaford[tmp1]<-i+0.5
                    scaford[tmp2]<-i
                }
            }else{
                stop("Unknown ordering option in sortby")
            }
        }else if(partial==1){
            sca<-genomemap$scaff[grepl(paste0("^",ordnames[i],ordsfx),
                                       genomemap$scaff)]
            sca<-c(sca,genomemap$scaff[grepl(paste0("^",ordpfx,ordnames[i],
                                                    ordsfx),genomemap$scaff)])
            ## order exact and partial matches
            if(sortby=="size"){
                tmp2<-names(sort(table(sca),decreasing=TRUE))
                if(length(tmp2)>0){
                    for(j in 1:length(tmp2)){
                        tmp<-which(grepl(paste0("^",tmp2[j],"$"),
                                         genomemap$scaff,perl=TRUE)==TRUE)
                        scaford[tmp]<-i + (j-1)/length(tmp2)
                    }
                }
            }else if(sortby=="name"){
                tmp2<-sort(names(table(sca)))
                if(length(tmp2)>0){
                    for(j in 1:length(tmp2)){
                        tmp<-which(grepl(paste0("^",tmp2[j],"$"),
                                         genomemap$scaff,perl=TRUE)==TRUE)
                        scaford[tmp]<-i + (j-1)/length(tmp2)
                    }
                }
            }else{
                stop("Unknown ordering option in sortby")
            }
        }else{
            stop("Unknown ordering option in partial")
        }
    }
    ## sort remaining scaffolds
    sca<-genomemap$scaff[which(scaford==0)]
    if(length(sca)>0){
        if(sortby=="size"){
            tmp2<-names(sort(table(sca),decreasing=TRUE))
            for(j in 1:length(tmp2)){
                tmp<-which(grepl(paste0("^",tmp2[j],"$"),genomemap$scaff,
                                 perl=TRUE)==TRUE)
                scaford[tmp]<-i + j
            }
        }else if(sortby=="name"){
            tmp2<-sort(names(table(sca)))
            for(j in 1:length(tmp2)){
                tmp<-which(grepl(paste0("^",tmp2[j],"$"),genomemap$scaff,
                                 perl=TRUE)==TRUE)
                scaford[tmp]<-i + j
            }
        }else{
            stop("Unknown ordering option in sortby")
        }
    }

    return(scaford)
}

## ------------------------------------------------------------------------
