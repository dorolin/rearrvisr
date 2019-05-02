## ------------------------------------------------------------------------
## main wrapper to compute rearrangements
## ------------------------------------------------------------------------

#' Compute Rearrangements
#'
#' Detect and classify rearrangements along a focal genome relative to an
#' ancestral genome reconstruction or an extant genome
#'
#' @param focalgenome Data frame representing the focal genome, containing the
#'   mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
#'   \code{$end}, and \code{$strand}, and optional further columns. Markers need
#'   to be ordered by their map position.
#' @param compgenome Data frame representing the compared genome (e.g., an
#'   ancestral genome reconstruction, or an extant genome), with the first three
#'   columns \code{$marker}, \code{$orientation}, and \code{$car}, followed by
#'   columns that alternate between node type and node element. Markers need to
#'   be ordered by their node elements.
#' @param doubled Logical. If \code{TRUE}, markers in the ancestral genome
#'   reconstruction contain information about their orientation.
#' @param remWgt A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Sets the tagging weight for the component of a rearrangement
#'   that is less parsimonious to have changed position relative to the
#'   alternative component to \code{remWgt}, and that of the alternative
#'   component to \code{1 - remWgt}.
#' @param splitnodes Logical. Split nodes into subnodes according to
#'   rearrangements that occurred one step further up the hierarchy during the
#'   rearrangement detection algorithm. \code{splitnodes = TRUE} prevents that
#'   the same rearrangement receives tags across multiple levels of the
#'   \emph{PQ-tree} hierarchy.
## @param chromLev Integer of value \code{0} or \code{1}. Indicates whether the
##   focal genome has been assembled to chromosome-level (\code{chromLev = 1})
##   or not (\code{chromLev = 0}). \code{chromLev = 1} is an experimental
##   setting and should be used with caution.
## @param bestHitOpt Integer of value \code{1}, \code{2}, or \code{3}. Sets the
##   strictness with which focal segment - CAR fragment associations are
##   identified as \emph{best hits}. Such focal segment - CAR fragment
##   \emph{best hits} will not be tagged as rearranged. The least strict level
##   is \code{bestHitOpt = 1}. \code{bestHitOpt = 2} and \code{bestHitOpt = 3}
##   are experimental settings and should be used with caution.
#'
#' @details
#'
#'   \code{focalgenome} must contain the column \code{$marker}, a vector of
#'   either characters or integers with unique ortholog IDs that can be matched
#'   to the values in the \code{$marker} column of \code{compgenome}. Values can
#'   be \code{NA} for markers that have no ortholog. \code{$scaff} must be a
#'   character vector giving the name of the focal genome segment (i.e.,
#'   chromosome or scaffold) of origin of each marker. \code{$start} and
#'   \code{$end} must be numeric vectors giving the location of each marker on
#'   its focal genome segment. \code{$strand} must be a vector of \code{"+"} and
#'   \code{"-"} characters giving the reading direction of each marker.
#'   Additional columns are ignored and may store custom information, such as
#'   marker names. Markers need to be ordered by their map position within each
#'   focal genome segment, for example by running the
#'   \code{\link{orderGenomeMap}} function. See Examples below for the
#'   \code{focalgenome} format.
#'
#'   \code{compgenome} must contain the column \code{$marker}, a vector of
#'   either characters or integers with unique ortholog IDs that can be matched
#'   to the values in the \code{$marker} column of \code{focalgenome}.
#'   \code{$orientation} must be a vector of \code{"+"} and \code{"-"}
#'   characters giving the reading direction of each marker in the compared
#'   genome. If \code{doubled = FALSE}, all values should be \code{"+"}.
#'   \code{$car} must be an integer vector giving the location of each marker on
#'   its compared genome segment (i.e., \emph{Contiguous Ancestral Region}, or
#'   CAR), analogous to contiguous sets of genetic markers on a chromosome,
#'   scaffold, or contig. Each CAR is represented by a \emph{PQ-tree} (Booth &
#'   Lueker 1976; Chauve & Tannier 2008). The \emph{PQ} structure of each CAR is
#'   defined by additional columns (at least two) that have to alternate between
#'   character vectors of node type (\code{"P"}, \code{"Q"}, or \code{NA}) in
#'   even columns, and integer vectors of node elements in odd columns (missing
#'   values are permitted past the fifth column). Every set of node type and
#'   node element column reflects the hierarchical structure of each
#'   \emph{PQ-tree}, with the rightmost columns representing the lowest level of
#'   the hierarchy. \emph{P-nodes} contain contiguous markers and/or nodes in
#'   arbitrary order, while \emph{Q-nodes} contain contiguous markers and/or
#'   nodes in fixed order (including their reversal). For additional details on
#'   \emph{PQ-trees} see Booth & Lueker 1976, Chauve & Tannier 2008, or the
#'   package vignette. See Examples below for the \code{compgenome} format.
#'
#'   \code{doubled = TRUE} indicates that orientation information for the
#'   markers in the ancestral genome reconstruction is available. (This is the
#'   case for example when the genome was reconstructed with the software ANGES,
#'   Jones \emph{et al.} 2012, using the option \code{markers_doubled 1}.)
#'   Orientation information facilitates detecting and classifying
#'   rearrangements as inversions or translocations, and can help determining
#'   whether \emph{PQ-tree} nodes are aligned to the focal genome in ascending
#'   (i.e., standard) or descending (i.e., inverted) direction.
#'
#'   \code{remWgt} provides the tagging weight for rearrangements consisting of
#'   alternative sets of markers, either of which may have caused an apparent
#'   translocation (e.g., a set of markers may have translocated upstream, or
#'   alternatively another set of markers may have translocated downstream). The
#'   set of markers that is more parsimonious to have changed position relative
#'   to the other set receives tag values equal \code{1 - remWgt}, while the
#'   alternative set of markers receives tag values equal \code{remWgt}. Setting
#'   this argument to non-default may require adjusting the \code{remThld}
#'   argument in the \code{genomeImagePlot} and \code{renomeRearrPlot} functions
#'   accordingly.
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
#'   Jones, B. R. \emph{et al.} (2012). ANGES: reconstructing ANcestral GEnomeS
#'   maps. \emph{Bioinformatics}, \strong{28}, 2388--2390. doi:
#'   \href{https://doi.org/10.1093/bioinformatics/bts457}{10.1093/bioinformatics/bts457}
#'
#' @return A list of matrices that store data on different classes of
#'   rearrangements and additional information on the structure of each
#'   \emph{PQ-tree} and its alignment to the focal genome. Markers are in rows,
#'   and the row names of each matrix correspond to the IDs in the
#'   \code{$marker} column of the \code{focalgenome} and \code{compgenome} data
#'   frames. The matrices contain all markers common to \code{focalgenome} and
#'   \code{compgenome}, and are ordered by their position in \code{focalgenome}.
#'
#'   The list elements \code{$TLBS}, \code{$TLWS}, \code{$TLWC}, and \code{$IV}
#'   are numeric matrices that store identified rearrangements. \code{$TLBS}
#'   stores \code{T}rans\code{L}ocations between CARs \code{B}etween focal
#'   \code{S}egments; \code{$TLWS} stores \code{T}rans\code{L}ocations between
#'   CARs \code{W}ithin focal \code{S}egments; \code{$TLWC} stores
#'   \code{T}rans\code{L}ocations within CARs \code{W}ithin focal
#'   \code{S}egments; \code{$IV} stores \code{I}n\code{V}ersions within CARs
#'   within focal segments. See the package vignette for a detailed explanation
#'   of these classes of rearrangements.
#'
#'   Each rearrangement is represented by a separate column. Except for
#'   \code{TLBS}, which are identified across focal segments, columns for
#'   individual focal segments are joined across rows to save space (i.e., for
#'   \code{TLWS}, \code{TLWC}, and \code{IV}, which are identified within focal
#'   segments). To preserve the tabular format, these matrices are filled by
#'   zeros for focal segments with a non-maximal number of rearrangements, if
#'   necessary. If no rearrangements were detected for a certain class, the
#'   matrix has zero columns. Markers that are part of a rearrangement have a
#'   tag value of \code{>0} within their respective column. Tagged markers
#'   within a column are not necessarily consecutive, for example, when a
#'   rearrangement is split into several parts through an insertion of a
#'   different CAR, or when a rearrangement has an upstream and a downstream
#'   component (i.e., when alternative sets of markers may have caused an
#'   apparent translocation). Note that some columns in \code{$TLWS} or
#'   \code{$TLWC} may be duplicated for a particular focal segment due to the
#'   underlying algorithm in \code{\link{computeRearrs}}; although corresponding
#'   to the same rearrangement, these duplicated columns are nevertheless
#'   included for completeness.
#'
#'   For \code{TLBS}, markers part of a translocation have a value of \code{0.5}
#'   if non of the involved CAR fragments is a focal segment - CAR fragment
#'   \emph{best hit}. Otherwise, markers part of the CAR fragment that is
#'   assigned as focal segment - CAR fragment \emph{best hit} have a value of
#'   \code{0}, while markers part of all other non-\emph{best hit} CAR fragments
#'   have a value of \code{1}. For \code{TLWS} and \code{TLWC}, markers part of
#'   a rearrangement with an upstream and a downstream component have a value of
#'   \code{1 - remWgt} (or \code{remWgt}) when they are part of the component
#'   that is more (or less) parsimonious to have changed position; if either
#'   component is equally parsimonious to have changed position, both have a
#'   value of \code{0.5}; all other markers part of a rearrangement have a value
#'   of \code{1}. For \code{IV}, markers part of an inversion have a value of
#'   \code{1}.
#'
#'   The list elements \code{$TLBSbS}, \code{$TLBSbE}, \code{$TLWSbS},
#'   \code{$TLWSbE}, \code{$TLWCbS}, \code{$TLWCbE}, \code{$IVbS}, and
#'   \code{$IVbE} are numeric matrices that tag markers that denote the start
#'   (\code{$*bS}) and end (\code{$*bE}) elements for the four classes of
#'   rearrangements (i.e., the markers adjacent to rearrangement breakpoints).
#'   Each rearrangement is represented by a separate column, but columns for
#'   individual focal segments are joined for all matrices across rows
#'   (including \code{$TLBSbS} and \code{$TLBSbE}) to save space. Tag values
#'   correspond to the ones in \code{$TLBS}, \code{$TLWS}, \code{$TLWC}, and
#'   \code{$IV}.
#'
#'   The list elements \code{$nodeori}, \code{$blockori}, \code{$blockid},
#'   \code{$premask}, and \code{$subnode} are matrices that store information on
#'   the structure of each \emph{PQ-tree}, its alignment to the focal genome,
#'   and internal data. The first column of each matrix corresponds to the CAR
#'   level, and the following columns correspond to the hierarchical structure
#'   of each \emph{PQ-tree}, with information on the lowest level stored in the
#'   last column. \code{$nodeori} is a numeric matrix that stores the alignment
#'   direction of each \emph{Q-node} to the focal genome, with \code{1}
#'   indicating ascending (i.e., standard), and \code{-1} descending (i.e.,
#'   inverted) alignment. \emph{Q-nodes} that have no alignment direction (e.g.,
#'   single-marker nodes) have a value of \code{9}, and \emph{P-nodes} are
#'   \code{NA}. \code{$blockori} is a numeric matrix that stores the orientation
#'   of each synteny block, with \code{1} indicating ascending (i.e., standard),
#'   and \code{-1} descending (i.e., inverted) orientation. Blocks that have no
#'   orientation (e.g., blocks containing a single marker, or a single
#'   \emph{PQ-tree} branch) have a value of \code{9}, and blocks that are part
#'   of \emph{P-nodes} are \code{NA}. \code{$blockid} is a character matrix that
#'   stores the ID of each synteny block within its node. For \emph{Q-nodes},
#'   IDs are consecutive and start at \code{1}, separately for each node and
#'   each hierarchy level, and reflect the order of synteny blocks. Block IDs
#'   with \code{".1"} or \code{".2"} suffixes (in arbitrary order) indicate
#'   blocks that were subject to an additional subdivision step. For
#'   \emph{P-nodes}, IDs are \code{0} unless the node is part of a
#'   rearrangement, in which case IDs indicate different rearrangements, but not
#'   block order. \code{$premask} and \code{$subnode} are numeric matrices that
#'   store internal data used for the alignment and identification of
#'   rearrangements. Integers \code{>0} in \code{$subnode} indicate subdivisions
#'   of the corresponding \emph{PQ-tree} due to translocations. All subdivisions
#'   have been searched separately for rearrangements one step further down the
#'   hierarchy. This is of main relevance when \code{splitnodes = TRUE}.
#'
#'   The returned data can be visualized with the \code{\link{genomeImagePlot}}
#'   function, or summarized and visualized with the
#'   \code{\link{summarizeBlocks}} and \code{\link{genomeRearrPlot}} functions.
#'
#' @seealso \code{\link{genomeImagePlot}}, \code{\link{summarizeBlocks}},
#'   \code{\link{genomeRearrPlot}}; \code{\link{orderGenomeMap}} to order the
#'   \code{focalgenome} data frame; \code{\link{convertPQtree}} or
#'   \code{\link{genome2PQtree}} to generate the \code{compgenome} data frame.
#'
#' @examples
#' computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' \dontrun{
#'
#' ## focalgenome format:
#' TOY24_focalgenome
#'
#' ## compgenome format:
#' TOY24_compgenome
#' }
#'
#' @export


## computeRearrs<-function(focalgenome, compgenome, doubled, remWgt = 0.05,
##                         splitnodes = TRUE, chromLev = 0, bestHitOpt = 1){
## computeRearrs<-function(focalgenome, compgenome, doubled, remWgt = 0.05,
##                         splitnodes = TRUE, chromLev = 0){
computeRearrs<-function(focalgenome, compgenome, doubled, remWgt = 0.05,
                        splitnodes = TRUE){

    ## removed from function arguments
    chromLev <- 0
    bestHitOpt <- 1

    ## checks
    ## -------------------------------------------

    ## ERRORS
    ## check compgenome
    checkInfile(compgenome, myclass="compgenome", checkorder = TRUE)
    ## check focalgenome
    checkInfile(focalgenome, myclass="focalgenome", checkorder = TRUE)
    ## class of marker in tree and marker file must be the same
    ## (either character or integer)
    if(class(focalgenome$marker) != class(compgenome$marker)){
        stop("class for '$marker' in focalgenome and compgenome need to be identical")
    }
    ## check settings
    if(length(doubled)!=1 | !is.logical(doubled)){
        stop("doubled must be 'TRUE' or 'FALSE'")
    }
    if(length(remWgt)!=1){
        stop("remWgt needs to be a single numeric value within [0,0.5)")
    }
    if(remWgt<0 | remWgt>=0.5){
        stop("remWgt needs to be a single numeric value within [0,0.5)")
    }
    if(length(splitnodes)!=1 | !is.logical(splitnodes)){
        stop("splitnodes must be 'TRUE' or 'FALSE'")
    }
    if(length(chromLev)!=1){
        stop("chromLev needs to be a single integer")
    }
    if(!is.element(chromLev,c(0,1))){
        stop("chromLev needs to be 0 or 1")
    }
    if(length(bestHitOpt)!=1){
        stop("bestHitOpt needs to be a single integer")
    }
    if(!is.element(bestHitOpt,c(1,2,3))){
        stop("bestHitOpt needs to be 1, 2, or 3")
    }

    ## WARNINGS
    if(sum(compgenome$orientation=="-")>0 & doubled==FALSE){
        warning("sure with having no orientation information in tree?",immediate.=TRUE)
    }
    if(sum(compgenome$orientation=="-")==0 & doubled==TRUE){
        warning("sure with having correct orientation information in tree?",immediate.=TRUE)
    }
    ## experimental settings
    if(!is.element(chromLev,0)){
        warning(paste0("chromLev=",chromLev," is experimental setting"),immediate.=TRUE)
    }
    if(!is.element(bestHitOpt,1)){
        warning(paste0("bestHitOpt=",bestHitOpt," is experimental setting"),immediate.=TRUE)
    }


    ## initial processing
    ## -------------------------------------------

    ntreecol<-ncol(compgenome)
    ## levels of hierarchy in tree (including CARs)
    nhier<-1+(ntreecol-3)/2

    treepos<-match(focalgenome$marker,compgenome$marker)
    treepos<-treepos[!is.na(treepos)]
    ## exclude markers from tree that are not in focalgenome
    tree<-compgenome[treepos,,drop=FALSE]
    markerpos<-match(tree$marker,focalgenome$marker)
    markers<-focalgenome[markerpos,,drop=FALSE]
    ## exclude markers (incl. NAs) from focalgenome that are not in tree
    if(sum(sum(markers$marker != tree$marker)) > 0){
        stop("markers in compgenome and focalgenome could not correctly be subset")
    }

    ## check original tree for nodes containing only one marker
    ## (compgenome needs to be in original order)
    allsinglemarkers<-getSingles(compgenome)
    singlemarkers<-allsinglemarkers[treepos]

    ## markers have same (1) or opposite orientation (-1) with tree
    orientation<-rep(NA,nrow(markers))
    if(doubled==TRUE){
        orientation[tree$orientation==markers$strand]<-1
        orientation[tree$orientation!=markers$strand]<- -1
        orientation[singlemarkers==1]<-NA
    }

    ##myscafs<-sort(unique(markers$scaff))
    myscafs<-unique(markers$scaff)

    mycars<-sort(unique(tree$car))


    ## set up main storage
    ## -------------------------------------------
    tmp<-matrix(0,ncol=nhier,nrow=nrow(markers))
    rownames(tmp)<-tree$marker

    SYNT<-list(TLBS=matrix(NA,ncol=0,nrow=0),
               TLWS=matrix(NA,ncol=0,nrow=0),
               TLWC=matrix(NA,ncol=0,nrow=0),
               IV=matrix(NA,ncol=0,nrow=0),
               TLBSbS=matrix(NA,ncol=0,nrow=0),
               TLBSbE=matrix(NA,ncol=0,nrow=0),
               TLWSbS=matrix(NA,ncol=0,nrow=0),
               TLWSbE=matrix(NA,ncol=0,nrow=0),
               TLWCbS=matrix(NA,ncol=0,nrow=0),
               TLWCbE=matrix(NA,ncol=0,nrow=0),
               IVbS=matrix(NA,ncol=0,nrow=0),
               IVbE=matrix(NA,ncol=0,nrow=0),
               nodeori=matrix(NA,ncol=nhier,nrow=0),
               blockori=matrix(NA,ncol=nhier,nrow=0),
               blockid=matrix(NA,ncol=nhier,nrow=0),
               premask=matrix(NA,ncol=nhier,nrow=0),
               subnode=tmp
               )


    ## main functions wrapper
    ## -------------------------------------------

    ## identify best hits for scaffold - CAR pairs
    scafcarbest<-getBestHits(markers,tree,myscafs,mycars)

    ## check if CAR - scaffold assignments indicate translocations
    TLL<-tagTLcar2(markers,tree,myscafs,mycars,chromLev=chromLev)
    ## returns $TLbetween, $TLwithin

    ## post-processing to filter out best CAR-scaffold hits and to
    ##  filter out large inserts at CAR level
    SYNT<-filterCars3(SYNT,markers,tree,nhier,myscafs,mycars,TLL,
                      scafcarbest,bestHitOpt=bestHitOpt,remWgt=remWgt)

    ## identify intra-scaffold rearrangements
    SYNT<-tagRearr(SYNT,markers,tree,orientation,nhier,myscafs,
                   splitnodes=splitnodes,remWgt=remWgt)


    colnames(SYNT$TLBS)<-NULL
    colnames(SYNT$TLWS)<-NULL
    colnames(SYNT$TLWC)<-NULL
    colnames(SYNT$IV)<-NULL
    colnames(SYNT$TLBSbS)<-NULL
    colnames(SYNT$TLBSbE)<-NULL
    colnames(SYNT$TLWSbS)<-NULL
    colnames(SYNT$TLWSbE)<-NULL
    colnames(SYNT$TLWCbS)<-NULL
    colnames(SYNT$TLWCbE)<-NULL
    colnames(SYNT$IVbS)<-NULL
    colnames(SYNT$IVbE)<-NULL

    return(SYNT)
}

## ------------------------------------------------------------------------



