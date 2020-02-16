## ------------------------------------------------------------------------
## plot blocks and rearrangements
## ------------------------------------------------------------------------

#' Genome Rearrangements Plot
#'
#' Generate a plot that shows synteny blocks between a focal and a compared
#' genome in columns, and information on their alignment and rearrangements in
#' rows, for a given set of focal genome segments
#'
#' @param BLOCKS A list of lists generated with the
#'   \code{\link{summarizeBlocks}} function. The top-level list elements of
#'   \code{BLOCKS} are focal genome segments, and the lower-level list elements
#'   contain information on synteny blocks and rearrangements for each focal
#'   genome segment.
#' @param compgenome Data frame representing the compared genome (e.g., an
#'   ancestral genome reconstruction, or an extant genome), with the first three
#'   columns \code{$marker}, \code{$orientation}, and \code{$car}, followed by
#'   columns alternating node type and node element. Markers need to be ordered
#'   by their node elements. It must be the same data frame that was used to
#'   generate the list \code{BLOCKS} with the \code{\link{summarizeBlocks}}
#'   function.
#' @param ordfocal Character vector with the IDs of the focal genome segments
#'   that will be plotted. Have to match (a subset of) names of the top-level
#'   list elements of \code{BLOCKS}.
#' @param remstr String that should be removed from the IDs in \code{ordfocal}
#'   to simplify the y-axis labels. Only relevant when \code{yaxlab = NULL}.
#' @param main Title of the plot.
#' @param remThld A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Controls whether components of rearrangements that are less
#'   parsimonious to have changed position relative to the alternative
#'   components will be plotted. To plot all components, \code{remThld} needs to
#'   be smaller than \code{remWgt} used in the \code{\link{computeRearrs}}
#'   function.
#' @param mar A numerical vector of the form \code{c(bottom, left, top, right)}
#'   that specifies the margins on the four sides of the plot. The default
#'   \code{mar = NULL} sets the margins automatically.
#' @param pad A numeric value of \code{0} or greater that sets the amount of
#'   space between all plot margins and the actual plot area.
#' @param y0pad A numeric value of \code{0} or greater that sets the amount of
#'   additional space between the bottom plot margin and the bottom plot area.
#'   Setting this value too small may result in some rearrangements for the
#'   bottom-most focal genome segment to be outside the bottom plot area (and
#'   thus to be invisible in the plot).
#' @param uniqueCarColor Logical. If \code{TRUE}, CARs are uniquely colored
#'   across all focal genome segments. If \code{FALSE}, CARs are colored
#'   separately for each focal genome segment based on the number of markers per
#'   CAR (forces \code{sortColsBySize = TRUE}).
#' @param sortColsBySize Logical. If \code{TRUE}, \code{carColors} and
#'   \code{carTextColors} are assigned based on the number of markers per CAR,
#'   so that the first color is allocated to the largest CAR.
#' @param plotelem A numerical vector of the form \code{c(nor, bid, bor, eid,
#'   rea)} that determines which synteny block information is visualized.
#'   \code{nor} is the alignment orientation of the node, \code{bid} is the
#'   block ID, \code{bor} is the block orientation within its node, \code{eid}
#'   is the element ID within its node, and \code{rea} are rearrangements. The
#'   information is plotted when the value is \code{1}, and omitted when
#'   \code{0}.
#' @param simplifyTags Logical. If \code{TRUE}, duplicated rearrangement tags,
#'   if any, are excluded from the plot. Note that this will work properly only
#'   when \code{remWgt} used in the \code{\link{computeRearrs}} function was set
#'   to a value \code{>0}.
#' @param blockwidth A numeric value that specifies the relative width of
#'   synteny blocks (i.e., columns) in the plot.
#' @param yaxlab Annotations of the y-axis. Must be the same length as
#'   \code{ordfocal}. The default \code{yaxlab = NULL} uses as annotations the
#'   names in \code{ordfocal}, simplified through removal of the string in
#'   \code{remstr}.
#' @param carColors Character vector with the color names used for coloring
#'   CARs. If the number of CARs is greater than \code{length(carColors)},
#'   remaining CARs are colored in grayscale. \code{carColors} must have the
#'   same length as \code{carTextColors}.
#' @param carTextColors Character vector with the color names used for
#'   coloring CAR IDs. \code{carTextColors} must have the same length as
#'   \code{carColors}.
#' @param cex.main Numerical value that specifies the magnification of the main
#'   title.
#' @param cex.axis Numerical value that specifies the magnification of the axis
#'   annotation.
#' @param cex.text Numerical value that specifies the magnification of text
#'   within the plot area.
#' @param font.main Font for the title of the plot. \code{1} corresponds to
#'   plain text, \code{2} to bold face, \code{3} to italic, and \code{4} to bold
#'   italic.
#' @param makepdf Logical. Save plot as PDF. See \code{filename} and
#'   \code{colormodel}.
#' @param newdev Logical. Opens a new default graphics device (but not
#'   "RStudioGD") via \code{dev.new}. Only relevant when \code{makepdf = FALSE}.
#' @param filename Character string that gives the name of the PDF file when
#'   \code{makepdf = TRUE}.
#' @param colormodel Character string that gives the color model for the PDF
#'   when \code{makepdf = TRUE}. Allowed values are \code{"srgb"},
#'   \code{"cmyk"}, \code{"gray"}, or \code{"grey"}.
#'
#' @details
#'
#' Parameters \code{BLOCKS}, \code{compgenome}, and \code{ordfocal} need to be
#' specified, all other parameters have default or automatic settings.
#'
#' When \code{makepdf = TRUE} or \code{newdev = TRUE}, the width and height of
#' the graphic will be set automatically. The dimensions are determined in
#' inches, thus \code{makepdf = FALSE} and \code{newdev = TRUE} will produce an
#' error or not work correctly when the default units of the default graphics
#' device are not inches (such as \code{\link[grDevices]{bmp}},
#' \code{\link[grDevices]{jpeg}}, \code{\link[grDevices]{png}}, or
#' \code{\link[grDevices]{tiff}}). This can be avoided by setting the default
#' graphics device to one that has inches as default units. Setting both
#' \code{makepdf = FALSE} and \code{newdev = FALSE} will allow to specify
#' alternative, user-defined dimensions of the graphic. See examples below.
#'
#' Colors are assigned to CARs by size, unless \code{sortColsBySize = FALSE}.
#' When \code{carColors = NULL}, 47 easily distinguishable default colors are
#' used for coloring CARs. The first 14 colors are color blindness friendly and
#' were obtained from
#' \href{http://mkweb.bcgsc.ca/biovis2012}{mkweb.bcgsc.ca/biovis2012}.
#' \code{carTextColors} are either black or white dependent on the hue of the
#' default \code{carColors}.
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
#' @return A plot to the default graphics device (but not "RStudioGD") or a PDF
#'   file.
#'
#'   The plot visualizes the data contained in \code{BLOCKS} for each focal
#'   genome segment in \code{ordfocal}, arranged along the y-axis. Each synteny
#'   block is represented by a column (corresponding to rows in \code{BLOCKS}),
#'   and information on each block is visualized in rows (corresponding to
#'   columns in \code{BLOCKS}). Note that separate blocks are also generated
#'   when the hierarchical structure of the underlying \emph{PQ-tree} changes,
#'   therefore not all column boundaries are caused by a rearrangement. Rows
#'   provide information on the structure of each \emph{PQ-tree} and its
#'   alignment to the focal genome, and whether blocks are part of different
#'   classes of rearrangements. For details on \emph{PQ-trees} see the
#'   description of the \code{"compgenome"} class in the Details section of the
#'   \code{\link{checkInfile}} function, Booth & Lueker 1976, Chauve & Tannier
#'   2008, or the package vignette.
#'
#'   For each focal genome segment, the top row gives the number of markers
#'   within each block, followed by a row that gives the IDs of the CARs. The
#'   remaining rows are optional and are controlled by the values in
#'   \code{plotelem}.
#'
#'   Up to four rows (depending on the values in the argument \code{plotelem})
#'   describe the structure of each \emph{PQ-tree} and its alignment to the
#'   focal genome, and are stacked for each level of the \emph{PQ-tree}
#'   hierarchy. The first row (\code{nor}) gives the alignment orientation of
#'   the \emph{PQ-tree} node to the focal genome, with white rectangles and
#'   \code{"+"} indicating ascending (i.e., standard), and black rectangles and
#'   \code{"-"} indicating descending (i.e., inverted) alignment. Nodes that
#'   have no alignment direction (e.g., single-marker nodes) are in gray, and
#'   \emph{P-nodes} (which have no fixed ordering and thus no direction) have
#'   gray shaded rectangles. The second row (\code{bid}) gives the block ID. For
#'   \emph{Q-nodes}, IDs are consecutive and start at \code{1}, separately for
#'   each node and each hierarchy level, and reflect the order of synteny
#'   blocks. Identical IDs mean that blocks might be joined, but are split
#'   either by an insertion of another CAR or because of a change in the
#'   hierarchical structure of the underlying \emph{PQ-tree}. Block IDs with
#'   \code{".1"} or \code{".2"} suffixes indicate rearrangements that may arise
#'   by either an inversion or a syntenic move between adjacent blocks, but that
#'   were classified as syntenic move for the sake of parsimony. For
#'   \emph{P-nodes}, the \code{bid} row is empty unless the node is part of a
#'   rearrangement, in which case IDs indicate different rearrangements, but not
#'   block order. The third row (\code{bor}) gives the block orientation within
#'   its node. It has the same color and symbol coding as \code{nor} above. For
#'   example, a \code{"-"} block within a \code{"+"} node indicates either an
#'   inversion or a syntenic move between adjacent blocks. The fourth row
#'   (\code{eid}) gives the range of element IDs for each block within its node
#'   and for its level of hierarchy. These IDs correspond to the node elements
#'   in the odd columns of \code{compgenome} (note that some IDs within blocks
#'   or in-between might be missing when markers in the compared genome are
#'   absent from the focal genome).
#'
#'   The final set of rows (\code{rea}) indicates whether blocks are part of
#'   different classes of rearrangements. Horizontal lines that are at identical
#'   height denote the same rearrangement (potentially disrupted by inserted
#'   CARs). Green are class I nonsyntenic moves (NM1); blue are class II
#'   nonsyntenic moves (NM2); purple are syntenic moves (SM); maroon are
#'   inversions (IV). Inversions that involve only a single marker (i.e.,
#'   markers with switched orientation) are indicated by a short vertical rather
#'   than a horizontal line. Lighter coloration denotes smaller weights for
#'   rearrangement tags in the respective matrices in \code{BLOCKS}. Unless the
#'   argument \code{remThld} is set to a value smaller than that of
#'   \code{remWgt} used in the \code{\link{computeRearrs}} function, only lines
#'   for blocks that are more parsimonious to have changed position relative to
#'   alternative blocks are plotted. If \code{simplifyTags = FALSE}, all tags
#'   for NM2 and SM will be plotted for completeness, i.e., including those that
#'   are duplicated due to the functioning of the underlying algorithm in
#'   \code{\link{computeRearrs}}.
#'
#' @seealso \code{\link{checkInfile}}, \code{\link{computeRearrs}},
#'   \code{\link{summarizeBlocks}}, \code{\link{genomeImagePlot}}. For more
#'   information about arguments that are passed to other functions, see
#'   \code{\link[grDevices]{dev.new}}, \code{\link[grDevices]{pdf}},
#'   \code{\link[graphics]{plot}}, \code{\link[graphics]{par}}.
#'
#' @examples
#' SYNT<-computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#' BLOCKS<-summarizeBlocks(SYNT, TOY24_focalgenome, TOY24_compgenome,
#'                         c("1","2","3"))
#'
#' genomeRearrPlot(BLOCKS, TOY24_compgenome, c("1", "2", "3"), main = "TOY24")
#'
#' \dontrun{
#'
#' ## make PDF (automatically determine the width and height of the graphic)
#' genomeRearrPlot(BLOCKS, TOY24_compgenome, c("1", "2", "3"), main = "TOY24",
#'                 makepdf = TRUE, newdev = FALSE, filename = "rearr.pdf")
#'
#' ## make PDF (default dimensions, i.e., square format)
#' pdf(file = "rearr.pdf")
#' genomeRearrPlot(BLOCKS, TOY24_compgenome, c("1", "2", "3"), main = "TOY24",
#'                 makepdf = FALSE, newdev = FALSE)
#' dev.off()
#'
#' ## plot in R Studio window
#' op <- options(device = "RStudioGD")
#' genomeRearrPlot(BLOCKS, TOY24_compgenome, c("1", "2", "3"), main = "TOY24",
#'                 newdev = FALSE)
#' options(op)
#'
#' ## make EPS, and set user-specified dimensions
#' setEPS()
#' postscript("rearr.eps", width=4.5,height=6.0,pointsize=9)
#' genomeRearrPlot(BLOCKS, TOY24_compgenome, c("1", "2"), main = "TOY24",
#'                 pad = 1, y0pad = 1, makepdf = FALSE, newdev = FALSE)
#' dev.off()
#' }
#'
#' @export
#' @importFrom grDevices dev.new dev.off pdf colorRampPalette
#' @importFrom graphics par plot text rect segments axis

genomeRearrPlot<-function(BLOCKS,compgenome,ordfocal,
                          remstr = "",
                          main = "",
                          remThld = 0.05,
                          mar = NULL,
                          pad = 0,
                          y0pad = 5,
                          uniqueCarColor = TRUE,
                          sortColsBySize = TRUE,
                          plotelem = c(1,1,1,1,1),
                          simplifyTags = TRUE,
                          blockwidth = 1,
                          yaxlab = NULL,
                          carColors = NULL,
                          carTextColors = NULL,
                          cex.main = 2,
                          cex.axis = 1.2,
                          cex.text = 1,
                          font.main = 3,
                          makepdf = FALSE,
                          newdev = TRUE,
                          filename = "rearr.pdf",
                          colormodel = "srgb"){

    ## removed from function arguments
    space <- NULL
    mylim <- NULL

    ## checks
    ## -------------------------------------------

    ## ERRORS

    ## check BLOCKS
    checkInfile(BLOCKS, myclass="BLOCKS")
    ## check compgenome
    checkInfile(compgenome, myclass="compgenome", checkorder = FALSE)
    ## check settings
    if(length(ordfocal) < 1){
        stop("require at least one focal segment in ordfocal")
    }
    if(class(ordfocal) != "character"){
        stop("class for focal segments in ordfocal needs to be 'character'")
    }
    if(anyDuplicated(ordfocal) > 0){
        stop("some focal segments in ordfocal are duplicated")
    }
    if(sum(is.na(match(ordfocal,names(BLOCKS)))) > 0){
        stop(paste("focal segment", ordfocal[is.na(match(ordfocal,names(BLOCKS)))][1], "is absent in BLOCKS"))
    }
    if(length(remstr) != 1){
        stop("remstr must be unique string of length one")
    }
    if(class(remstr) != "character"){
        stop("class for remstr needs to be 'character'")
    }
    if(length(remThld)!=1){
        stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if(remThld<0 | remThld>=0.5){
        stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if(length(pad) != 1 | sum(pad < 0) > 0 | !is.numeric(pad)){
        stop("pad needs to be a single positive numeric value")
    }
    if(length(y0pad) != 1 | sum(y0pad < 0) > 0 | !is.numeric(y0pad)){
        stop("y0pad needs to be a single positive numeric value")
    }
    if(length(uniqueCarColor) != 1 | !is.logical(uniqueCarColor)){
        stop("uniqueCarColor must be 'TRUE' or 'FALSE'")
    }
    if(length(sortColsBySize) != 1 | !is.logical(sortColsBySize)){
        stop("sortColsBySize must be 'TRUE' or 'FALSE'")
    }
    if(length(plotelem) != 5 | !is.vector(plotelem) | !is.numeric(plotelem)){
        stop("plotelem must be a numeric vector of length five")
    }
    if(sum(plotelem != 0 & plotelem != 1) > 0){
        stop("values of plotelem need to be 0 or 1")
    }
    if(length(simplifyTags) != 1 | !is.logical(simplifyTags)){
        stop("simplifyTags must be 'TRUE' or 'FALSE'")
    }
    if(length(blockwidth) != 1 | sum(blockwidth < 0) > 0 |
       !is.numeric(blockwidth)){
        stop("blockwidth needs to be a single positive numeric value")
    }
    if(length(makepdf) != 1 | !is.logical(makepdf)){
        stop("makepdf must be 'TRUE' or 'FALSE'")
    }
    if(length(newdev) != 1 | !is.logical(newdev)){
        stop("newdev must be 'TRUE' or 'FALSE'")
    }
    if(!is.element(colormodel,c("srgb","cmyk","gray","grey"))){
        stop("colormodel needs to be 'srgb', 'cmyk', 'gray', or 'grey'")
    }
    if(length(cex.axis) !=1 | sum(cex.axis < 0) > 0 | !is.numeric(cex.axis)){
        stop("cex.axis must be a single positive numeric value")
    }
    ## cex.main, cex.axis, font.main, filename should be checked by
    ##  'par()' or 'pdf()'/'dev.new()'

    ## check automatic settings
    ## heights of each plot row
    if(!is.null(space)){
        if(!is.data.frame(space)){
            stop("space needs to be a data frame")
        }
        if(any(!is.element(c("nmark","carid","ndori","blkid","blkori",
                             "elem","stat","rearr","scaff"),
                           colnames(space)))){
            stop("column names in space need to include 'nmark', 'carid',\n    'ndori', 'blkid', 'blkori', 'elem', 'stat', 'rearr', 'scaff'")
        }
        if(nrow(space) != 1){
            stop("space must have exactly one row")
        }
        if(!is.numeric(as.matrix(space)[1,]) | sum(space[1,] < 0) > 0){
            stop("values in space cannot be negative")
        }
        space$ndori<-space$ndori*plotelem[1]
        space$blkid<-space$blkid*plotelem[2]
        space$blkori<-space$blkori*plotelem[3]
        space$elem<-space$elem*plotelem[4]
        space$stat<-space$stat*plotelem[5]
        space$rearr<-space$rearr*plotelem[5]
    }
    ## simplified y-axis labels
    if(!is.null(yaxlab)){
        if(length(yaxlab) != length(ordfocal) |
           !is.vector(yaxlab)){
            stop("yaxlab must have the same length as ordfocal")
        }
    }
    ## car colors
    if(!is.null(carColors)){
        if(length(carColors) != length(carTextColors) |
           !is.vector(carColors)){
            stop("carColors and carTextColors must have the same length")
        }
    }
    ## car text colors
    if(!is.null(carTextColors)){
        if(length(carColors) != length(carTextColors) |
           !is.vector(carTextColors)){
            stop("carColors and carTextColors must have the same length")
        }
    }
    ## plot margins
    if(!is.null(mar)){
        if(length(mar) != 4 | !is.vector(mar) | !is.numeric(mar)){
            stop("mar must be a numeric vector of length four")
        }
        if(sum(mar < 0) > 0){
            stop("values of mar cannot be negative")
        }
    }
    ## dimensions of the plot
    if(!is.null(mylim)){
        if(!is.list(mylim)){
            stop("mylim needs to be a list")
        }
        if(any(!is.element(c("x","y","height","width"),names(mylim)))){
            stop("objects in mylim need to include $x, $y, $height, $width")
        }
        ## check $x, $y, $height, $width
        if(length(mylim$x) != 2 | !is.vector(mylim$x) |
           !is.numeric(mylim$x)){
            stop("mylim$x must be a numeric vector of length two")
        }
        if(length(mylim$y) != 2 | !is.vector(mylim$y) |
           !is.numeric(mylim$y)){
            stop("mylim$y must be a numeric vector of length two")
        }
        if(length(mylim$height) != 1 | length(mylim$width) != 1 |
           !is.numeric(mylim$height) | !is.numeric(mylim$width) |
           sum(c(mylim$height,mylim$width) < 0) > 0){
            stop("mylim$height and mylim$width must be a positive numeric\n    value of length one")
        }
    }
    ## main label: can't think of any serious errors


    ## initial processing
    ## -------------------------------------------

    ## simplify BLOCK tags for NM2 and SM
    if(simplifyTags==TRUE){
        BLOCKS<-simplifyBlockTags(BLOCKS,remThld)
    }

    ## make automatic settings
    ## heights of each plot row
    if(is.null(space)){
        space<-data.frame(nmark=0.65,
                          carid=1,
                          ndori=0.4*plotelem[1],
                          blkid=1*plotelem[2],
                          blkori=0.4*plotelem[3],
                          elem=0.65*plotelem[4],
                          stat=0.75*plotelem[5],
                          rearr=0.25*plotelem[5],
                          scaff=2.25)
    }
    ## simplified y-axis labels
    if(is.null(yaxlab)){
        yaxlab<-gsub(remstr,"",ordfocal)
    }
    ## plot margins
    if(is.null(mar)){
        mar<-c(0.4,0.6+max(nchar(yaxlab))*0.42*cex.axis,
               2.2*cex.main,0.4)
    }
    ## dimensions of the plot
    pad<-pad+0.7
    ## some rearrangement tags of last scaffold might be invisible
    ##  when not having some padding (because yaxs="i")
    if(is.null(mylim)){
        mylim<-setRearrPlot(BLOCKS,ordfocal,space,blockwidth,
                            mar,remThld,pad,y0pad)
        ## returns $x, $y, $height, $width
    }

    ## make colors for cars
    ## --------------------
    if(is.null(carColors)){
        ## from http://mkweb.bcgsc.ca/biovis2012/color-blindness-palette.png,
        ##  coverted to hex code with grDevices::rgb()
        mycolors<-c("#490092","#920000","#004949","#006DDB","#924900",
                    "#009292","#B66DFF","#DB6D00","#FF6DB6","#6DB6FF",
                    "#24FF24","#FFB6DB","#B6DBFF","#FFFF6D")
        ## many of the following colors are from
        ##  http://r-statistics.co/ggplot2-cheatsheet.html,
        ##  but rearranged and excluded/replaced some, and converted to
        ##  hex code with gplots::col2hex()
        mycolorsB<-c("#0000CD","#FF4500","#6B8E23","#4682B4","#BC8F8F",
                     "#228B22","#F4A460","#FF1493","#1E90FF","#ADFF2F",
                     "#D8BFD8","#F5F5F5","#191970","#CD2626","#2E8B57",
                     "#87CEEB","#BDB76B","#FFDEAD","#6A5ACD","#D02090",
                     "#696969","#FF6347","#EEE8AA","#FFF5EE","#0000FF",
                     "#B03060","#00FF7F","#DDA0DD","#FA8072","#32CD32",
                     "#EEDD82","#FFD700","#F0FFFF")
        carColors<-c(mycolors,mycolorsB)
    }
    ## potentially add grayscale from light "gray90" to dark "gray10"
    mygrays<-colorRampPalette(c("#E5E5E5","#1A1A1A"))

    ## textcolor within blocks
    if(is.null(carTextColors)){
        mytextcolors<-rep(c("white","black"),c(8,6))
        mytextcolorsB<-rep("white",33)
        mytextcolorsB[c(5,7,10:12,16:18,23,24,27:33)]<-"black"
        carTextColors<-c(mytextcolors,mytextcolorsB)
    }

    if(sortColsBySize==TRUE){
        mycarsall<-sort(table(compgenome$car),decreasing=TRUE)
    }else{
        mycarsall<-table(compgenome$car)
    }
    mycarsall<-as.numeric(names(mycarsall))

    if(uniqueCarColor==TRUE){
        if(length(mycarsall)>length(carColors)){
            ## add grayscale
            toadd<-length(mycarsall)-length(carColors)
            mycolors2<-c(carColors,mygrays(toadd))
            mytextcolors2<-c(carTextColors,rep(c("black","white"),
              c(ceiling(toadd/2),floor(toadd/2))))
        }else{
            mycolors2<-c(carColors)[1:length(mycarsall)]
            mytextcolors2<-c(carTextColors)[1:length(mycarsall)]
        }
    }

    ## colors for rearrangements
    ## -------------------------
    colorbreaks<-seq(-0.5,1.5,0.2) ## vals are only [0,1], so colors are less extreme
    ## reds for IV
    rPal<-colorRampPalette(c("#FFAAAA","#800000")) ## "#FF8080"
    ## purples for SM
    pPal<-colorRampPalette(c("#CCAAFF","#330080")) ## "#B380FF"
    ## blues for NM2
    bPal<-colorRampPalette(c("#AACCFF","#003380")) ## "#80B3FF"
    ## greens for NM1
    ##gPal<-colorRampPalette(c("#AAFFEE","#008066")) ## "#80FFE6"
    gPal<-colorRampPalette(c("#AAFFAA","#008000")) ## "#2AFF2A"


    ## PLOT
    ## -------------------------------------------

    if(makepdf==TRUE){
        pdf(filename,
            height=mylim$height,width=mylim$width,
            pointsize=12,paper="special",
            colormodel=colormodel)
    }else if(newdev==TRUE){
        dev.new(height=mylim$height,width=mylim$width,
                pointsize=12,noRStudioGD=TRUE)
    }

    ## ------
    op<-par(mar=mar,cex=1,cex.axis=cex.axis,cex.main=cex.main,
            lwd=1,xaxs="i",yaxs="i",font.main=font.main)
    ## ------


    plot(1,1,xlim=mylim$x,ylim=mylim$y,xlab="",ylab="",
         axes=FALSE,type="n",main=main)
    ##mtext(main,side=3,line=1.0,cex=1,font=3)

    myy2<-mylim$y[2]-pad
    scafftop<-numeric()

    mysubset<-match(ordfocal,names(BLOCKS))

    for(s in mysubset){

        ## make colors for cars
        if(uniqueCarColor==TRUE){ ## colors for cars stay the same across
            ## scaffolds - identify correct colors to use
            mycars<-BLOCKS[[s]]$blocks[,5]
            blockcolsCar<-mycolors2[match(mycars,mycarsall)]
            blockcolsText<-mytextcolors2[match(mycars,mycarsall)]
        }else if(uniqueCarColor==FALSE){ ## make colors based on number of
            ## markers for each scaffold separately
            mycars<-unique(BLOCKS[[s]]$blocks[,5])
            nm<-integer(length(mycars))
            for(i in 1:length(mycars)){
                tmp<-which(BLOCKS[[s]]$blocks[,5]==mycars[i])
                nm[i]<-sum(apply(BLOCKS[[s]]$blocks[tmp,c(1,2),drop=FALSE],1,
                                 function(x) diff(as.numeric(x))+1))
            }
            mycars<-mycars[order(nm,decreasing=TRUE)]
            if(length(mycars)>length(carColors)){
                ## add grayscale
                toadd<-length(mycars)-length(carColors)
                mycolors2<-c(carColors,mygrays(toadd))
                mytextcolors2<-c(carTextColors,
                                 rep(c("black","white"),
                                     c(ceiling(toadd/2),floor(toadd/2))))
            }else{
                mycolors2<-c(carColors)[1:length(mycars)]
                mytextcolors2<-c(carTextColors)[1:length(mycars)]
            }
            blockcolsCar<-character(nrow(BLOCKS[[s]]$blocks))
            blockcolsText<-character(nrow(BLOCKS[[s]]$blocks))
            for(q in 1:length(mycars)){
                tmp<-which(BLOCKS[[s]]$blocks[,5]==mycars[q])
                blockcolsCar[tmp]<-mycolors2[q]
                blockcolsText[tmp]<-mytextcolors2[q]
            }
        }

        scafftop<-c(scafftop,myy2-space$nmark)
        ## position for y-axis

        tmpx<-seq(blockwidth,blockwidth*nrow(BLOCKS[[s]]$blocks),blockwidth)
        ## positions on x-axis

        ## summary of alignment markers - tree
        for(k in 1:nrow(BLOCKS[[s]]$blocks)){

            tmpy2<-myy2
            ## number of markers in block
            tmpy1<-tmpy2-space$nmark
            text(x=tmpx[k]-blockwidth/2,y=mean(c(tmpy1,tmpy2)),
                 labels=1+diff(as.numeric(BLOCKS[[s]]$blocks[k,1:2])),adj=0.5,
                 col="black",cex=cex.text*0.7)
            ## cars
            tmpy2<-tmpy1
            tmpy1<-tmpy2-space$carid
            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                 xright=tmpx[k],ytop=tmpy2,col=blockcolsCar[k],
                 border=blockcolsCar[k])
            text(x=tmpx[k]-blockwidth/2,y=mean(c(tmpy1,tmpy2)),
                 labels=BLOCKS[[s]]$blocks[k,5],
                 adj=0.5,col=blockcolsText[k],cex=cex.text*1.0)
            ## summary for blocks in hierarchy
            tmpy2<-tmpy1
            for(n in 1:((ncol(BLOCKS[[s]]$blocks)-5)/9)){
                blockcol<-6+(n-1)*9
                ## start of columns for that hierarchy level
                if(is.element("Q",BLOCKS[[s]]$blocks[k,blockcol])){
                    ## nodeori
                    if(plotelem[1]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$ndori
                        if(BLOCKS[[s]]$blocks[k,blockcol+4]=="-1"){
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="black",border=blockcolsCar[k])
                            text(x=tmpx[k]-blockwidth/2,
                                 y=mean(c(tmpy1,tmpy2)),labels="-",
                                 adj=0.5,col="white",cex=cex.text*0.8)
                        }else if(BLOCKS[[s]]$blocks[k,blockcol+4]=="1"){
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="white",border=blockcolsCar[k])
                            text(x=tmpx[k]-blockwidth/2,
                                 y=mean(c(tmpy1,tmpy2)),labels="+",
                                 adj=0.5,col="black",cex=cex.text*0.8)
                        }else{
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="gray50",border=blockcolsCar[k])
                        }
                    }
                    ## blockid
                    if(plotelem[2]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$blkid
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col=NA,border=blockcolsCar[k])
                        text(x=tmpx[k]-blockwidth/2,y=mean(c(tmpy1,tmpy2)),
                             labels=BLOCKS[[s]]$blocks[k,blockcol+6],
                             adj=0.5,col="black",cex=cex.text*1.0)
                    }
                    ## blockori
                    if(plotelem[3]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$blkori
                        if(BLOCKS[[s]]$blocks[k,blockcol+7]=="-1"){
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="black",border=blockcolsCar[k])
                            text(x=tmpx[k]-blockwidth/2,
                                 y=mean(c(tmpy1,tmpy2)),labels="-",
                                 adj=0.5,col="white",cex=cex.text*0.8)
                        }else if(BLOCKS[[s]]$blocks[k,blockcol+7]=="1"){
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="white",border=blockcolsCar[k])
                            text(x=tmpx[k]-blockwidth/2,
                                 y=mean(c(tmpy1,tmpy2)),labels="+",
                                 adj=0.5,col="black",cex=cex.text*0.8)
                        }else{
                            rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                                 xright=tmpx[k],ytop=tmpy2,
                                 col="gray50",border=blockcolsCar[k])
                        }
                    }
                    ## block elements
                    if(plotelem[4]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$elem
                        if(BLOCKS[[s]]$blocks[k,blockcol+1]==BLOCKS[[s]]$blocks[k,blockcol+2]){
                            mylabels<-BLOCKS[[s]]$blocks[k,blockcol+1]
                        }else{
                            mylabels<-paste0(BLOCKS[[s]]$blocks[k,blockcol+1],
                                             ":",
                                             BLOCKS[[s]]$blocks[k,blockcol+2])
                        }
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col=NA,border=blockcolsCar[k])
                        text(x=tmpx[k]-blockwidth/2,y=mean(c(tmpy1,tmpy2)),
                             labels=mylabels,adj=0.5,col="black",
                             cex=cex.text*0.7)
                    }
                }else if(is.element("P",BLOCKS[[s]]$blocks[k,blockcol])){
                    ## nodeori
                    if(plotelem[1]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$ndori
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col="gray50",angle=45,density=12,lwd=1.6,border=NA)
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             border=blockcolsCar[k])
                    }
                    ## blockid
                    if(plotelem[2]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$blkid
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col=NA,border=blockcolsCar[k])
                        if(BLOCKS[[s]]$blocks[k,blockcol+6]!="0"){
                            text(x=tmpx[k]-blockwidth/2,y=mean(c(tmpy1,tmpy2)),
                                 labels=BLOCKS[[s]]$blocks[k,blockcol+6],
                                 adj=0.5,col="black",cex=cex.text*1.0)
                        }
                    }
                    ## blockori
                    if(plotelem[3]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$blkori
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col="gray50",angle=45,density=12,lwd=1.6,border=NA)
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             border=blockcolsCar[k])
                    }
                    ## block elements
                    if(plotelem[4]==1){
                        tmpy2<-tmpy1
                        tmpy1<-tmpy2-space$elem
                        rect(xleft=tmpx[k]-blockwidth,ybottom=tmpy1,
                             xright=tmpx[k],ytop=tmpy2,
                             col=NA,border=blockcolsCar[k])
                    }
                }
            }
        }


        ## rearrangements for blocks
        myy2<-myy2-((ncol(BLOCKS[[s]]$blocks)-5)/9)*
            (space$ndori+space$blkid+space$blkori+space$elem)-
                space$nmark-space$carid-space$stat
        myy1<-myy2-space$rearr

        if(plotelem[5]==1){
            ## NM1
            if(ncol(BLOCKS[[s]]$NM1)>0){
                for(l in 1:ncol(BLOCKS[[s]]$NM1)){
                    if(sum(BLOCKS[[s]]$NM1[,l]>remThld)>0){
                        for(k in 1:nrow(BLOCKS[[s]]$NM1)){
                            if(BLOCKS[[s]]$NM1[k,l]>remThld){
                                rect(xleft=tmpx[k]-blockwidth,ybottom=myy1,
                                     xright=tmpx[k],ytop=myy2,
                                     col=gPal(15)[cut(BLOCKS[[s]]$NM1[k,l],
                                         breaks=colorbreaks)],border=NA)
                            }
                        }
                        myy2<-myy1
                        myy1<-myy2-space$rearr
                    }
                }
            }

            ## NM2
            if(ncol(BLOCKS[[s]]$NM2)>0){
                for(l in 1:ncol(BLOCKS[[s]]$NM2)){
                    if(sum(BLOCKS[[s]]$NM2[,l]>remThld)>0){
                        for(k in 1:nrow(BLOCKS[[s]]$NM2)){
                            if(BLOCKS[[s]]$NM2[k,l]>remThld){
                                rect(xleft=tmpx[k]-blockwidth,ybottom=myy1,
                                     xright=tmpx[k],ytop=myy2,
                                     col=bPal(15)[cut(BLOCKS[[s]]$NM2[k,l],
                                         breaks=colorbreaks)],border=NA)
                            }
                        }
                        myy2<-myy1
                        myy1<-myy2-space$rearr
                    }
                }
            }

            ## SM
            if(ncol(BLOCKS[[s]]$SM)>0){
                for(l in 1:ncol(BLOCKS[[s]]$SM)){
                    if(sum(BLOCKS[[s]]$SM[,l]>remThld)>0){
                        for(k in 1:nrow(BLOCKS[[s]]$SM)){
                            if(BLOCKS[[s]]$SM[k,l]>remThld){
                                rect(xleft=tmpx[k]-blockwidth,ybottom=myy1,
                                     xright=tmpx[k],ytop=myy2,
                                     col=pPal(15)[cut(BLOCKS[[s]]$SM[k,l],
                                         breaks=colorbreaks)],border=NA)
                            }
                        }
                        myy2<-myy1
                        myy1<-myy2-space$rearr
                    }
                }
            }

            ## IV - blocks
            if(ncol(BLOCKS[[s]]$IV)>0){
                for(l in 1:ncol(BLOCKS[[s]]$IV)){
                    if(sum(BLOCKS[[s]]$IV[,l]>remThld)>0){
                        for(k in 1:nrow(BLOCKS[[s]]$IV)){
                            if(BLOCKS[[s]]$IV[k,l]>remThld){
                                rect(xleft=tmpx[k]-blockwidth,ybottom=myy1,
                                     xright=tmpx[k],ytop=myy2,
                                     col=rPal(15)[cut(BLOCKS[[s]]$IV[k,l],
                                         breaks=colorbreaks)],border=NA)
                            }
                        }
                        myy2<-myy1
                        myy1<-myy2-space$rearr
                    }
                }
            }

            ## IV - single markers
            if(ncol(BLOCKS[[s]]$IVsm)>0){
                for(l in 1:ncol(BLOCKS[[s]]$IVsm)){
                    for(k in 1:nrow(BLOCKS[[s]]$IVsm)){
                        if(BLOCKS[[s]]$IVsm[k,l]>0){
                            nmarker<-1+diff(as.numeric(BLOCKS[[s]]$blocks[k,1:2]))
                            xposns<-(tmpx[k]-blockwidth)+seq(0,blockwidth,length.out=nmarker+1)
                            xpos<-mean(xposns[(BLOCKS[[s]]$IVsm[k,l]):(BLOCKS[[s]]$IVsm[k,l]+1)])
                            segments(x0=xpos,x1=xpos,y0=myy1,y1=myy2,
                                     col=rPal(2)[2],lwd=4.0,lend=2)
                            ## rect(xleft=xposns[BLOCKS[[s]]$IVsm[k,l]],
                            ##      ybottom=myy1,
                            ##      xright=xposns[BLOCKS[[s]]$IVsm[k,l]+1],
                            ##      ytop=myy2,col=rPal(2)[2],border=NA)
                        }
                    }
                    myy2<-myy1
                    myy1<-myy2-space$rearr
                }
            }
        } ## close if(plotelem[5]==1)

        myy2<-myy1+space$rearr-space$scaff
    } ## close loop for(s in mysubset)

    axis(2,at=scafftop-space$carid/2,labels=yaxlab,las=2,padj=0.5,lwd=0,
         font=1,line= -1.2)

    ## ------
    par(op) ## set back to default
    ## ------

    if(makepdf==TRUE){
        dev.off()
        return(paste("Wrote",filename))
    }

    return(invisible(0))
}

## ------------------------------------------------------------------------
