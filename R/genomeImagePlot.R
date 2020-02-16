## ------------------------------------------------------------------------
## make image plot showing rearrangements along focal genome
## ------------------------------------------------------------------------

#' Genome Image Plot
#'
#' Generate a plot that shows different classes of rearrangements along a given
#' set of focal genome segments
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
#'   that will be plotted. Have to match (a subset of) IDs in
#'   \code{focalgenome$scaff}.
#' @param remstr String that should be removed from the IDs in \code{ordfocal}
#'   to simplify the y-axis labels. Only relevant when \code{yaxlab = NULL}.
#' @param main Title of the plot.
#' @param remThld A numeric value between \code{0} (inclusive) and \code{0.5}
#'   (exclusive). Controls whether components of rearrangements that are less
#'   parsimonious to have changed position relative to the alternative
#'   components will be plotted. To plot all components, \code{remThld} needs to
#'   be smaller than \code{remWgt} used in the \code{\link{computeRearrs}}
#'   function.
#' @param vspace A numeric value between \code{0} and \code{5}. Controls the
#'   amount of vertical space between focal genome segments.
#' @param hspace A numeric value between \code{0} and \code{1} that controls the
#'   proportion of horizontal space that is added to the last marker of a focal
#'   genome segment, calculated from the size of the largest focal genome
#'   segment.
#' @param mar A numerical vector of the form \code{c(bottom, left, top, right)}
#'   that specifies the margins on the four sides of the plot. The default
#'   \code{mar = NULL} sets the margins automatically.
#' @param xaxlab Data frame with two columns that gives the annotations of the
#'   x-axis. The column \code{$pos} specifies the position of tick-marks, and
#'   the column \code{$lab} specifies the labels of tick-marks. The default
#'   \code{xaxlab = NULL} makes annotations every 5 Mb.
#' @param xlab Title for the x-axis.
#' @param yaxlab Annotations of the y-axis. Must be the same length as
#'   \code{ordfocal}. The default \code{yaxlab = NULL} uses as annotations the
#'   names in \code{ordfocal}, simplified through removal of the string in
#'   \code{remstr}.
#' @param cex Numerical value that specifies the amount of text magnification.
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
#'   Parameters \code{SYNT}, \code{focalgenome}, and \code{ordfocal} need to be
#'   specified, all other parameters have default or automatic settings.
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
#'   shared markers being in the same order. Having additional markers in
#'   \code{focalgenome} can be useful for example to include additional
#'   non-orthologous markers in the plot.
#'
#'   When \code{makepdf = TRUE} or \code{newdev = TRUE}, the width and height of
#'   the graphic will be set automatically. The dimensions are determined in
#'   inches, thus \code{makepdf = FALSE} and \code{newdev = TRUE} will produce
#'   an error or not work correctly when the default units of the default
#'   graphics device are not inches (such as \code{\link[grDevices]{bmp}},
#'   \code{\link[grDevices]{jpeg}}, \code{\link[grDevices]{png}}, or
#'   \code{\link[grDevices]{tiff}}). This can be avoided by setting the default
#'   graphics device to one that has inches as default units. Setting both
#'   \code{makepdf = FALSE} and \code{newdev = FALSE} will allow to specify
#'   alternative, user-defined dimensions of the graphic. See examples below.
#'
#' @return A plot to the default graphics device (but not "RStudioGD") or a PDF
#'   file.
#'
#'   The plot visualizes the data contained in \code{SYNT} for each focal genome
#'   segment in \code{ordfocal}, arranged along the y-axis. The x-axis gives the
#'   map position of markers on their focal genome segment. Thin vertical lines
#'   (ticks) indicate the positions of markers, rearrangement breakpoints, and
#'   different classes of rearrangements within six rows per focal genome
#'   segment.
#'
#'   Starting at the bottom, the first row gives the map positions of all
#'   markers in \code{focalgenome} and of the subset of markers present in
#'   \code{SYNT} by gray and black ticks, respectively. The second row gives the
#'   positions of rearrangement breakpoints by red ticks. Ticks are drawn at the
#'   midpoints between the positions of the two markers present in \code{SYNT}
#'   that are adjacent to the breakpoints. Only markers with black ticks can
#'   receive colored ticks in the following rows. Maroon indicates markers that
#'   are part of inversions (IV); purple indicates markers that are part of
#'   syntenic moves (SM); blue indicates markers that are part of class II
#'   nonsyntenic moves (NM2); green indicates markers that are part of class I
#'   nonsyntenic moves (NM1).
#'
#'   Unless the argument \code{remThld} is set to a value smaller than that of
#'   \code{remWgt} in the \code{\link{computeRearrs}} function, only markers
#'   that are more parsimonious to have changed position relative to alternative
#'   markers are highlighted. Lighter tints denote markers that are part of
#'   large rearrangements, while darker shades denote markers that are part of
#'   small rearrangements. To distinguish between individual large
#'   rearrangements versus multiple short adjacent rearrangements, it may be
#'   helpful to take the position of breakpoints into account.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{filterRearrs}},
#'   \code{\link{getBreakpoints}}, \code{\link{genomeRearrPlot}}. For more
#'   information about arguments that are passed to other functions, see
#'   \code{\link[grDevices]{dev.new}}, \code{\link[grDevices]{pdf}},
#'   \code{\link[graphics]{plot}}, \code{\link[graphics]{par}}.
#'
#' @examples
#' SYNT<-computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)
#'
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2", "3"), main = "TOY24")
#'
#' ## change x-axis labels
#' op<-options(scipen=999)
#' myxaxlab<-data.frame(pos=seq(2*10^6,10^7,2*10^6),
#'                      lab=as.character(seq(2*10^3,10^4,2*10^3)))
#' options(op)
#' myxlab<-"Position [kb]"
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2", "3"), main = "TOY24",
#'                 xaxlab = myxaxlab, xlab = myxlab)
#'
#' \dontrun{
#'
#' ## make PDF (automatically determine the width and height of the graphic)
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2", "3"), main = "TOY24",
#'                 makepdf = TRUE, newdev = FALSE, filename = "genome.pdf")
#'
#' ## make PDF (default dimensions, i.e., square format)
#' pdf(file = "genome.pdf")
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2", "3"), main = "TOY24",
#'                 makepdf = FALSE, newdev = FALSE)
#' dev.off()
#'
#' ## plot in R Studio window
#' op <- options(device = "RStudioGD")
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2", "3"), main = "TOY24",
#'                 newdev = FALSE)
#' options(op)
#'
#' ## make EPS, and set user-specified dimensions
#' setEPS()
#' postscript("genome.eps", width=4.5,height=3.42,pointsize=9)
#' genomeImagePlot(SYNT, TOY24_focalgenome, c("1", "2"), main = "TOY24",
#'                 vspace = 1, makepdf = FALSE, newdev = FALSE)
#' dev.off()
#' }
#'
#' @export
#' @importFrom grDevices dev.new dev.off pdf colorRampPalette
#' @importFrom graphics par plot axis segments mtext rect


genomeImagePlot<-function(SYNT,focalgenome,ordfocal,
                          remstr = "",
                          main = "",
                          remThld = 0.05,
                          vspace = 0,
                          hspace = 0.01,
                          mar = NULL,
                          xaxlab = NULL,
                          xlab = "Position [Mb]",
                          yaxlab = NULL,
                          cex = 1, font.main = 3,
                          makepdf = FALSE,
                          newdev = TRUE,
                          filename = "genome.pdf",
                          colormodel = "srgb"){


    ## checks
    ## -------------------------------------------

    ## ERRORS
    ## check SYNT
    checkInfile(SYNT, myclass="SYNT")
    ## check focalgenome
    checkInfile(focalgenome, myclass="focalgenome", checkorder = TRUE)
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
    if(length(remstr) != 1){
        stop("remstr must be unique string of length one")
    }
    if(class(remstr) != "character"){
        stop("class for remstr needs to be 'character'")
    }
    if (length(remThld) != 1) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if (remThld < 0 | remThld >= 0.5) {
      stop("remThld needs to be a single numeric value within [0,0.5)")
    }
    if (length(xlab) != 1) {
      stop("label for x-axis must be string of length one")
    }
    if (length(vspace) != 1) {
      stop("vspace needs to be numeric value within [0,5]")
    }
    if (vspace < 0 | vspace > 5) {
      stop("vspace needs to be numeric value within [0,5]")
    }
    if (length(hspace) != 1) {
      stop("hspace needs to be a single numeric value within [0,1]")
    }
    if (hspace < 0 | hspace > 1) {
      stop("hspace needs to be a single numeric value within [0,1]")
    }
    if (length(makepdf) != 1 | !is.logical(makepdf)) {
      stop("makepdf must be 'TRUE' or 'FALSE'")
    }
    if (length(newdev) != 1 | !is.logical(newdev)) {
      stop("newdev must be 'TRUE' or 'FALSE'")
    }
    if (!is.element(colormodel, c("srgb", "cmyk", "gray", "grey"))) {
      stop("colormodel needs to be 'srgb', 'cmyk', 'gray', or 'grey'")
    }
    ## cex, font.main, filename should be checked by
    ##  'par()' or 'pdf()'/'dev.new()'

    ## check automatic settings
    ## simplified y-axis labels
    if(!is.null(yaxlab)){
        if(length(yaxlab) != length(ordfocal) |
           !is.vector(yaxlab)){
            stop("yaxlab must have the same length as ordfocal")
        }
    }
    ## x-axis
    if(!is.null(xaxlab)){
        if(!is.data.frame(xaxlab)){
            stop("xaxlab needs to be a data frame")
        }
        if(any(!is.element(c("pos","lab"),colnames(xaxlab)))){
            stop("column names in xaxlab need to include 'pos' and 'lab'")
        }
    }
    ## main label: can't think of any serious errors
    ## margins: should be checked by 'par()'


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

    ## colors for color gradient
    colorbreaks<-seq(0,1,0.01)
    ## reds for IV
    rPal<-colorRampPalette(c("#FFAAAA","#800000")) ## "#FF8080"
    ## purples for SM
    pPal<-colorRampPalette(c("#CCAAFF","#330080")) ## "#B380FF"
    ## blues for NM2
    bPal<-colorRampPalette(c("#AACCFF","#003380")) ## "#80B3FF"
    ## greens for NM1
    ##gPal<-colorRampPalette(c("#AAFFEE","#008066")) ## "#80FFE6"
    gPal<-colorRampPalette(c("#AAFFAA","#008000")) ## "#2AFF2A"


    ssize<-max(focalgenome$end[is.element(focalgenome$scaff,ordfocal)])
    xend<-ssize*hspace
    ssize<-ssize+xend
    ## add a few percent of padding to the right side of the x-axis
    xstart<-min(0,min(focalgenome$start[is.element(focalgenome$scaff,ordfocal)]) - xend)
    ## potentially add padding to the left side too

    ## make automatic settings
    ## simplified y-axis labels
    if(is.null(yaxlab)){
        yaxlab<-gsub(remstr,"",ordfocal)
    }
    ## x-axis
    if(is.null(xaxlab)){
        xaxlab<-data.frame(pos=seq(0,ssize,5000000),
                           lab=seq(0,ssize,5000000)/1000000)
    }
    ## margins
    if(is.null(mar)){
        mar<-c(3,0.6+max(nchar(yaxlab))*0.41,2,1)
    }


    ## PLOT
    ## -------------------------------------------

    if(makepdf==TRUE){
        pdf(filename,
            height=(length(ordfocal)+
                        (vspace/5)*(length(ordfocal)-1))*0.7+
                            sum(mar[c(1,3)])*0.2,
            width=7+sum(mar[c(2,4)])*0.2,
            pointsize=12,paper="special",
            colormodel=colormodel)
    }else if(newdev==TRUE){
        dev.new(height=(length(ordfocal)+
                            (vspace/5)*(length(ordfocal)-1))*0.7+
                                sum(mar[c(1,3)])*0.2,
                width=7+sum(mar[c(2,4)])*0.2,
                pointsize=12,noRStudioGD=TRUE)
    }
    ## mar*0.2 converts units in lines to inches

    ## ------
    op<-par(mar=mar,cex=cex,cex.axis=0.8,cex.main=1.2,
            lwd=1,xaxs="i",yaxs="i",font.main=font.main)
    ## ------


    plot(1,1,xlim=c(xstart,ssize),ylim=c(0,length(ordfocal)*(vspace+5)),
         xlab="",ylab="",axes=FALSE,type="n",main=main)
    axis(1,at=xaxlab$pos,labels=xaxlab$lab,
         line=0.0, lwd=0,lwd.ticks=1.0,padj=-1.6,tcl= -0.3)
    segments(x0=xstart,x1=ssize,y0=0,y1=0,lwd=1)
    mtext(xlab,side=1,line=1.6,cex=cex)
    ## cex within mtext is absolute value, not scaled by par(cex)
    ##mtext(main,side=3,line=0.5,cex=cex,font=3)


    ## for each scaffold
    for(s in 1:length(ordfocal)){

        mpos<-which(markers$scaff==ordfocal[s])
        ampos<-which(focalgenome$scaff==ordfocal[s])
        mysize<-max(focalgenome$end[focalgenome$scaff==ordfocal[s]])
        mysize<-mysize+xend
        ## myy<-seq(0,(length(ordfocal)-1)*(5+vspace),
        ##          5+vspace)[s]+vspace ## start plot at bottom
        myy<-rev(seq(0,(length(ordfocal)-1)*(5+vspace),
                     5+vspace))[s]+vspace ## start plot on top

        ## positions of all genes in focal genome
        mypos<-(focalgenome[ampos,3]+focalgenome[ampos,4])/2
        segments(x0=mypos,x1=mypos,y0=myy,y1=myy+0.5,col="gray",lwd=0.3)
        ## positions of genes included in SYNT
        mypos<-(markers[mpos,3]+markers[mpos,4])/2
        segments(x0=mypos,x1=mypos,y0=myy+0.25,y1=myy+0.5,col="black",lwd=0.3)
        ## breakpoints (combined for all rearrangement classes)
        tmpS<-cbind(SYNT$IVbS[mpos,],SYNT$SMbS[mpos,],
                    SYNT$NM2bS[mpos,],SYNT$NM1bS[mpos,])
        tmpE<-cbind(SYNT$IVbE[mpos,],SYNT$SMbE[mpos,],
                    SYNT$NM2bE[mpos,],SYNT$NM1bE[mpos,])
        if(sum(tmpS)>remThld | sum(tmpE)>remThld){
            brpts<-getBreakpnts2BP(tmpS,tmpE,markers$start[mpos],
                                   markers$end[mpos],remThld)
            if(nrow(brpts)>0){
                segments(x0=brpts$bptmid,x1=brpts$bptmid,
                         y0=myy+0.5,y1=myy+1,col="red",lwd=0.7)
            }
        }
        ## IV
        tmp<-summarizeTags(SYNT$IV[mpos,],mpos,
                           rownames(SYNT$IV)[mpos],ordfocal,s,remThld)
        if(sum(tmp)>0){
            tags<-which(tmp>0)
            mypos<-(markers[mpos[tags],3]+markers[mpos[tags],4])/2
            mycolors<-rPal(100)[as.numeric(cut(tmp[tags],breaks=colorbreaks))]
            segments(x0=mypos,x1=mypos,y0=myy+1,y1=myy+2,col=mycolors,lwd=0.7)
        }
        ## SM
        tmp<-summarizeTags(SYNT$SM[mpos,],mpos,
                           rownames(SYNT$SM)[mpos],ordfocal,s,remThld)
        if(sum(tmp)>0){
            tags<-which(tmp>0)
            mypos<-(markers[mpos[tags],3]+markers[mpos[tags],4])/2
            mycolors<-pPal(100)[as.numeric(cut(tmp[tags],breaks=colorbreaks))]
            segments(x0=mypos,x1=mypos,y0=myy+2,y1=myy+3,col=mycolors,lwd=0.7)
        }
        ## NM2
        tmp<-summarizeTags(SYNT$NM2[mpos,],mpos,
                           rownames(SYNT$NM2)[mpos],ordfocal,s,remThld)
        if(sum(tmp)>0){
            tags<-which(tmp>0)
            mypos<-(markers[mpos[tags],3]+markers[mpos[tags],4])/2
            mycolors<-bPal(100)[as.numeric(cut(tmp[tags],breaks=colorbreaks))]
            segments(x0=mypos,x1=mypos,y0=myy+3,y1=myy+4,col=mycolors,lwd=0.7)
        }
        ## NM1
        tmp<-summarizeTags(SYNT$NM1[mpos,],mpos,
                           rownames(SYNT$NM1)[mpos],ordfocal,s,remThld)
        if(sum(tmp)>0){
            tags<-which(tmp>0)
            mypos<-(markers[mpos[tags],3]+markers[mpos[tags],4])/2
            mycolors<-gPal(100)[as.numeric(cut(tmp[tags],breaks=colorbreaks))]
            segments(x0=mypos,x1=mypos,y0=myy+4,y1=myy+5,col=mycolors,lwd=0.7)
        }

        axis(2,at=myy+2.5,labels=yaxlab[s],las=2,padj=0.5,lwd=0,
             font=1,line= -0.6)
        rect(xleft=xstart,xright=mysize,ybottom=myy,ytop=myy+5,border="gray70",
             lwd=0.9)

    }

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
