#' rearrvisr: A package to detect, classify, and visualize genome
#' rearrangements
#'
#' The rearrvisr package provides functions to identify and visualize inter- and
#' intrachromosomal translocations and inversions between a focal genome and an
#' ancestral genome reconstruction, or two extant genomes. Rearrangements,
#' breakpoints, and synteny blocks are identified along the focal genome and
#' output in a tabular format. Rearrangements and synteny blocks can be
#' visualized along the focal genome by two graphical functions.
#'
#' @section Functions to convert and verify input files:
#'
#'   \itemize{
#'
#'   \item \code{\link{convertPQtree}} converts linearly encoded \emph{PQ-trees}
#'   of an ancestral genome reconstruction, for example as output by the
#'   software ANGES (Jones \emph{et al.} 2012), into a two-dimensional
#'   \emph{PQ-structure} (i.e., a data frame representing the compared genome).
#'
#'   \item \code{\link{orderGenomeMap}} orders a one-dimensional genome map by
#'   genome segments (i.e., chromosomes or scaffolds) and by the position of
#'   markers within genome segments.
#'
#'   \item \code{\link{genome2PQtree}} converts an ordered, one-dimensional
#'   genome map into a two-dimensional \emph{PQ-structure} (i.e., a data frame
#'   representing the compared genome).
#'
#'   \item \code{\link{checkInfile}} checks input data to ensure that their file
#'   formats are correct.
#'
#'   }
#'
#' @section Functions to identify and summarize rearrangements:
#'
#'   \itemize{
#'
#'   \item \code{\link{computeRearrs}} detects and classifies rearrangements
#'   along a focal genome relative to an ancestral genome reconstruction or an
#'   extant genome. The focal genome has to be provided as an ordered,
#'   one-dimensional genome map, and the compared genome as a two-dimensional
#'   \emph{PQ-structure}.
#'
#'   \item \code{\link{filterRearrs}} filters detected rearrangements by their
#'   size.
#'
#'   \item \code{\link{getBreakpoints}} extracts breakpoint coordinates.
#'
#'   \item \code{\link{summarizeBlocks}} summarizes rearrangements and
#'   information on the alignment between the focal genome and the compared
#'   genome for each synteny block.
#'
#'   }
#'
#' @section Functions to visualize rearrangements:
#'
#'   \itemize{
#'
#'   \item \code{\link{genomeImagePlot}} generates a plot that shows different
#'   classes of rearrangements along a given set of focal genome segments.
#'
#'   \item \code{\link{genomeRearrPlot}} generates a plot that shows synteny
#'   blocks between a focal and a compared genome in columns, and information on
#'   their alignment and rearrangements in rows, for a given set of focal genome
#'   segments.
#'
#'   }
#'
#' @section Data:
#'
#'   The package includes example data (\code{MEL_markers}, \code{SIM_markers},
#'   \code{YAK_markers}, and \code{MSSYE_PQTREE_HEUR}) generated from 12
#'   publicly available \emph{Drosophila} genome assemblies downloaded from
#'   Ensemble Release 91
#'   \href{http://dec2017.archive.ensembl.org}{(http://dec2017.archive.ensembl.org)}
#'   and Ensemble Metazoa Release 37
#'   \href{http://oct2017-metazoa.ensembl.org}{(http://oct2017-metazoa.ensembl.org)}.
#'   Orthologs were identified with OMA standalone v2.2.0 (Altenhoff \emph{et
#'   al.} 2015). \code{MEL_markers}, \code{SIM_markers}, and \code{YAK_markers}
#'   are maps of extant genomes from \emph{D. melanogaster}, \emph{D. simulans},
#'   and \emph{D. yakuba}, respectively. \code{MSSYE_PQTREE_HEUR} is a genome
#'   reconstruction of the ancestor of the \emph{melanogaster} subgroup
#'   (Drosophila 12 Genomes Consortium 2007), computed with the software ANGES
#'   (Jones \emph{et al.} 2012). See the package vignette for details.
#'
#' @section References:
#'
#'   Altenhoff, A.M. \emph{et al.} (2015). The OMA orthology database in 2015:
#'   function predictions, better plant support, synteny view and other
#'   improvements. \emph{Nucleic Acids Research}, \strong{43}, D240--D249. doi:
#'   \href{https://doi.org/10.1093/nar/gku1158}{10.1093/nar/gku1158}.
#'
#'   Drosophila 12 Genomes Consortium (2007). Evolution of genes and genomes on
#'   the \emph{Drosophila} phylogeny. \emph{Nature}, \strong{450}, 203--218.
#'   doi: \href{https://doi.org/10.1038/nature06341}{10.1038/nature06341}.
#'
#'   Jones, B. R. \emph{et al.} (2012). ANGES: reconstructing ANcestral GEnomeS
#'   maps. \emph{Bioinformatics}, \strong{28}, 2388--2390. doi:
#'   \href{https://doi.org/10.1093/bioinformatics/bts457}{10.1093/bioinformatics/bts457}
#'
#' @examples
#' ## verify input format of focal genome:
#' checkInfile(MEL_markers, "focalgenome", checkorder = TRUE)
#' ## convert ancestral genome reconstruction to PQ-structure:
#' Comp_genome <- convertPQtree(MSSYE_PQTREE_HEUR)
#' ## alternatively, convert extant genome to PQ-stucture:
#' ## (note that minor scaffolds need to be excluded)
#' Comp_genome <- genome2PQtree(SIM_markers[is.element(SIM_markers$scaff,
#'                                                     c("2L", "2R", "3L", "3R", "X")), ])
#' ## verify input format of compared genome:
#' checkInfile(Comp_genome, "compgenome", checkorder = TRUE)
#'
#' ## identify and summarize rearrangements:
#' SYNT <- computeRearrs(MEL_markers, Comp_genome, doubled = TRUE)
#' BLOCKS <- summarizeBlocks(SYNT, MEL_markers, Comp_genome,
#'                           c("2L", "2R", "3L", "3R", "X"))
#'
#' ## visualize rearrangements:
#' genomeImagePlot(SYNT, MEL_markers, c("2L", "2R", "3L", "3R", "X"))
#' genomeRearrPlot(BLOCKS, Comp_genome, c("2L", "2R", "3L", "3R", "X"),
#'                 blockwidth = 1.15, y0pad = 3)
#'
#' @docType package
#' @name rearrvisr
NULL
