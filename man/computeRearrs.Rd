% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeRearrs.R
\name{computeRearrs}
\alias{computeRearrs}
\title{Compute Rearrangements}
\usage{
computeRearrs(focalgenome, compgenome, doubled, remWgt = 0.05,
  splitnodes = TRUE, testlim = 100)
}
\arguments{
\item{focalgenome}{Data frame representing the focal genome, containing the
mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
\code{$end}, and \code{$strand}, and optional further columns. Markers need
to be ordered by their map position.}

\item{compgenome}{Data frame representing the compared genome (e.g., an
ancestral genome reconstruction, or an extant genome), with the first three
columns \code{$marker}, \code{$orientation}, and \code{$car}, followed by
columns that alternate between node type and node element. Markers need to
be ordered by their node elements.}

\item{doubled}{Logical. If \code{TRUE}, markers in the ancestral genome
reconstruction contain information about their orientation.}

\item{remWgt}{A numeric value between \code{0} (inclusive) and \code{0.5}
(exclusive). Sets the tagging weight for the component of a rearrangement
that is less parsimonious to have changed position relative to the
alternative component to \code{remWgt}, and that of the alternative
component to \code{1 - remWgt}.}

\item{splitnodes}{Logical. Split nodes into subnodes according to
rearrangements that occurred one step further up the hierarchy during the
rearrangement detection algorithm. \code{splitnodes = TRUE} prevents that
the same rearrangement receives tags across multiple levels of the
\emph{PQ-tree} hierarchy.}

\item{testlim}{A positive integer specifying the maximum number of tests 
performed to detect markers part of complex rearrangements. A lower value 
can improve speed, but might lead to less optimal results. Set to 
\code{Inf} for exhaustive testing (not recommended for highly rearranged
genomes).}
}
\value{
A list of matrices that store data on different classes of
  rearrangements and additional information on the structure of each
  \emph{PQ-tree} and its alignment to the focal genome. Markers are in rows,
  and the row names of each matrix correspond to the IDs in the
  \code{$marker} column of the \code{focalgenome} and \code{compgenome} data
  frames. The matrices contain all markers common to \code{focalgenome} and
  \code{compgenome}, and are ordered by their position in \code{focalgenome}.

  The list elements \code{$NM1}, \code{$NM2}, \code{$SM}, and \code{$IV}
  are numeric matrices that store identified rearrangements. \code{$NM1}
  stores \code{T}rans\code{L}ocations between CARs \code{B}etween focal
  \code{S}egments; \code{$NM2} stores \code{T}rans\code{L}ocations between
  CARs \code{W}ithin focal \code{S}egments; \code{$SM} stores
  \code{T}rans\code{L}ocations within CARs \code{W}ithin focal
  \code{S}egments; \code{$IV} stores \code{I}n\code{V}ersions within CARs
  within focal segments. See the package vignette for a detailed explanation
  of these classes of rearrangements.

  Each rearrangement is represented by a separate column. Except for
  \code{NM1}, which are identified across all focal segments, columns for
  individual focal segments are joined across rows to save space (i.e., for
  \code{NM2}, \code{SM}, and \code{IV}, which are identified within focal
  segments). To preserve the tabular format, these matrices are padded by
  zeros for focal segments with a non-maximal number of rearrangements, if
  necessary. If no rearrangements were detected for a certain class, the
  matrix has zero columns. Markers that are part of a rearrangement have a
  tag value of \code{>0} within their respective column. Tagged markers
  within a column are not necessarily consecutive, for example, when a
  rearrangement is split into several parts through an insertion of a
  different CAR, or when a rearrangement has an upstream and a downstream
  component (i.e., when alternative sets of markers may have caused an
  apparent nonsyntenic or syntenic move). Note that some columns in
  \code{$NM2} or \code{$SM} may be duplicated for a particular focal segment
  due to the functioning of the underlying algorithm; although corresponding
  to the same rearrangement, these duplicated columns are nevertheless
  included for completeness.

  For \code{NM1}, markers part of a class I nonsyntenic move have a value of
  \code{0.5} if non of the involved CAR fragments is a focal segment - CAR
  fragment \emph{best hit}. Otherwise, markers part of the CAR fragment that
  is assigned as focal segment - CAR fragment \emph{best hit} have a value of
  \code{0}, while markers part of all other non-\emph{best hit} CAR fragments
  have a value of \code{1}. For \code{NM2} and \code{SM}, markers part of a
  rearrangement with an upstream and a downstream component have a value of
  \code{1 - remWgt} (or \code{remWgt}) when they are part of the component
  that is more (or less) parsimonious to have changed position; if either
  component is equally parsimonious to have changed position, both have a
  value of \code{0.5}; all other markers part of a rearrangement have a value
  of \code{1}. For \code{IV}, markers part of an inversion have a value of
  \code{1}.

  The list elements \code{$NM1bS}, \code{$NM1bE}, \code{$NM2bS},
  \code{$NM2bE}, \code{$SMbS}, \code{$SMbE}, \code{$IVbS}, and \code{$IVbE}
  are numeric matrices that tag markers that denote the start (\code{$*bS})
  and end (\code{$*bE}) elements for the four classes of rearrangements
  (i.e., the markers adjacent to rearrangement breakpoints). Each
  rearrangement is represented by a separate column, but columns for
  individual focal segments are joined for all matrices across rows
  (including \code{$NM1bS} and \code{$NM1bE}) to save space. Tag values
  correspond to the ones in \code{$NM1}, \code{$NM2}, \code{$SM}, and
  \code{$IV}.

  The list elements \code{$nodeori}, \code{$blockori}, \code{$blockid},
  \code{$premask}, and \code{$subnode} are matrices that store information on
  the structure of each \emph{PQ-tree}, its alignment to the focal genome,
  and internal data. The first column of each matrix corresponds to the CAR
  level, and the following columns correspond to the hierarchical structure
  of each \emph{PQ-tree}, with information on the lowest level stored in the
  last column. \code{$nodeori} is a numeric matrix that stores the alignment
  direction of each \emph{Q-node} to the focal genome, with \code{1}
  indicating ascending (i.e., standard), and \code{-1} descending (i.e.,
  inverted) alignment. \emph{Q-nodes} that have no alignment direction (e.g.,
  single-marker nodes) have a value of \code{9}, and \emph{P-nodes} are
  \code{NA}. \code{$blockori} is a numeric matrix that stores the orientation
  of each synteny block, with \code{1} indicating ascending (i.e., standard),
  and \code{-1} descending (i.e., inverted) orientation. Blocks that have no
  orientation (e.g., blocks containing a single marker, or a single
  \emph{PQ-tree} branch) have a value of \code{9}, and blocks that are part
  of \emph{P-nodes} are \code{NA}. \code{$blockid} is a character matrix that
  stores the ID of each synteny block within its node. For \emph{Q-nodes},
  IDs are consecutive and start at \code{1}, separately for each node and
  each hierarchy level, and reflect the order of synteny blocks. Block IDs
  with \code{".1"} or \code{".2"} suffixes (in arbitrary order) indicate
  blocks that were subject to an additional subdivision step. For
  \emph{P-nodes}, IDs are \code{0} unless the node is part of a
  rearrangement, in which case IDs indicate different rearrangements, but not
  block order. \code{$premask} and \code{$subnode} are numeric matrices that
  store internal data used for the alignment and identification of
  rearrangements. Integers \code{>0} in \code{$subnode} indicate subdivisions
  of the corresponding \emph{PQ-tree} due to nonsyntenic or syntenic moves.
  All subdivisions have been searched separately for rearrangements one step
  further down the hierarchy. This is of main relevance when \code{splitnodes
  = TRUE}.

  The returned data can be visualized with the \code{\link{genomeImagePlot}}
  function, or summarized and visualized with the
  \code{\link{summarizeBlocks}} and \code{\link{genomeRearrPlot}} functions.
  The returned rearrangements can be filtered by size with the
  \code{\link{filterRearrs}} function. Breakpoint coordinates of
  rearrangements can be extracted with the \code{\link{getBreakpoints}}
  function.
}
\description{
Detect and classify rearrangements along a focal genome relative to an
ancestral genome reconstruction or an extant genome
}
\details{
\code{focalgenome} must contain the column \code{$marker}, a vector of
  either characters or integers with unique ortholog IDs that can be matched
  to the values in the \code{$marker} column of \code{compgenome}. Values can
  be \code{NA} for markers that have no ortholog. \code{$scaff} must be a
  character vector giving the name of the focal genome segment (i.e.,
  chromosome or scaffold) of origin of each marker. \code{$start} and
  \code{$end} must be numeric vectors giving the location of each marker on
  its focal genome segment. \code{$strand} must be a vector of \code{"+"} and
  \code{"-"} characters giving the reading direction of each marker.
  Additional columns are ignored and may store custom information, such as
  marker names. Markers need to be ordered by their map position within each
  focal genome segment, for example by running the
  \code{\link{orderGenomeMap}} function. See Examples below for the
  \code{focalgenome} format.

  \code{compgenome} must contain the column \code{$marker}, a vector of
  either characters or integers with unique ortholog IDs that can be matched
  to the values in the \code{$marker} column of \code{focalgenome}.
  \code{$orientation} must be a vector of \code{"+"} and \code{"-"}
  characters giving the reading direction of each marker in the compared
  genome. If \code{doubled = FALSE}, all values should be \code{"+"}.
  \code{$car} must be an integer vector giving the location of each marker on
  its compared genome segment (i.e., \emph{Contiguous Ancestral Region}, or
  CAR), analogous to contiguous sets of genetic markers on a chromosome,
  scaffold, or contig. Each CAR is represented by a \emph{PQ-tree} (Booth &
  Lueker 1976; Chauve & Tannier 2008). The \emph{PQ} structure of each CAR is
  defined by additional columns (at least two) that have to alternate between
  character vectors of node type (\code{"P"}, \code{"Q"}, or \code{NA}) in
  even columns, and integer vectors of node elements in odd columns (missing
  values are permitted past the fifth column). Every set of node type and
  node element column reflects the hierarchical structure of each
  \emph{PQ-tree}, with the rightmost columns representing the lowest level of
  the hierarchy. \emph{P-nodes} contain contiguous markers and/or nodes in
  arbitrary order, while \emph{Q-nodes} contain contiguous markers and/or
  nodes in fixed order (including their reversal). For additional details on
  \emph{PQ-trees} see Booth & Lueker 1976, Chauve & Tannier 2008, or the
  package vignette. See Examples below for the \code{compgenome} format.

  \code{doubled = TRUE} indicates that orientation information for the
  markers in the ancestral genome reconstruction is available. (This is the
  case for example when the genome was reconstructed with the software ANGES,
  Jones \emph{et al.} 2012, using the option \code{markers_doubled 1}.)
  Orientation information facilitates detecting and classifying
  rearrangements as inversions or syntenic moves, and can help determining
  whether \emph{PQ-tree} nodes are aligned to the focal genome in ascending
  (i.e., standard) or descending (i.e., inverted) direction.

  \code{remWgt} provides the tagging weight for rearrangements consisting of
  alternative sets of markers, either of which may have caused an apparent
  nonsyntenic or syntenic move (e.g., a set of markers may have moved
  upstream, or alternatively another set of markers may have moved
  downstream). The set of markers that is more parsimonious to have changed
  position relative to the other set receives tag values equal \code{1 -
  remWgt}, while the alternative set of markers receives tag values equal
  \code{remWgt}. Setting this argument to non-default may require adjusting
  the \code{remThld} argument in the \code{genomeImagePlot} and
  \code{renomeRearrPlot} functions accordingly.
}
\section{Algorithm}{


  A detailed description of the implemented algorithm can be found in the
  Supplementary information of the manuscript associated with the package
}

\section{References}{


  Booth, K.S. & Lueker, G.S. (1976). Testing for the consecutive ones
  property, interval graphs, and graph planarity using \emph{PQ}-Tree
  algorithms. \emph{Journal of Computer and System Sciences}, \strong{13},
  335--379. doi:
  \href{https://doi.org/10.1016/S0022-0000(76)80045-1}{10.1016/S0022-0000(76)80045-1}.

  Chauve, C. & Tannier, E. (2008). A methodological framework for the
  reconstruction of contiguous regions of ancestral genomes and its
  application to mammalian genomes. \emph{PLOS Computational Biology},
  \strong{4}, e1000234. doi:
  \href{https://doi.org/10.1371/journal.pcbi.1000234}{10.1371/journal.pcbi.1000234}.

  Jones, B. R. \emph{et al.} (2012). ANGES: reconstructing ANcestral GEnomeS
  maps. \emph{Bioinformatics}, \strong{28}, 2388--2390. doi:
  \href{https://doi.org/10.1093/bioinformatics/bts457}{10.1093/bioinformatics/bts457}
}

\examples{
computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)

\dontrun{

## focalgenome format:
TOY24_focalgenome

## compgenome format:
TOY24_compgenome
}

}
\seealso{
\code{\link{filterRearrs}}, \code{\link{genomeImagePlot}},
  \code{\link{getBreakpoints}}, \code{\link{summarizeBlocks}},
  \code{\link{genomeRearrPlot}}, \code{\link{summarizeRearrs}}; 
  \code{\link{orderGenomeMap}} to order the
  \code{focalgenome} data frame; \code{\link{convertPQtree}} or
  \code{\link{genome2PQtree}} to generate the \code{compgenome} data frame.
}
