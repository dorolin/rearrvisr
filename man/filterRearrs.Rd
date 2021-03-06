% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterRearrs.R
\name{filterRearrs}
\alias{filterRearrs}
\title{Filter Rearrangements}
\usage{
filterRearrs(SYNT, focalgenome, filterMin = c(NA, NA, NA, NA),
  filterMax = c(NA, NA, NA, NA))
}
\arguments{
\item{SYNT}{A list of matrices that store data on different classes of
rearrangements and additional information. \code{SYNT} must have been
generated with the \code{\link{computeRearrs}} function.}

\item{focalgenome}{Data frame representing the focal genome, containing the
mandatory columns \code{$marker}, \code{$scaff}, \code{$start},
\code{$end}, and \code{$strand}, and optional further columns. Markers need
to be ordered by their map position.}

\item{filterMin}{A numerical vector of the form \code{c(nm1, nm2, sm,
iv)} that specifies the minimum number of markers a rearrangement has to
comprise to be retained. \code{nm1} is the minimum number of markers in
\code{SYNT$NM1}, \code{nm2} is the minimum number of markers in
\code{SYNT$NM2}, \code{sm} is the minimum number of markers in
\code{SYNT$SM}, and \code{iv} is the minimum number of markers in
\code{SYNT$IV}.}

\item{filterMax}{A numerical vector of the form \code{c(nm1, nm2, sm,
iv)} that specifies the maximum number of markers a rearrangement is
allowed to comprise to be retained. \code{nm1} is the maximum number of
markers in \code{SYNT$NM1}, \code{nm2} is the maximum number of markers
in \code{SYNT$NM2}, \code{sm} is the maximum number of markers in
\code{SYNT$SM}, and \code{iv} is the maximum number of markers in
\code{SYNT$IV}.}
}
\value{
A filtered version of \code{SYNT}. An additional list element
  \code{$filter} is created that specifies the applied filter.

  Note that for rearrangements that have more than one component, only the
  component that falls in the specified filter range is removed. This may
  result in an overestimation of the number of breakpoints when a filtered
  version of \code{SYNT} is used as input for the
  \code{\link{summarizeRearrs}} function.
}
\description{
Remove rearrangements that comprise less than a minimum or more than a
maximum number of markers
}
\details{
Parameters \code{SYNT} and \code{focalgenome} need to be
  specified.

  \code{focalgenome} must contain the column \code{$marker}, a vector of
  either characters or integers with unique ortholog IDs that can be matched
  to the values in the rownames of \code{SYNT}. Values can be \code{NA} for
  markers that have no ortholog. \code{$scaff} must be a character vector
  giving the name of the focal genome segment (e.g., chromosome or scaffold)
  of origin of each marker. \code{$start} and \code{$end} must be numeric
  vectors giving the location of each marker on its focal genome segment.
  \code{$strand} must be a vector of \code{"+"} and \code{"-"} characters
  giving the reading direction of each marker. Additional columns are ignored
  and may store custom information, such as marker names. Markers need to be
  ordered by their map position within each focal genome segment, for example
  by running the \code{\link{orderGenomeMap}} function. \code{focalgenome}
  may contain additional rows that were absent when running the
  \code{\link{computeRearrs}} function. However, all markers present in
  \code{SYNT} need to be contained in \code{focalgenome}, with the subset of
  shared markers being in the same order.

  Rearrangements are stored in \code{SYNT} and include the following
  rearrangement classes: NM1 are class I nonsyntenic moves; NM2 are class II
  nonsyntenic moves; SM are syntenic moves; IV are inversions.
}
\examples{
SYNT <- computeRearrs(TOY24_focalgenome, TOY24_compgenome, doubled = TRUE)

## only retain inversions comprising at least two markers
SYNT_filt<-filterRearrs(SYNT, TOY24_focalgenome, filterMin = c(0, 0, 0, 2))

}
\seealso{
\code{\link{computeRearrs}}, \code{\link{genomeImagePlot}},
  \code{\link{summarizeRearrs}}.
}
