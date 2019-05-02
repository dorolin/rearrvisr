#' \code{compgenome} example for 24 markers, created for illustrative purposes
#'
#' A data set illustrating a compared (ancestral or extant) genome
#' representation with 24 markers on four \emph{Contiguous Ancestral Regions}
#' (CARs), where each CAR is represented by a \emph{PQ-tree}.
#'
#' @format A data frame with markers in rows and nine columns that specify the
#'   structure of a compared genome representation:
#' \tabular{ll}{
#'   \code{$marker} \tab unique ortholog ID of marker \cr
#'   \code{$orientation} \tab reading direction of marker \cr
#'   \code{$car} \tab CAR where marker is located \cr
#'   \code{$type1} \tab node type at highest level of \emph{PQ-tree} \cr
#'   \code{$elem1} \tab node element at highest level of \emph{PQ-tree} \cr
#'   \code{$type2} \tab node type at second level of \emph{PQ-tree} \cr
#'   \code{$elem2} \tab node element at second level of \emph{PQ-tree} \cr
#'   \code{$type3} \tab node type at third level of \emph{PQ-tree} \cr
#'   \code{$elem3} \tab node element at third level of \emph{PQ-tree} \cr
#' }
#' @seealso \code{\link{computeRearrs}}, \code{\link{TOY24_focalgenome}}
"TOY24_compgenome"

#' \code{focalgenome} example for 24 markers, created for illustrative purposes
#'
#' A data set illustrating a focal (extant) genome map with 24 markers on three
#' genome segments (e.g., scaffolds).
#'
#' @format A data frame with markers in rows and five columns that specify the
#'   structure of a focal genome map:
#' \tabular{ll}{
#'   \code{$marker} \tab unique ortholog ID of marker \cr
#'   \code{$scaff} \tab genome segment where marker is located \cr
#'   \code{$start} \tab start position of marker on genome segment \cr
#'   \code{$end} \tab end position of marker on genome segment \cr
#'   \code{$strand} \tab reading direction of marker \cr
#' }
#' @seealso \code{\link{computeRearrs}}, \code{\link{TOY24_compgenome}}
"TOY24_focalgenome"




