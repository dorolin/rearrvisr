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
#'
#' @examples
#' \dontrun{
#'
#' ## recreate the data set
#' TOY24_rawtree <- matrix(
#'   c(">TOY",
#'     "#CAR1",
#'     "_Q 1 2 3 4 5 -6 7 8 Q_",
#'     "#CAR2",
#'     "_Q 9 _Q 10 Q_ _Q 11 12 13 Q_ Q_",
#'     "#CAR3",
#'     "_Q 14 Q_",
#'     "#CAR4",
#'     "_Q _P 15 16 17 18 P_ _Q _Q -19 20 21 Q_ _Q -22 -23 24 Q_ Q_ Q_"),
#'   nrow = 9)
#' TOY24_compgenome <- convertPQtree(TOY24_rawtree)
#' }
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
#'
#' @examples
#' \dontrun{
#'
#' ## recreate the data set
#' TOY24_focalgenome <- data.frame(
#'   marker = as.integer(c(1,7,2,6:4,8:10,3,13:11,14,17:15,18,21,20,22:24,19)),
#'   scaff = as.character(rep(c(1,2,3), times = c(7,11,6))),
#'   start = as.integer(c(seq(10^6, by = 10^6, length.out = 7),
#'                        seq(10^6, by = 10^6, length.out = 11),
#'                        seq(10^6, by = 10^6, length.out = 6))),
#'   end = as.integer(c(seq(10^6+2, by = 10^6, length.out = 7),
#'                      seq(10^6+2, by = 10^6, length.out = 11),
#'                      seq(10^6+2, by = 10^6, length.out = 6))),
#'   strand = rep(rep(c("+", "-"), 4),
#'                times = c(4,2,4,3,8,1,1,1)),
#'   stringsAsFactors = FALSE)
#' }
"TOY24_focalgenome"




