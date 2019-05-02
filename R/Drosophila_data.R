#' Linearly encoded \emph{PQ-tree} of the \emph{Drosophila} ancestral genome
#' \emph{MSSYE}
#'
#' A data set representing the \emph{Drosophila} ancestral genome \emph{MSSYE}
#' with 8,973 markers on 20 \emph{Contiguous Ancestral Regions} (CARs),
#' reconstructed with the software ANGES (Jones \emph{et al.} 2012).
#' \emph{MSSYE} represents the ancestor of the \emph{melanogaster} subgroup
#' (Drosophila 12 Genomes Consortium 2007; i.e., separating the group containing
#' \emph{D. melanogaster}, \emph{D. simulans}, \emph{D. sechellia}, \emph{D.
#' yakuba}, and \emph{D. erecta} from the remainder of the \emph{Drosophila}
#' species). The genomes used to reconstruct the ancestral genome \emph{MSSYE}
#' were downloaded from Ensemble Release 91
#' \href{http://dec2017.archive.ensembl.org}{(http://dec2017.archive.ensembl.org)}
#' and Ensemble Metazoa Release 37
#' \href{http://oct2017-metazoa.ensembl.org}{(http://oct2017-metazoa.ensembl.org)}.
#' Orthologs were identified with OMA standalone v2.2.0 (Altenhoff \emph{et al.}
#' 2015). See the package vignette for details.
#'
#' @format A data frame with a single column encoding an ancestral genome
#'   reconstruction as strings. The first row gives the name of the ancestor
#'   (preceded by \code{>}), followed by two rows for each ancestral genome
#'   segment (CAR). The first row of each set gives the CAR ID (preceded by
#'   \code{#CAR}), and the second row provides the \emph{PQ-tree} structure of
#'   that CAR. The children of each node in the \emph{PQ-tree} are enclosed by
#'   \code{_P} and \code{P_} markups for \emph{P-nodes}, or by \code{_Q} and
#'   \code{Q_} markups for \emph{Q-nodes}. Markers (i.e., ortholog IDs) that
#'   belong to a particular node are located between its corresponding markups.
#'   Markers with reversed orientation are preceded by a \code{-} sign. The
#'   opening (\code{_P} or \code{_Q}) and closing (\code{P_} or \code{Q_})
#'   markups can be nested to allow the representation of the hierarchical
#'   structure of the \emph{PQ-tree}. For details on \emph{PQ-trees} see Booth &
#'   Lueker 1976, Chauve & Tannier 2008, or the package vignette.
#'
#' @section References:
#'
#'   Altenhoff, A.M. \emph{et al.} (2015). The OMA orthology database in 2015:
#'   function predictions, better plant support, synteny view and other
#'   improvements. \emph{Nucleic Acids Research}, \strong{43}, D240--D249. doi:
#'   \href{https://doi.org/10.1093/nar/gku1158}{10.1093/nar/gku1158}.
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
#'   Drosophila 12 Genomes Consortium (2007). Evolution of genes and genomes on
#'   the \emph{Drosophila} phylogeny. \emph{Nature}, \strong{450}, 203--218.
#'   doi: \href{https://doi.org/10.1038/nature06341}{10.1038/nature06341}.
#'
#'   Jones, B. R. \emph{et al.} (2012). ANGES: reconstructing ANcestral GEnomeS
#'   maps. \emph{Bioinformatics}, \strong{28}, 2388--2390. doi:
#'   \href{https://doi.org/10.1093/bioinformatics/bts457}{10.1093/bioinformatics/bts457}
#'
#' @seealso \code{\link{convertPQtree}}, \code{\link{computeRearrs}},
#'   \code{\link{MEL_markers}}, \code{\link{SIM_markers}},
#'   \code{\link{YAK_markers}}
"MSSYE_PQTREE_HEUR"

#' Focal genome map (i.e., \code{focalgenome}) of \emph{Drosophila melanogaster}
#'
#' A data set representing the \emph{D. melanogaster} genome with 13,818 markers
#' on 12 genome segments (i.e., chromosomes and scaffolds). The genome map was
#' generated based on publicly available genome assemblies downloaded from
#' Ensemble Release 91
#' \href{http://dec2017.archive.ensembl.org}{(http://dec2017.archive.ensembl.org)}
#' and Ensemble Metazoa Release 37
#' \href{http://oct2017-metazoa.ensembl.org}{(http://oct2017-metazoa.ensembl.org)}.
#' Orthologs were identified with OMA standalone v2.2.0 (Altenhoff \emph{et al.}
#' 2015). See the package vignette for details.
#'
#' @format A data frame with markers in rows and six columns that specify the
#'   structure of the \emph{D. melanogaster} focal genome map:
#' \tabular{ll}{
#'   \code{$marker} \tab unique ortholog ID of marker \cr
#'   \code{$scaff} \tab genome segment where marker is located \cr
#'   \code{$start} \tab start position of marker on genome segment \cr
#'   \code{$end} \tab end position of marker on genome segment \cr
#'   \code{$strand} \tab reading direction of marker \cr
#'   \code{$name} \tab marker name \cr
#' }
#' @section References:
#'
#'   Altenhoff, A.M. \emph{et al.} (2015). The OMA orthology database in 2015:
#'   function predictions, better plant support, synteny view and other
#'   improvements. \emph{Nucleic Acids Research}, \strong{43}, D240--D249. doi:
#'   \href{https://doi.org/10.1093/nar/gku1158}{10.1093/nar/gku1158}.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{SIM_markers}},
#'   \code{\link{YAK_markers}}, \code{\link{MSSYE_PQTREE_HEUR}}
"MEL_markers"

#' Focal genome map (i.e., \code{focalgenome}) of \emph{Drosophila simulans}
#'
#' A data set representing the \emph{D. simulans} genome with 15,355 markers on
#' 1,038 genome segments (i.e., chromosomes and scaffolds). The genome map was
#' generated based on publicly available genome assemblies downloaded from
#' Ensemble Release 91
#' \href{http://dec2017.archive.ensembl.org}{(http://dec2017.archive.ensembl.org)}
#' and Ensemble Metazoa Release 37
#' \href{http://oct2017-metazoa.ensembl.org}{(http://oct2017-metazoa.ensembl.org)}.
#' Orthologs were identified with OMA standalone v2.2.0 (Altenhoff \emph{et al.}
#' 2015). See the package vignette for details.
#'
#' @format A data frame with markers in rows and six columns that specify the
#'   structure of the \emph{D. melanogaster} focal genome map:
#' \tabular{ll}{
#'   \code{$marker} \tab unique ortholog ID of marker \cr
#'   \code{$scaff} \tab genome segment where marker is located \cr
#'   \code{$start} \tab start position of marker on genome segment \cr
#'   \code{$end} \tab end position of marker on genome segment \cr
#'   \code{$strand} \tab reading direction of marker \cr
#'   \code{$name} \tab marker name \cr
#' }
#' @section References:
#'
#'   Altenhoff, A.M. \emph{et al.} (2015). The OMA orthology database in 2015:
#'   function predictions, better plant support, synteny view and other
#'   improvements. \emph{Nucleic Acids Research}, \strong{43}, D240--D249. doi:
#'   \href{https://doi.org/10.1093/nar/gku1158}{10.1093/nar/gku1158}.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{MEL_markers}},
#'   \code{\link{YAK_markers}}, \code{\link{MSSYE_PQTREE_HEUR}}
"SIM_markers"

#' Focal genome map (i.e., \code{focalgenome}) of \emph{Drosophila yakuba}
#'
#' A data set representing the \emph{D. yakuba} genome with 16,039 markers on
#' 605 genome segments (i.e., chromosomes and scaffolds). The genome map was
#' generated based on publicly available genome assemblies downloaded from
#' Ensemble Release 91
#' \href{http://dec2017.archive.ensembl.org}{(http://dec2017.archive.ensembl.org)}
#' and Ensemble Metazoa Release 37
#' \href{http://oct2017-metazoa.ensembl.org}{(http://oct2017-metazoa.ensembl.org)}.
#' Orthologs were identified with OMA standalone v2.2.0 (Altenhoff \emph{et al.}
#' 2015). See the package vignette for details.
#'
#' @format A data frame with markers in rows and six columns that specify the
#'   structure of the \emph{D. melanogaster} focal genome map:
#' \tabular{ll}{
#'   \code{$marker} \tab unique ortholog ID of marker \cr
#'   \code{$scaff} \tab genome segment where marker is located \cr
#'   \code{$start} \tab start position of marker on genome segment \cr
#'   \code{$end} \tab end position of marker on genome segment \cr
#'   \code{$strand} \tab reading direction of marker \cr
#'   \code{$name} \tab marker name \cr
#' }
#' @section References:
#'
#'   Altenhoff, A.M. \emph{et al.} (2015). The OMA orthology database in 2015:
#'   function predictions, better plant support, synteny view and other
#'   improvements. \emph{Nucleic Acids Research}, \strong{43}, D240--D249. doi:
#'   \href{https://doi.org/10.1093/nar/gku1158}{10.1093/nar/gku1158}.
#'
#' @seealso \code{\link{computeRearrs}}, \code{\link{MEL_markers}},
#'   \code{\link{SIM_markers}}, \code{\link{MSSYE_PQTREE_HEUR}}
"YAK_markers"
