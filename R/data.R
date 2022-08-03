#' Genomic Relationship matrix of 502 lines
#'
#' The genomic relationship matrix is a symmetric, square matrix
#' giving the average sharing of genotypes across the genome of pairs 
#' of lines, and used to estimate heritability and do genomic prediction.
#' This GRM is based on lines used in the GenomesToFields initiative,
#' but line names have been anonymized.
#'
#' @format A square matrix with row and column names
"K"

#' Yield BLUPs for 502 lines across 19 trials
#'
#' @format A data frame with 3318 rows and 4 columns
#' \describe{
#'   \item{Line}{Line identifier, corresponds to K}
#'   \item{Population}{Population identifier}
#'   \item{Env}{Trial identifier}
#'   \item{Yield}{Yield BLUP per Line:Env}
#' }
"yield_data"