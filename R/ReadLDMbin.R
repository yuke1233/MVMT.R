#' Read LD matrix from GCTB binary format
#' @param prefix File prefix (without .info or .bin extension).
#' @return An \code{m × m} numeric matrix.
#' @keywords internal
#' @noRd
ReadLDMbin <- function(prefix) {
  InfoFile <- paste0(prefix, ".info")
  BinFile  <- paste0(prefix, ".bin")
  info <- utils::read.table(InfoFile, header = TRUE)
  m <- as.numeric(nrow(info))
  ldm <- matrix(
    readBin(BinFile, n = m * m, what = numeric(0), size = 4),
    nrow = m
  )
  return(ldm)
}
