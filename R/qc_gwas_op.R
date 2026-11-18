#' Filter GWAS data to QC-passed SNPs
#'
#' @param gwas A data frame of GWAS summary statistics. Must contain a \code{snp} column.
#' @param qc_snp_list A data frame of SNPs passing QC. Must contain a \code{snp} column.
#' @return A subset of \code{gwas} containing only SNPs present in \code{qc_snp_list}.
#' @keywords internal
#' @noRd
qc_gwas_op <- function(gwas, qc_snp_list) {
  snplist <- intersect(gwas$snp, qc_snp_list$snp)
  gwas_qced <- gwas[match(snplist, gwas$snp), ]
  rownames(gwas_qced) <- NULL
  return(gwas_qced)
}