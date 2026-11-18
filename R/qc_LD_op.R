#' Filter and reorder LD matrix to QC-passed SNPs
#'
#' @param gwas Data frame with \code{snp} column.
#' @param LD Square numeric matrix with SNP IDs as row/column names.
#' @param qc_snp_list Data frame with \code{snp} column, or character vector of SNP IDs.
#' @return Filtered and reordered LD matrix.
#' @keywords internal
#' @noRd
qc_LD_op <- function(gwas, LD, qc_snp_list) {
  
  # 统一处理 vector 和 data.frame
  qc_snps <- if (is.data.frame(qc_snp_list)) qc_snp_list$snp else qc_snp_list
  
  snplist <- intersect(gwas$snp, qc_snps)
  if (length(snplist) == 0) {
    stop("No common SNPs found between gwas, qc_snp_list, and LD matrix.")
  }
  
  # 检查是否都在 LD 矩阵中
  missing_in_ld <- setdiff(snplist, colnames(LD))
  if (length(missing_in_ld) > 0) {
    stop("SNPs not found in LD matrix: ", paste(missing_in_ld, collapse = ", "))
  }
  
  LD_qced <- LD[snplist, snplist, drop = FALSE]
  return(LD_qced)
}