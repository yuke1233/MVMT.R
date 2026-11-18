#' MVMT Statistic via Simultaneous Whitening
#'
#' Computes the multi-trait gene-based test statistic by decorrelating
#' SNPs (LD whitening) and traits (trait whitening), then summing squared
#' elements. Equivalent to the quadratic form but avoids explicit
#' inversion of the full Kronecker product.
#'
#' @param gwas_common_list Named list of GWAS summary data frames per trait.
#' @param Z_matrix m × K matrix of SNP-level Z-statistics (SNPs × traits).
#' @param eig_R Output of \code{eigen(R, symmetric = TRUE)} for LD matrix R.
#' @param cor_e K × K between-trait correlation matrix.
#' @param prop Proportion of eigenvalue variance to retain (default 0.9).
#' @param chr_name Chromosome name.
#' @param ID Gene identifier.
#' @param gene_annotate Data frame with columns gene, chr, start, end.
#'
#' @return One-row data.frame: Gene, Chr, Start, End, No.SNPs, No.Traits,
#'   SNP_start, SNP_end, No.Eigenvalues, Chisq_mBAT_multi, P_MVMT.
#'
#' @details
#' Three steps: (1) truncate LD eigencomponents, (2) LD whitening
#' \code{Λ_r^{-1/2} U_r^T Z}, (3) trait whitening \code{Z_star V D^{-1/2}}.
#' Statistic is sum of squared elements, distributed as χ²_{rK}.
#'
#' @importFrom dplyr %>% filter pull
#' @importFrom stats pchisq
#' @export
MVMT_SVD <- function(gwas_common_list,
                           Z_matrix,
                           eig_R,
                           cor_e,
                           prop,
                           chr_name,
                           ID,
                           gene_annotate) {

  n_traits <- length(gwas_common_list)
  trait_names <- names(gwas_common_list)
  n_snps   <- nrow(gwas_common_list[[1]])

  # ---------------- LD truncation ----------------
  values  <- eig_R$values
  vectors <- eig_R$vectors

  posval <- values[values > 1e-10]

  if (prop == 1) {
    r <- length(posval)
  } else {
    r <- which(cumsum(posval / sum(posval)) > prop)[1]
  }

  lambda_r <- values[1:r]
  U_r      <- vectors[, 1:r]

  # ---------- LD whitening ----------
  Z_star <- crossprod(U_r, Z_matrix)
  Z_star <- sweep(Z_star, 1, sqrt(lambda_r), "/")

  # ---------- Trait whitening ----------
  eig_T <- eigen(cor_e, symmetric = TRUE)   # cor_e是n_traits×n_traits的小矩阵，不是瓶颈
  V <- eig_T$vectors
  D <- eig_T$values
  if (any(D <= 1e-10)) warning("Trait correlation matrix near-singular")
  
  Z_final <- Z_star %*% V
  Z_final <- sweep(Z_final, 2, sqrt(D), "/")

  # ---------- Test statistic ----------
  chisq_opt <- sum(Z_final^2)
  df_opt    <- r * n_traits

  P_mvmt <- pchisq(chisq_opt, df=df_opt, lower.tail=FALSE)

  res_mvmt <- data.frame(
    Gene = ID,
    Chr = chr_name,
    Start = as.numeric(gene_annotate %>% dplyr::filter(gene==ID, chr==chr_name) %>% dplyr::pull(start)),
    End   = as.numeric(gene_annotate %>% dplyr::filter(gene==ID, chr==chr_name) %>% dplyr::pull(end)),
    No.SNPs = n_snps,
    No.Traits = n_traits,
    SNP_start = gwas_common_list[[1]]$snp[1],
    SNP_end   = gwas_common_list[[1]]$snp[n_snps],
    No.Eigenvalues = df_opt,
    Chisq_mBAT_multi = chisq_opt,
    P_MVMT = P_mvmt
  )

  return(res_mvmt)
}