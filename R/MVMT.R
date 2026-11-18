#' MVMT: Multivariate Multi-Trait Gene-Based Association Test
#'
#' Performs joint gene-based association testing across multiple correlated
#' traits using SNP-level summary statistics and per-gene LD matrices.
#' The method accounts for both linkage disequilibrium (LD) among variants
#' and phenotypic correlation among traits via separable Kronecker covariance
#' structure and simultaneous whitening.
#'
#' @param bim_file Path to PLINK *.bim file (columns: chr, snp, bp, pos, A1, A2).
#' @param map_file Gene annotation file without header
#'   (columns: chr, start, end, gene).
#' @param assoc_files Character vector of GWAS summary file paths,
#'   one per trait. Each file must contain columns:
#'   snp, A1, A2, AF1, beta, se, P, N.
#' @param trait_names Character vector of trait names, corresponding to
#'   \code{assoc_files}.
#' @param LD_path Directory containing per-gene GCTB LD files
#'   (prefix: \code{<gene_id>.ldm.full}).
#' @param result_path Output directory for results.
#' @param prop Proportion of total eigen-variance to retain in truncated
#'   spectral decomposition of the LD matrix (default 0.9).
#' @param n_cores Number of parallel cores for gene-wise computation
#'   (passed to \code{parallel::mclapply}).
#' @param af_threshold Allele frequency difference threshold for
#'   allele-flip QC (default 0.2).
#'
#' @return A \code{data.frame} with one row per gene and columns:
#'   \item{Gene}{Gene identifier}
#'   \item{Chr}{Chromosome}
#'   \item{Start, End}{Gene coordinates}
#'   \item{No.SNPs}{Number of SNPs in the gene}
#'   \item{No.Traits}{Number of traits analyzed}
#'   \item{SNP_start, SNP_end}{First and last SNP in the region}
#'   \item{No.Eigenvalues}{Degrees of freedom (retained LD PCs × traits)}
#'   \item{Chisq_MVMT}{Chi-square test statistic}
#'   \item{P_MVMT}{P-value from chi-squared distribution}
#'
#' @importFrom data.table fread setnames as.data.table setkey
#' @importFrom dplyr inner_join filter arrange select slice pull mutate
#' @importFrom parallel mclapply
#' @export
#'
#' @examples
#' \dontrun{
#' results <- MVMT(
#'   bim_file    = "data/g1000_eur_maf001hwe1e10.bim",
#'   map_file    = "data/gene_coords.txt",
#'   assoc_files = c("data/HDL.ma", "data/LDL.ma"),
#'   trait_names = c("HDL", "LDL"),
#'   LD_path     = "data/LD",
#'   result_path = "results",
#'   prop        = 0.9,
#'   n_cores     = 10
#' )
#' }
#---------------MVMT-------------------#
MVMT <- function(bim_file,
                 map_file,
                 assoc_files,
                 trait_names,
                 LD_path,
                 result_path,
                 prop,
                 n_cores,
                 af_threshold = 0.2) {

  # ---------- 读数据 ----------
  pheno_str <- paste(trait_names, collapse = ",")
  gene_annotate <- fread(map_file, header = FALSE,
                         col.names = c("chr", "start", "end", "gene"))

  bim <- fread(bim_file,
               col.names = c("chr", "snp", "bp", "pos", "A1", "A2")) %>%
    mutate(chr = as.character(chr))

  snp2gene <- map_snp_to_gene(bim, gene_annotate,
                               window_start = 50, window_end = 50)$map

  dt <- inner_join(bim, snp2gene %>% select(snp, chr, gene),
                   by = c("snp", "chr")) %>%
    na.omit()

  # ---------- GWAS ----------
  gwas_list <- setNames(
    lapply(assoc_files, function(f) {
      g <- fread(f)
      setnames(g, 1:8, c("snp","A1","A2","AF1","beta","se","P","N"))
      g
    }),
    trait_names
  )

  # ---------- common SNP ----------
  common_snps <- Reduce(intersect, lapply(gwas_list, `[[`, "snp"))
  dt <- dt %>% filter(snp %in% common_snps)

  gene_list <- unique(dt$gene)

  # ---------- cor_e ----------
  null_snps <- Reduce(intersect,
    lapply(gwas_list, function(g) {
      g$snp[g$snp %in% common_snps & g$P > 0.01]
    })
  )

  Z_null <- sapply(gwas_list, function(g) {
    idx <- match(null_snps, g$snp)
    g$beta[idx] / g$se[idx]
  })

  cor_e <- cor(Z_null, use = "pairwise.complete.obs")

  gwas_dt_list <- lapply(gwas_list, function(g) {
    g <- as.data.table(g)
    setkey(g, snp)
    g
  })
  dt <- as.data.table(dt)
  setkey(dt, gene)

  # ---------- 并行计算 ----------
  all_res_list <- mclapply(gene_list, function(gene_id) {

    tryCatch({

      chr_name <- dt[gene == gene_id, chr][1]

      inLD_prefix <- file.path(LD_path,
                               paste0(gene_id, ".ldm.full"))

      if (!file.exists(paste0(inLD_prefix, ".info"))) {
        message("Skip ", gene_id, " : LD file not found")
        return(NULL)
      }

      LD_in   <- ReadLDMbin(inLD_prefix)
      info_in <- fread(paste0(inLD_prefix, ".info"), na.strings = "NA")

      colnames(LD_in) <- rownames(LD_in) <- info_in$ID

      gene_snps  <- dt[gene == gene_id, snp]
      common_snps_gene <- intersect(gene_snps, colnames(LD_in))

      LD_order   <- colnames(LD_in)[match(common_snps_gene, colnames(LD_in))]
      info_check <- info_in[match(LD_order, info_in$ID), ]

      # 速度优化：用data.table键值查找替代全量match
      gwas_noqc_list <- lapply(gwas_dt_list, function(g)
        g[common_snps_gene]
      )

      # 内存优化：用完info_in立刻释放
      rm(info_in, gwas_noqc_list)

      gwas_flipped_list <- lapply(gwas_dt_list, function(g) {
        g_sub <- g[common_snps_gene]
        allele_flip_gwas_beta(g_sub, info_check, no_AF = FALSE, af_threshold = 0.2)
      })

      if (any(sapply(gwas_flipped_list, nrow) <= 10)) {
        message("Too few SNPs in ", gene_id)
        return(NULL)
      }

      snplist_final <- Reduce(intersect,
                              c(lapply(gwas_flipped_list, `[[`, "snp"),
                                list(colnames(LD_in))))

      snplist_final <- intersect(snplist_final, common_snps_gene)

      LD_idx <- colnames(LD_in)[match(snplist_final, colnames(LD_in))]
      R <- as.matrix(LD_in[LD_idx, LD_idx])
      rm(LD_in)
      gc(verbose = FALSE)
      gwas_common_list <- lapply(gwas_flipped_list, function(g)
        g[match(snplist_final, g$snp), ])

      n_snp <- nrow(R)
      eig_R <- eigen(R, symmetric = TRUE)
      rm(R)

      Z_matrix <- sapply(trait_names, function(trait) {
        gwas_common_list[[trait]]$beta /
          gwas_common_list[[trait]]$se
      })

      rownames(Z_matrix) <- LD_idx
      res_pre_prop <- MVMT_SVD(
        gwas_common_list, Z_matrix, eig_R, cor_e,
        prop, chr_name, gene_id, gene_annotate
      )

      return(res_pre_prop)

    }, error = function(e) {
      message("ERROR in gene: ", gene_id)
      message("Reason: ", e$message)
      return(NULL)
    })

  }, mc.cores = n_cores, mc.preschedule = TRUE)

  all_res <- do.call(rbind, all_res_list)
  if (is.null(all_res) || nrow(all_res) == 0) {
    warning("No valid gene results produced.")
    return(invisible(NULL))
   }
  out_file <- file.path(result_path,
                        paste0("MVMT_", pheno_str, ".txt"))
 if (!dir.exists(result_path)) dir.create(result_path, recursive = TRUE)
data.table::fwrite(all_res, file = out_file, sep = "\t", quote = FALSE)

return(all_res)
}