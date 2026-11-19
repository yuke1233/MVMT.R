# MVMT power simulation
#
# Parameter grid (n_traits; shared_prop; rg; re):
# power: 2 5 7; 0.2; 0.2; 0.3 0.5 0.7
#
# Usage: Rscript MVMT_power.R <n_traits> <shared_prop> <rg> <re>
# Examples:
#   Rscript MVMT_power.R 2 0.2 0.2 0.3
#   Rscript MVMT_power.R 5 0.2 0.2 0.5
#   Rscript MVMT_power.R 7 0.2 0.2 0.7
# Submit all nine combinations: bash submit_power.sh
#
# n_traits: number of traits; shared_prop: proportion of causal SNPs shared
# by all traits; rg/re: genetic/residual correlation parameters.
# pi: causal SNP proportion per trait; h: heritability parameter for each trait.
# n_sim: replicates per gene; alpha: nominal P-value threshold.
# explained_threshold: cumulative positive LD eigenvalue proportion retained.
# n_cores: parallel::mclapply workers (Linux/macOS).
#
# Equal h = 0.0005 follows the earlier power script and supports 2, 5 or 7 traits.
# The simulation and statistical functions follow the newer power script.
# MTaSPUsSet uses 1,000 permutations. Summary means and standard errors are
# calculated across gene-specific rejection rates. Edit the paths before use.

library(data.table)
library(mvtnorm)
library(MASS)
library(genio)
library(aSPU)
library(parallel)
library(metaCCA)
library(MSKAT)

# ====================
# Parameter settings
# ====================
set.seed(123)

plink_prefix <- "data/plink_chr22_100gene"
gene_bed     <- "data/chr22_100gene.bed"
out_dir      <- "results/power"

pi <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 4) {
  n_traits    <- as.integer(args[1])
  shared_prop <- as.numeric(args[2])
  rg          <- as.numeric(args[3])
  re          <- as.numeric(args[4])
} else {
  n_traits    <- 2
  shared_prop <- 0.2
  rg          <- 0.2
  re          <- 0.3
}

h <- rep(0.0005, n_traits)
n_sim <- 100
alpha <- 0.05
explained_threshold <- 0.9
n_cores <- 30

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ====================
# Read genotypes and map SNPs to genes
# ====================

cat("Loading genotype data...\n")
geno_all <- read_plink(plink_prefix)
X_all    <- t(geno_all$X)
bim_all  <- geno_all$bim

X_all <- apply(X_all, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})
cat(sprintf("Genotype: %d samples x %d SNPs\n", nrow(X_all), ncol(X_all)))

if (!"pos" %in% colnames(bim_all)) bim_all$pos <- bim_all$bp
if (!"chr" %in% colnames(bim_all)) bim_all$chr <- bim_all$chrom

# ====================
# Read gene coordinates and map SNPs to genes
# ====================
gene_bed <- fread(gene_bed, header = FALSE,
                  col.names = c("chr", "start", "end", "gene"))

gene_snp_idx <- lapply(seq_len(nrow(gene_bed)), function(i) {
  g <- gene_bed[i,]
  which(bim_all$chr %in% c("22","chr22") &
        bim_all$pos >= g$start &
        bim_all$pos <= g$end)
})
names(gene_snp_idx) <- gene_bed$gene

keep <- sapply(gene_snp_idx, length) >= 20
gene_snp_idx <- gene_snp_idx[keep]
gene_list    <- names(gene_snp_idx)
cat(sprintf("Genes to simulate: %d\n\n", length(gene_list)))

# ====================
# Helper functions
# ====================

make_positive_definite <- function(mat, epsilon = 1e-6) {
  if (!isSymmetric(mat)) mat <- (mat + t(mat)) / 2
  mat <- mat + diag(epsilon, nrow(mat))
  eigen_vals <- eigen(mat, symmetric = TRUE)$values
  if (min(eigen_vals) < 0)
    mat <- as.matrix(Matrix::nearPD(mat, corr = FALSE)$mat)
  return(mat)
}

select_representative_snps <- function(R_mat, var_threshold = 0.95) {
  M <- nrow(R_mat)

  selected   <- integer(0)      # Selected SNP indices
  remaining  <- 1:M             # Remaining SNP indices
  resid_var  <- rep(1.0, M)     # Residual variance for each SNP, initially one
  partial_R  <- R_mat           # Current partial-correlation matrix

  repeat {
    # LD score for each candidate SNP, weighted by residual variances
    ld_scores <- sapply(remaining, function(g) {
      sum(partial_R[g, remaining]^2 * resid_var[remaining])
    })

    # Select the SNP with the largest LD score
    best_local  <- which.max(ld_scores)
    best_global <- remaining[best_local]
    selected    <- c(selected, best_global)
    remaining   <- setdiff(remaining, best_global)

    if (length(remaining) == 0) break

    # Remove the selected SNP effect from the remaining SNPs
    # Schur-complement update
    r_sg <- partial_R[best_global, remaining]   # Correlations with the remaining SNPs
    r_ss <- partial_R[best_global, best_global] # Variance of the selected SNP

    # Update residual variances
    resid_var[remaining] <- resid_var[remaining] *
                            (1 - (r_sg^2 / r_ss))

    # Update the remaining partial-correlation submatrix
    partial_R[remaining, remaining] <- partial_R[remaining, remaining] -
      outer(r_sg, r_sg) / r_ss

    # Stop when the median residual variance meets the threshold
    # Median residual variance <= 1 - var_threshold
    if (length(remaining) > 0) {
      median_resid <- median(resid_var[remaining])
      if (median_resid <= (1 - var_threshold)) break
    }
  }

  return(sort(selected))
}

# ---- MGAS approximation ----

effective_tests <- function(eigenvalues) {
  j <- length(eigenvalues)
  j - sum(pmax(eigenvalues - 1, 0))
}

poly6 <- function(x) {
  -0.3867*x^6 + 0.0021*x^5 - 0.1347*x^4 -
   0.0104*x^3 + 0.7276*x^2 + 0.0068*x
}

extended_simes_fast <- function(p_sorted, U_sub) {
  j_total <- length(p_sorted)
  if (j_total == 1) return(min(p_sorted[1], 1.0))

  # Compute eigenvalues once
  eig_full <- eigen(U_sub, symmetric = TRUE, only.values = TRUE)$values
  q_e <- effective_tests(eig_full)

  # Precompute q_ej by linear interpolation
  # Interpolate between q_e1 = 1 and q_e_j_total = q_e
  q_ej_vec <- 1 + (q_e - 1) * (seq_len(j_total) - 1) / (j_total - 1)
  q_ej_vec <- pmax(q_ej_vec, 1)

  # Calculate weighted P values
  p_weighted <- pmin((q_e / q_ej_vec) * p_sorted, 1.0)
  p_simes <- min(p_weighted)

  return(p_simes)
}

mgas_pvalue <- function(pmat, R_snp, R_phe, cutoff = 0.5) {
  m  <- nrow(pmat)
  n  <- ncol(pmat)
  mn <- m * n

  p_vec <- as.vector(pmat)
  ord   <- order(p_vec)

  # Construct the full correlation matrix once
  Kron_base <- kronecker(R_snp, R_phe)
  U_full    <- poly6(Kron_base)
  diag(U_full) <- 1
  U_full <- pmin(pmax(U_full, -1), 1)

  # ---- Define correlation blocks ----
  assigned <- rep(0L, mn)
  block_id <- 0L
  for (seed in ord) {
    if (assigned[seed] != 0L) next
    block_id <- block_id + 1L
    candidates <- which(abs(U_full[seed, ]) > cutoff & assigned == 0L)
    assigned[candidates] <- block_id
  }
  K <- block_id

  # ---- Within-block extended Simes approximation ----
  block_p    <- numeric(K)
  key_member <- integer(K)

  for (b in seq_len(K)) {
    mem   <- which(assigned == b)
    p_mem <- p_vec[mem]
    ord_b <- order(p_mem)
    mem_s <- mem[ord_b]
    p_s   <- p_mem[ord_b]
    U_b   <- U_full[mem_s, mem_s, drop = FALSE]
    block_p[b]    <- extended_simes_fast(p_s, U_b)
    key_member[b] <- mem_s[1]
  }

  # ---- Between-block extended Simes approximation ----
  ord_block  <- order(block_p)
  bp_sorted  <- block_p[ord_block]
  key_sorted <- key_member[ord_block]
  U_key      <- U_full[key_sorted, key_sorted, drop = FALSE]
  p_mgas     <- extended_simes_fast(bp_sorted, U_key)

  return(p_mgas)
}

# ---- Precompute gene-level LD eigendecomposition ----

precompute_LD <- function(X, explained_threshold) {
  R <- cor(X)
  eig <- eigen(R, symmetric = TRUE)

  ord     <- order(eig$values, decreasing = TRUE)
  values  <- eig$values[ord]
  vectors <- eig$vectors[, ord]
  posval  <- values[values > 1e-8]
  R_psd <- make_positive_definite(R)

  r <- which(cumsum(posval / sum(posval)) >= explained_threshold)[1]

  list(R = R, U_r = vectors[, 1:r], lambda_r = values[1:r], r = r, R_psd = R_psd)
}

# ====================
# Simulate one replicate for one gene
# ====================

run_one_sim <- function(i, X, n_traits, h, pi, rg, re, U_r, lambda_r, r, R, R_psd,shared_prop) {

  N <- nrow(X); M <- ncol(X)

  # ---- Select causal SNPs ----
  m_total <- ceiling(pi * M)

  m_shared <- ceiling(
    shared_prop * m_total
  )

  m_specific <- m_total - m_shared
# --------------------------------------------------
# 1. shared causal SNP
# --------------------------------------------------
  shared_idx <- sample(
    seq_len(M),
    m_shared
  )
# --------------------------------------------------
# 2. specific causal SNP
# --------------------------------------------------
  available_idx <- setdiff(
    seq_len(M),
    shared_idx
  )
  needed_specific <- n_traits * m_specific

  if (length(available_idx) < needed_specific) {
    stop(
      "Not enough SNPs for trait-specific causal variants. ",
      "Please decrease pi or n_traits."
    )
  }
  specific_idx <- sample(
    available_idx,
    needed_specific
  )
# --------------------------------------------------
# 3. Assign trait-specific causal SNPs
# --------------------------------------------------
  specific_idx_list <- split(
    specific_idx,
    rep(
      seq_len(n_traits),
      each = m_specific
    )
  )
# --------------------------------------------------
# 4. Combine shared and trait-specific causal SNPs for each trait
# --------------------------------------------------
  causal_idx_list <- vector(
    "list",
    n_traits
  )
  for (k in seq_len(n_traits)) {
    causal_idx_list[[k]] <- c(
      shared_idx,
      specific_idx_list[[k]]
    )
  }
  # ---- Genetic correlation matrix ----
  rg_matrix <- matrix(rg, n_traits, n_traits); diag(rg_matrix) <- 1
  sigmaG <- diag(sqrt(h)) %*% rg_matrix %*% diag(sqrt(h))
  sigmaG <- make_positive_definite(sigmaG)


  B <- matrix(0, nrow = M, ncol = n_traits)
  m = ceiling(pi * M)
  bmat <- mvrnorm(m , mu = rep(0, n_traits), Sigma = sigmaG)

  for (k in 1:n_traits) {
  idx <- causal_idx_list[[k]]
  B[idx, k] <-  bmat[,k] }

  G <- X %*% B

  # ---- Residual correlation matrix ----
  ve <- apply(G, 2, var) * (1 / h - 1)
  re_matrix <- matrix(re, n_traits, n_traits); diag(re_matrix) <- 1
  sigmaE <- diag(sqrt(ve)) %*% re_matrix %*% diag(sqrt(ve))
  sigmaE <- make_positive_definite(sigmaE)
  E <- mvrnorm(n = N, mu = rep(0, n_traits), Sigma = sigmaE)

  Y <- G + E
  Y <- scale(Y)

  # ---- Single-SNP GWAS ----
  # Store effect estimates, standard errors, Z statistics and P values
 beta_list <- list(); se_list <- list(); z_list <- list(); p_list <- list()
  for (k in 1:n_traits) {
    y_k <- Y[, k]
    gwas_k <- apply(X, 2, function(col) {
      fit <- summary(lm(y_k ~ col))
      c(beta = fit$coefficients[2,1],
        se   = fit$coefficients[2,2],
        z    = fit$coefficients[2,3],
        p    = fit$coefficients[2,4])
    })
    beta_list[[k]] <- gwas_k["beta", ]
    se_list[[k]]   <- gwas_k["se", ]
    z_list[[k]]    <- gwas_k["z", ]
    p_list[[k]]    <- gwas_k["p", ]
  }
  Z <- do.call(cbind, z_list)
  colnames(Z) <- paste0("Trait", 1:n_traits)
  rownames(Z) <- colnames(X)
  E_resid <- Y - G
  cor_e   <- cor(E_resid)
  cor_rg <- cor(G[,1], G[,2])

  # ============================================================
  # Whitening-based test: reuse the precomputed U_r, lambda_r and r
  # ============================================================
  Z_star <- crossprod(U_r, Z)
  Z_star <- sweep(Z_star, 1, sqrt(lambda_r), "/")

  T_single <- colSums(Z_star^2)
  p_single <- pchisq(T_single, df = r, lower.tail = FALSE)

  eig_T <- eigen(cor_e, symmetric = TRUE)
  V <- eig_T$vectors
  D <- eig_T$values

  Z_final <- Z_star %*% V
  Z_final <- sweep(Z_final, 2, sqrt(D), "/")

  T_multi  <- sum(Z_final^2)
  df_multi <- r * n_traits
  p_multi  <- pchisq(T_multi, df = df_multi, lower.tail = FALSE)

  # ---- MTaSPUsSet ----
  p_MTaSPUsSet <- tryCatch({
    corPhe <- cor(Z)
    res <- MTaSPUsSet(Zs = Z, corSNP = R, corPhe = corPhe,
                      pow = c(1,2,4,8), pow2 = c(1,2,4,8),
                      n.perm = 1000, Ps = FALSE)
    as.numeric(res["MTaSPUsSet"])
  }, error = function(e) NA)

  # ============================================================
  # MGAS
  # ============================================================
   P_mat <- do.call(rbind, p_list)
  rownames(P_mat) <- paste0("Trait", 1:n_traits)
  colnames(P_mat) <- colnames(X)

  p_MGAS <- tryCatch({
    R_phe <- cor(Z)
    mgas_pvalue(pmat = P_mat, R_snp = R, R_phe = R_phe)
  }, error = function(e) {
    message(sprintf("  MGAS failed sim %d: %s", i, conditionMessage(e)))
    NA
  })


  # ---- Return one simulation replicate ----
result_row <- data.frame(sim = i)
  for (k in 1:n_traits) result_row[[paste0("p", k)]] <- p_single[k]
  result_row$p_multi      <- p_multi
  result_row$p_MTaSPUsSet <- p_MTaSPUsSet
  result_row$p_MGAS       <- p_MGAS
  return(result_row)    # Return the replicate results
}

# ====================
# Simulate each gene and summarize results
# ====================

all_results <- list()

for (gene_name in gene_list) {
  cat(sprintf("Simulating gene: %s  (SNPs: %d)\n",
              gene_name, length(gene_snp_idx[[gene_name]])))

  idx    <- gene_snp_idx[[gene_name]]
  X_gene <- X_all[, idx, drop = FALSE]
  X_gene <- scale(X_gene, center = TRUE, scale = TRUE)

  ld <- precompute_LD(X_gene, explained_threshold)

  sim_list <- mclapply(1:n_sim, function(i) {
    set.seed(i + which(gene_list == gene_name) * 10000)
    tryCatch(
      run_one_sim(i, X_gene, n_traits, h, pi, rg, re,
                  U_r = ld$U_r, lambda_r = ld$lambda_r, r = ld$r, R = ld$R,R_psd = ld$R_psd,shared_prop ),
      error = function(e) {
      message(sprintf("Gene %s sim %d error: %s", gene_name, i, conditionMessage(e)))
      NULL
    }
    )
  }, mc.cores = n_cores)

  sim_list <- Filter(Negate(is.null), sim_list)
  if (length(sim_list) == 0) next

  results_gene <- do.call(rbind, sim_list)

  results_gene$gene <- gene_name
  all_results[[gene_name]] <- results_gene
}

# ====================
# Summarize rejection rates across genes
# ====================
if (length(all_results) == 0) stop("all_results is empty")

cat(sprintf("\nTotal genes: %d\n", length(all_results)))

# Calculate the rejection rate separately for each gene
gene_power_mat <- t(sapply(all_results, function(g) {
  c(setNames(sapply(1:n_traits, function(k) mean(g[[paste0("p",k)]] < alpha, na.rm=TRUE)), paste0("Trait",1:n_traits)),
    multi     = mean(g$p_multi < alpha, na.rm=TRUE),
    MTaSPUsSet= mean(g$p_MTaSPUsSet < alpha, na.rm=TRUE),
    MGAS      = mean(g$p_MGAS < alpha, na.rm=TRUE))
}))

# Mean rejection rate across genes and its gene-level standard error
power_all <- colMeans(gene_power_mat, na.rm=TRUE)
se_all <- apply(gene_power_mat, 2, function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x))))

# Extract single-trait and multi-trait summaries
power_single <- power_all[paste0("Trait",1:n_traits)]
se_single <- se_all[paste0("Trait",1:n_traits)]

multi_methods <- c("multi","MTaSPUsSet","MGAS")
power_multi <- power_all[multi_methods]; names(power_multi) <- multi_methods
se_multi <- se_all[multi_methods]; names(se_multi) <- multi_methods

cat("\nPower summary:\n")
print(data.frame(Method=names(power_all), Power=as.numeric(power_all), SE=as.numeric(se_all)))

# ---- Summary table ----
test_names   <- c(paste0("mBAT_trait", 1:n_traits),
                  "mBAT_multi", "MTaSPUsSet", "MGAS")
method_names <- c(paste0("Trait", 1:n_traits), "multi", "MTaSPUsSet", "MGAS")

power_summary <- data.frame(
  Test  = test_names,
  Power = as.numeric(power_all[method_names]),
  SE    = as.numeric(se_all[method_names])
)
# ---- Include simulation parameters in output filenames ----
file_params <- paste0("pi=", pi,
                      ",ntraits=", n_traits,
                      ",h=", paste(round(h, 4), collapse="_"),
                      ",rg=", rg,
                      ",re=", re,
                      ",shared_prop=", shared_prop,
                      ",ngenes=", length(all_results),
                      ",nsim=", n_sim,
                      ",n_cores=", n_cores)

# ---- Combine and save results ----
results_df <- do.call(rbind, all_results)

write.csv(results_df,
          file.path(out_dir, paste0("p_", file_params, ".csv")),
          row.names = FALSE)

write.csv(power_summary,
          file.path(out_dir, paste0("power_", file_params, ".csv")),
          row.names = FALSE)

cat("\n=== Done ===\n")
print(power_summary)
