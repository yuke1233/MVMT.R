# MVMT

**MVMT: Multi-Trait Gene-Based Association Test**

MVMT is an R package for gene-level association testing across multiple correlated traits using GWAS summary statistics and an external LD reference. It accounts for LD among variants and correlations among trait-specific association statistics, and returns a joint chi-square statistic and P value for each successfully analyzed gene.

## Installation

Install the dependencies and the downloaded source package in R:

```r
install.packages(c("data.table", "dplyr"))

install.packages(
  "MVMT_0.1.0.9000.tar.gz",
  repos = NULL,
  type = "source"
)

library(MVMT)
```

Replace the archive path with the location of your downloaded file. The `parallel` package is included with R.

## Usage

```r
library(MVMT)

results <- MVMT(
  bim_file    = "data/reference.bim",
  map_file    = "data/gene_coords.txt",
  assoc_files = c("data/HDL.ma", "data/LDL.ma", "data/logTG.ma"),
  trait_names = c("HDL", "LDL", "logTG"),
  LD_path     = "data/LD",
  result_path = "results",
  prop        = 0.9,
  n_cores     = 1,
  af_threshold = 0.2
)

head(results)
```

Provide at least two traits. The order of `trait_names` must match `assoc_files`. The current `MVMT()` function maps SNPs to gene boundaries extended by 50 kb on each side.

| Argument | Description | Default |
|---|---|---|
| `bim_file` | PLINK BIM annotation for the LD reference | Required |
| `map_file` | Gene coordinates without a header | Required |
| `assoc_files` | Character vector of GWAS summary-statistic file paths | Required |
| `trait_names` | Trait names in the same order as `assoc_files` | Required |
| `LD_path` | Directory containing per-gene LD files | Required |
| `result_path` | Output directory; created if necessary | Required |
| `prop` | Cumulative proportion of positive LD eigenvalue variance to retain | `0.9` |
| `n_cores` | Number of cores for gene-wise analysis | `1` |
| `af_threshold` | Maximum absolute allele-frequency difference after allele alignment; differences must be below this value | `0.2` |

### Parallel execution

On Linux or macOS, increase `n_cores`, for example to `10`. On Windows, use `n_cores = 1`; gene-wise parallelism uses `parallel::mclapply()`.

When running multiple gene workers, set BLAS threads to one before starting R to avoid nested parallelism:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
R
```

## Input files

### PLINK BIM file

A standard six-column PLINK `.bim` file without a header:

```text
CHR  SNP  CM  BP  A1  A2
```

The fourth column contains the physical position in base pairs. The line above describes the column order and should not be included as a header.

### Gene coordinates

A whitespace-delimited file without a header, with columns in this order:

```text
chromosome  start  end  gene_id
```

Gene identifiers must match the prefixes of the corresponding LD files.

### GWAS summary statistics

Use the same eight-column input layout as mBAT. Each file must have a header, followed by columns in this order:

```text
snp  A1  A2  AF1  beta  se  P  N
```

Here, `A1` is the effect allele, `A2` is the other allele, `AF1` is the frequency of `A1`, `beta` is the effect estimate for `A1`, `se` is its standard error, `P` is the association P value, and `N` is the sample size. The package interprets the first eight columns in this order.

Supply genome-wide summary statistics, rather than significant variants alone, because MVMT estimates cross-trait correlations using shared variants with `P > 0.01` in every trait. SNP identifiers and chromosome labels must be consistent across inputs, and all coordinates must use the same reference assembly. Harmonize GWAS allele orientation across traits before analysis; the package additionally aligns gene-level inputs to the LD reference.

### LD reference files

For each gene, provide:

```text
LD/<gene_id>.ldm.full.info
LD/<gene_id>.ldm.full.bin
```

The `.info` file must contain a header with `ID`, `PhysPos`, `A1`, `A2`, and `A1Freq`. The `.bin` file contains the full square LD correlation matrix as 32-bit floating-point values in column-major order, with rows and columns in the same SNP order as the `.info` file.

## Example

Place the supplied `extdata.zip` in your working directory and run:

```r
library(MVMT)

unzip("extdata.zip", exdir = ".")
example_dir <- "extdata"

results <- MVMT(
  bim_file = file.path(example_dir, "1000G_eur_unrel_ukbcom22.bim"),
  map_file = file.path(example_dir, "chr22_gene.bed"),
  assoc_files = file.path(example_dir, c("test1.ma", "test2.ma")),
  trait_names = c("Trait1", "Trait2"),
  LD_path = file.path(example_dir, "LD_Chr22"),
  result_path = "results",
  prop = 0.9,
  n_cores = 1,
  af_threshold = 0.2
)

head(results)
```

This chromosome 22 example returns results for 10 genes. Some annotated genes lack LD files and are skipped.

## Output

Results are returned as a data frame and written to a tab-separated file in `result_path`. The filename follows the trait order, for example `MVMT_HDL,LDL,logTG.txt` or `MVMT_Trait1,Trait2.txt`.

| Column | Description |
|---|---|
| `Gene` | Gene identifier |
| `Chr` | Chromosome |
| `Start`, `End` | Gene coordinates |
| `No.SNPs` | Number of SNPs included in the test |
| `No.Traits` | Number of traits analyzed jointly |
| `SNP_start`, `SNP_end` | First and last included SNP identifiers |
| `No.Eigenvalues` | Test degrees of freedom: retained LD components multiplied by the number of traits |
| `Chisq_mBAT_multi` | MVMT chi-square test statistic |
| `P_MVMT` | MVMT association P value |

Genes without valid LD data or sufficient SNPs after quality control are skipped. Skipped genes are not included in the result file.

## Help

```r
?MVMT
?MVMT_SVD
?map_snp_to_gene
?allele_flip_gwas_beta
```

## License

MIT. See the package `LICENSE` file.

## Acknowledgements

This implementation builds on the original [mBAT-combo framework](https://github.com/Share-AL-work/mBAT.combo) developed by Ang Li and colleagues.
