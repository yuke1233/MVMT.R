#' Map SNPs to Genes with Windowed Annotation
#'
#' Annotate single-nucleotide polymorphisms (SNPs) to genes based on their genomic
#' positions. This function extends each gene's boundaries upstream and downstream
#' by specified kilobase windows, then assigns SNPs falling within those extended
#' regions.
#'
#' @param info_snp Data frame of SNP annotations. Must contain columns:
#'   \code{snp}, \code{chr}, \code{pos}.
#' @param info_gene Data frame of gene boundaries. Must contain columns:
#'   \code{gene}, \code{chr}, \code{start}, \code{end}.
#' @param window_start Integer. Upstream extension in kilobases (default 50).
#' @param window_end Integer. Downstream extension in kilobases (default 50).
#' @param only_sets Logical. If \code{TRUE}, return only SNP sets per gene.
#'
#' @return A list with elements:
#' \describe{
#'   \item{sets}{Named list of SNP IDs per gene.}
#'   \item{map}{Data frame of SNP-to-gene mappings (if \code{only_sets = FALSE}).}
#' }
#'
#' @importFrom data.table data.table foverlaps setnames setDF setorder
#' @export
#'
#' @examples
#' \dontrun{
#' snps <- data.frame(
#'   snp = c("rs1", "rs2", "rs3"),
#'   chr = c("1", "1", "2"),
#'   pos = c(100000, 150000, 200000),
#'   stringsAsFactors = FALSE
#' )
#'
#' genes <- data.frame(
#'   gene  = c("GENE1", "GENE2"),
#'   chr   = c("1", "2"),
#'   start = c(120000, 190000),
#'   end   = c(140000, 210000),
#'   stringsAsFactors = FALSE
#' )
#'
#' result <- map_snp_to_gene(
#'   info_snp     = snps,
#'   info_gene    = genes,
#'   window_start = 20L,
#'   window_end   = 20L,
#'   only_sets    = FALSE
#' )
#' }
map_snp_to_gene <- function(info_snp,
                            info_gene,
                            window_start = 50L,
                            window_end = 50L,
                            only_sets = FALSE) {

  info_snp_name <- deparse(substitute(info_snp))

  df_check(info_snp)
  df_check(info_gene)
  columns_check(info_snp, c("snp", "chr", "pos"))
  columns_check(info_gene, c("gene", "chr", "start", "end"))
  nonneg_num_check(window_start, "window_start")
  nonneg_num_check(window_end, "window_end")

  if (anyDuplicated(info_snp$snp) > 0L) {
    stop("SNP IDs in ", "'", info_snp_name, "'", " must be unique.",
         call. = FALSE)
  }

  snpdat <- data.table::data.table(
    snp = as.character(info_snp$snp),
    chr = as.character(info_snp$chr),
    start = as.integer(info_snp$pos),
    end = as.integer(info_snp$pos),
    key = c("chr", "start", "end")
  )

  genedat <- data.table::data.table(
    gene = as.character(info_gene$gene),
    chr = as.character(info_gene$chr),
    start = as.integer(
      ifelse(info_gene$start - window_start * 1000 < 1L,
             1,
             info_gene$start - window_start * 1000)
    ),
    end = as.integer(info_gene$end + window_end * 1000),
    gene.start = as.integer(info_gene$start),
    gene.end = as.integer(info_gene$end),
    key = c("chr", "start", "end")
  )

  mapped <- data.table::foverlaps(snpdat, genedat, type = "within")

  snp_sets <- lapply(split(mapped[!is.na(gene), .(snp, gene)],
                           by = "gene", keep.by = FALSE),
                     function(x) unique(unname(unlist(x))))

  if (only_sets) {
    list(sets = snp_sets)
  } else {
    data.table::setnames(mapped,
                         old = c("start", "end", "i.start"),
                         new = c("gene.adj.start", "gene.adj.end", "pos"))

    vkeep <- c("snp", "chr", "pos",
               "gene", "gene.start", "gene.end",
               "gene.adj.start", "gene.adj.end")

    list(sets = snp_sets,
         map = data.table::setDF(
           data.table::setorder(mapped[, vkeep, with = FALSE],
                                chr, gene.start, gene.end, pos,
                                na.last = TRUE)
         ))
  }
}