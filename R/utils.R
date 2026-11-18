# ========== globals.R ==========

#' Global variables for non-standard evaluation
#' @keywords internal
#' @noRd
utils::globalVariables(
  c(".", "snp", "chr", "pos", "A1", "A2",
    "gene", "gene.start", "gene.end",
    "pvalue", "chisq", "ID", "PhysPos", "A1Freq",
    "dif_AF", "i.start", "gene.adj.start", "gene.adj.end",
    "beta", "se", "P", "N", "AF1")
)

# ========== columns_check ==========

#' Check if data frame contains required columns
#' @param df A data frame.
#' @param columns Character vector of required column names.
#' @keywords internal
#' @noRd
columns_check <- function(df, columns) {
  df_name <- deparse(substitute(df))
  df_columns <- names(df)
  diff_columns <- setdiff(columns, df_columns)
  if (length(diff_columns) > 0) {
    stop("'", df_name, "' doesn't have column(s): ",
      paste(diff_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
}

# ========== df_check ==========

#' Check if object is a data frame
#' @param df Object to check.
#' @keywords internal
#' @noRd
df_check <- function(df) {
  df_name <- deparse(substitute(df))
  if (!is.data.frame(df)) {
    stop("'", df_name, "' is not a data frame.",
      call. = FALSE
    )
  }
}

# ========== nonneg_num_check ==========

#' Check if value is a non-negative number
#' @param x Value to check.
#' @param arg Name of the argument (for error message).
#' @keywords internal
#' @noRd
nonneg_num_check <- function(x, arg) {
  if (!is.numeric(x) || is.na(x) || x < 0L || length(x) != 1L) {
    stop("'", arg, "' must be a non-negative number of length 1.",
      call. = FALSE
    )
  }
}

# ========== num_between_check ==========

#' Check if value is within a numeric range
#' @param x Value to check.
#' @param left Lower bound.
#' @param right Upper bound.
#' @param arg Name of the argument (for error message).
#' @keywords internal
#' @noRd
num_between_check <- function(x, left, right, arg) {
  if (!is.numeric(x) || is.na(x) || x < left || x > right) {
    stop("'", arg, "' must be a number of length 1 between ",
      left, " and ", right, ".",
      call. = FALSE
    )
  }
}

# ========== integer_vector_check ==========

#' Check if values are whole numbers
#' @param x Numeric vector.
#' @param tol Tolerance for checking integer equality.
#' @keywords internal
#' @noRd
integer_vector_check <- function(x, tol = .Machine$double.eps) {
  is_wholenumber <- function(x, tol = tol) {
    abs(x - round(x)) < tol && !is.na(x)
  }
  tryCatch(all(sapply(x, is_wholenumber, tol = tol)),
    error = function(e) invisible(FALSE)
  )
}

# ========== pretty_num ==========

#' Format numbers with comma separators
#' @param x Numeric vector.
#' @param ... Additional arguments passed to \code{prettyNum}.
#' @keywords internal
#' @noRd
pretty_num <- function(x, ...) {
  prettyNum(x, big.mark = ",", scientific = FALSE, ...)
}

# ========== if_dup_ind ==========

#' Find duplicated indices
#' @param x Vector to check.
#' @keywords internal
#' @noRd
if_dup_ind <- function(x) {
  which(duplicated(x) | duplicated(x, fromLast = TRUE))
}