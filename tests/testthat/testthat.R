library(testthat)

library(MVMT)

test_that("MVMT runs end-to-end and produces output file",{
  
  data_dir <- "/public/home/yuke/MVMT/inst/extdata"
  
  skip_if_not(dir.exists(data_dir), "Test data directory not found")
  
  bim_file    <- file.path(data_dir, "1000G_eur_unrel_ukbcom22.bim")
  map_file    <- file.path(data_dir, "chr22_gene.bed")
  assoc1      <- file.path(data_dir, "test1.ma")
  assoc2      <- file.path(data_dir, "test2.ma")
  LD_path     <- file.path(data_dir, "LD_Chr22")
  
  res <- MVMT(
    bim_file    = bim_file,
    map_file    = map_file,
    assoc_files = c(assoc1, assoc2),
    trait_names = c("HDL", "LDL"),
    LD_path     = LD_path,
    result_path = result_path,
    prop        = 0.9,
    n_cores     = 2
  )
  expect_false(is.null(res)) 
})