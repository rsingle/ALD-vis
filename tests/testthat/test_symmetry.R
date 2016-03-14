library(asymLD)
context("Symmetry")

test_that("ALD values are identical for symmetric bi-allelic data", {
  
  dat <- data.frame(haplo.freq=c(0.3, 0.5, 0.2), locus1=rep('A', 3), locus2=rep('B', 3), 
    allele1=c('A1', 'A2', 'A2'), allele2=c('B1', 'B2', 'B2'))
  ald.results <- compute.ALD(dat)
  
  expect_equal(ald.results$ALD.1.2, ald.results$ALD.2.1)
})
