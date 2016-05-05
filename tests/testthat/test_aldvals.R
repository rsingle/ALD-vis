library(asymLD)
context("Correct ALD values")

test_that("ALD values returned are correct", {
  
  dat <- data.frame(haplo.freq=c(0.3, 0.5, 0.2), locus1=rep('A', 3), locus2=rep('B', 3), 
    allele1=c('A1', 'A2', 'A2'), allele2=c('B1', 'B2', 'B3'))
  ald.results <- compute.ALD(dat)

  expect_equal(ald.results$ALD.1.2, 1)
  expect_equal(round(ald.results$ALD.2.1, 2), 0.73)
})
