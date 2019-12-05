testthat::test_that("Test best_unique",{
  
  df <- data.frame(col1=c('a','a', 'a'), evalue=c(1,2,1), percent_identity=c(94,99,94))
  testthat::expect_equal(best_unique(df, 'col1'), data.frame(col1='a',evalue=1,percent_identity=94))
})