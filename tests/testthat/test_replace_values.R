testthat::test_that("Test replace_values",{
  
  df1 <- data.frame(col1=c(1,2,2,3),col2=c('a','b','a','c'))
  df2 <- data.frame(ID1=c('a','b'), ID2=c('c','d'))
  testthat::expect_equal(replace_values(df1,df2,2)$col2 == c("c","d","c",NA))
})
