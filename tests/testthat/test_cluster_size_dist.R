testthat::test_that("Test single_cluster_df",{
  
  df1 <- data.frame(size=c(1,2), freq=c(10,20))
  df2 <- data.frame(size=c(1,2), freq=c(5,5))
  methods <- c("method1", "method2")
  ids <- c("id1", "id1")
  single_df <- single_cluster_df(list(df1,df2),methods,ids)
  testthat::expect_equal(ncol(single_df),4)
  testthat::expect_equal(nrow(single_df),4)
  testthat::expect_equal(single_df$size,c(1,2,1,2))
  testthat::expect_equal(single_df$freq,c(10,20,5,5))
  testthat::expect_equal(single_df$Method,c("method1","method1","method2","method2"))
  testthat::expect_equal(single_df$ID,rep("id1",4))
})
