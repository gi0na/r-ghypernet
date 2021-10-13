library(ghypernet)

test_that("scm returns right df", {
  expect_equal(scm(adj_karate)$df,c('directed'=nrow(adj_karate)))
  expect_equal(scm(adj_karate[1:15,16:25])$df,c('directed'=25))
})
