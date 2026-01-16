library(ghypernet)

test_that("scm returns right df", {
  expect_equal(scm(adj_karate)$df,c('directed'=nrow(adj_karate)))
  expect_equal(scm(adj_karate[1:15,16:25])$df,c('directed'=25))
})

test_that("scm likelihood is computed correctly", {
  fit <- scm(adj_karate, directed = F, selfloops = T)
  ix <- mat2vec.ix(adj_karate, directed = F, selfloops = T)
  
  # lp_extradistr <- extraDistr::dmvhyper(x = adj_karate[ix], n = fit$xi[ix], k = fit$m, log = TRUE)
  # expect_equal(lp_extradistr, fit$loglikelihood)
  
  lp <- dmvhyper_base(x = adj_karate[ix], n = fit$xi[ix], log = TRUE)
  lp_manual <- sum(lchoose(fit$xi[ix],adj_karate[ix])) - lchoose(sum(fit$xi[ix]), fit$m)
  expect_equal(fit$loglikelihood, lp)
  expect_equal(lp_manual, lp)
})
