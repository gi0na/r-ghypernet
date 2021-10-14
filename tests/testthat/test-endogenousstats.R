library(ghypernet)

test_that("get_zero_dummy() works", {
  recip_stat <- reciprocity_stat(adj_karate)
  recip_stat_dummy <- get_zero_dummy(recip_stat)
  names(recip_stat_dummy) <- c('reciprocity', 'reciprocity_zeroes')
  expect_equal(all(recip_stat_dummy$reciprocity_zeroes[recip_stat==0]==exp(1)), TRUE)
})
