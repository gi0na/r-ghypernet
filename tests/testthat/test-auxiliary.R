library(ghypernet)

test_that("check_specs returns right values for unipartite adj", {
  expect_equal(ghypernet:::check_specs.matrix(adj_karate), c('directed'=FALSE,'selfloops'=FALSE,'bipartite'=FALSE))
})

test_that("check_specs returns right values for bipartite incidence", {
  expect_equal(ghypernet:::check_specs.matrix(adj_karate[1:15,16:34]), c('directed'=TRUE,'selfloops'=FALSE,'bipartite'=TRUE))
})
