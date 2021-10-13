library(ghypernet)

test_that("bccm with identical labels works", {
  m <- bccm(adj_karate, labels = rep(1,nrow(adj_karate)))
  expect_equal(class(m), 'ghype')
})

test_that("bccm works with bipartite graphs", {
  b <- adj_karate[sample(34,15),sample(34,19)]
  m <- bccm(adj = b, labels = list(vertexlabels[sample(34,15)],vertexlabels[sample(34,19)]))
  expect_equal(class(m), c("bccm","ghypeBlock","ghype"))
})
