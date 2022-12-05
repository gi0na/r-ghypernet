library(ghypernet)

test_that("bccm with identical labels works", {
  m <- bccm(adj_karate, labels = rep(1,nrow(adj_karate)))
  expect_equal(class(m), 'ghype')
})

test_that("bccm works with bipartite graphs", {
  b <- adj_karate[sample(34,15),sample(34,19)]
  m <- bccm(adj = b, labels = list(vertexlabels[sample(34,15)],vertexlabels[sample(34,19)]), multinomial = TRUE)
  expect_equal(class(m), c("bccm","ghypeBlock","ghype"))
})

test_that("bccm returns same matrix for blockOmega and Omega when every node is in a block of its own", {
  adj <- matrix(0, 3, 3)
  adj[1, 2] <- 1
  adj[2, 3] <- 1
  adj[3, 1] <- 1
  m <- bccm(adj = adj, labels = 3:1, directed = TRUE, selfloops = FALSE, directedBlocks = TRUE, regular = TRUE)
  expect_equal(all(m$blockOmega == m$omega), TRUE)
})

test_that("bccm runs through", {
  m <- bccm(adj = adj_karate, labels = vertexlabels, directed = F, selfloops = F, homophily = F)
  expect_equal(all(m$blockOmega %in% m$coef), TRUE)
})