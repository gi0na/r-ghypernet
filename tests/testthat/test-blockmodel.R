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
  # blockOmega labels are always sorted
  expect_equal(all(m$blockOmega == m$omega[3:1,3:1]), TRUE)
  
  m <- bccm(adj = adj, labels = 3:1, directed = TRUE, selfloops = FALSE, directedBlocks = FALSE, regular = TRUE)
  expect_equal(all(m$blockOmega == m$omega[3:1,3:1]), TRUE)
})

test_that("bccm runs through", {
  m <- bccm(adj = adj_karate, labels = vertexlabels, directed = F, selfloops = F, homophily = F)
  expect_equal(all(m$blockOmega %in% m$coef), TRUE)
  expect_equal(all(m$coef %in% m$blockOmega), TRUE)
  
  m <- bccm(adj = adj_karate, labels = vertexlabels, directed = F, selfloops = F, homophily = T)
  expect_equal(m$coef['homologue']==1, c(homologue=TRUE))
  
  m <- bccm(adj = adj_karate, labels = vertexlabels, directed = F, selfloops = F, homophily = F, inBlockOnly = TRUE)
  expect_equal(m$coef['1<->1']==1, c('1<->1'=TRUE))
})

test_that("bccm returns same matrix for blockOmega and Omega when every node is in a block of its own, again", {
labels <- c("B1", "C2", "D3", "D4", "B5", "C6")
adj <- matrix(0, 6, 6)

adj[2, 1] <- 1
adj[2, 6] <- 10
adj[3, 1] <- 2
adj[3, 6] <- 2
adj[4, 3] <- 5
adj[5, 4] <- 2
adj[6, 1] <- 2

m <- bccm(adj, labels, directed = TRUE, selfloops = FALSE, directedBlocks = TRUE)
expect_equal(all(m$blockOmega[labels,labels] == m$omega), TRUE)
})


test_that("  ", {
  labels <- c("B", "C", "D", "D", "B", "C")
  adj <- matrix(0, 6, 6)
  
  adj[2, 1] <- 1
  adj[2, 6] <- 10
  adj[3, 1] <- 2
  adj[3, 6] <- 2
  adj[4, 3] <- 5
  adj[5, 4] <- 2
  adj[6, 1] <- 2
  
  m <- bccm(adj, labels, directed = TRUE, selfloops = FALSE, directedBlocks = TRUE)
  expect_equal(all(m$blockOmega == m$omega[1:3,c(5,6,4)]), TRUE)
})