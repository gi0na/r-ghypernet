test_that("lr.test runs on small regular vs configuration models (Beta approximation)", {
  adj <- matrix(c(
    0, 1, 0, 0, 1,
    1, 0, 2, 0, 0,
    0, 2, 0, 1, 0,
    0, 0, 1, 0, 1,
    1, 0, 0, 1, 0
  ), nrow = 5, byrow = TRUE)
  
  directed <- FALSE
  selfloops <- FALSE
  
  nullmodel <- regularm(graph = adj, directed = directed, selfloops = selfloops)
  altmodel  <- scm(graph = adj, directed = directed, selfloops = selfloops)
  
  res <- lr.test(
    nullmodel = nullmodel,
    altmodel = altmodel,
    Beta = TRUE,
    nempirical = 20,   # keep tiny
    parallel = FALSE,
    seed = 1
  )
  
  expect_true(is.list(res))
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(is.finite(unname(res$statistic)))
  
  # df should match model degrees of freedom difference
  expect_true(is.finite(unname(res$parameter)))
  expect_identical(as.integer(unname(res$parameter)), as.integer(altmodel$df - nullmodel$df))
})

test_that("lr.test works with Beta parameters vector (numeric Beta shortcut)", {
  adj <- matrix(c(
    0, 1, 1,
    1, 0, 0,
    1, 0, 0
  ), nrow = 3, byrow = TRUE)
  
  nullmodel <- regularm(graph = adj, directed = FALSE, selfloops = FALSE)
  altmodel  <- scm(graph = adj, directed = FALSE, selfloops = FALSE)
  
  # Provide arbitrary but valid Beta params (a,b,mm) > 0
  Beta <- c(2, 5, 10)
  
  res <- lr.test(
    nullmodel = nullmodel,
    altmodel = altmodel,
    Beta = Beta,
    parallel = FALSE
  )
  
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(is.finite(unname(res$statistic)))
})

test_that("gof.test runs on a small fitted model and returns valid p-value", {
  adj <- matrix(c(
    0, 2, 0, 1,
    2, 0, 1, 0,
    0, 1, 0, 2,
    1, 0, 2, 0
  ), nrow = 4, byrow = TRUE)
  
  model <- scm(graph = adj, directed = FALSE, selfloops = FALSE)
  
  res <- gof.test(
    model = model,
    Beta = TRUE,
    nempirical = 20,
    parallel = FALSE,
    seed = 1
  )
  
  expect_true(is.list(res))
  expect_true(is.numeric(res$p.value))
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(is.finite(unname(res$statistic)))
})

test_that("link_significance returns a matrix of correct dimension with finite probabilities", {
  adj <- matrix(c(
    0, 1, 0, 0,
    1, 0, 2, 0,
    0, 2, 0, 1,
    0, 0, 1, 0
  ), nrow = 4, byrow = TRUE)
  
  model <- scm(graph = adj, directed = FALSE, selfloops = FALSE)
  
  # Over-representation p-values
  P_over <- link_significance(graph = adj, model = model, under = FALSE, log.p = FALSE)
  expect_true(is.matrix(P_over))
  expect_identical(dim(P_over), dim(adj))
  expect_true(all(is.finite(P_over)))
  expect_true(all(P_over >= 0 & P_over <= 1))
  
  # Under-representation p-values
  P_under <- link_significance(graph = adj, model = model, under = TRUE, log.p = FALSE)
  expect_true(is.matrix(P_under))
  expect_identical(dim(P_under), dim(adj))
  expect_true(all(is.finite(P_under)))
  expect_true(all(P_under >= 0 & P_under <= 1))
  
  # Log-probabilities should be <= 0 (since probabilities in (0,1])
  LP_over <- link_significance(graph = adj, model = model, under = FALSE, log.p = TRUE)
  expect_true(all(is.finite(LP_over)))
  expect_true(all(LP_over <= 0))
})

