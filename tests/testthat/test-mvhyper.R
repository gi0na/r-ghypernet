test_that("dmvhyper_base returns 0 / -Inf for invalid x (support checks)", {
  n <- c(10L, 5L, 7L)
  k <- 8L
  
  x1 <- c(-1L, 2L, 7L)
  expect_identical(dmvhyper_base(x1, n, k, log = FALSE), 0)
  expect_identical(dmvhyper_base(x1, n, k, log = TRUE), -Inf)
  
  x2 <- c(11L, 0L, -3L)
  expect_identical(dmvhyper_base(x2, n, k, log = FALSE), 0)
  expect_identical(dmvhyper_base(x2, n, k, log = TRUE), -Inf)
  
  x3 <- c(1L, 1L, 1L) # sum=3, k=8
  expect_identical(dmvhyper_base(x3, n, k, log = FALSE), 0)
  expect_identical(dmvhyper_base(x3, n, k, log = TRUE), -Inf)
  
  expect_identical(dmvhyper_base(c(1L, 1L, 1L), n, k = 100L, log = TRUE), -Inf)
  expect_identical(dmvhyper_base(c(1L, 1L, 1L), n, k = 100L, log = FALSE), 0)
})

test_that("dmvhyper_base matches the closed-form formula on valid inputs", {
  n <- c(10L, 5L, 7L)
  x <- c(3L, 2L, 3L)
  k <- sum(x)
  
  lp_expected <- sum(lchoose(n, x)) - lchoose(sum(n), k)
  expect_equal(dmvhyper_base(x, n, k, log = TRUE), lp_expected, tolerance = 0)
  expect_equal(dmvhyper_base(x, n, k, log = FALSE), exp(lp_expected), tolerance = 0)
})

test_that("dmvhyper_base handles edge cases: k=0 and single-category", {
  # k = 0
  n <- c(4L, 9L, 2L)
  x <- c(0L, 0L, 0L)
  expect_equal(dmvhyper_base(x, n, k = 0L, log = TRUE), 0, tolerance = 0)
  expect_equal(dmvhyper_base(x, n, k = 0L, log = FALSE), 1, tolerance = 0)
  
  # single category => deterministic
  n1 <- 10L
  k1 <- 7L
  x1 <- 7L
  expect_equal(dmvhyper_base(x1, n1, k = k1, log = FALSE), 1, tolerance = 0)
  expect_equal(dmvhyper_base(x1, n1, k = k1, log = TRUE), 0, tolerance = 0)
  
  # invalid for single category
  expect_identical(dmvhyper_base(8L, n1, k = k1, log = FALSE), 0)
  expect_identical(dmvhyper_base(8L, n1, k = k1, log = TRUE), -Inf)
})

test_that("rmvhyper_base outputs valid samples (nonnegative, within n, correct row sums)", {
  set.seed(1)
  
  n <- c(5L, 10L, 20L, 0L, 7L)
  k <- 12L
  nn <- 500L
  
  x <- rmvhyper_base(nn = nn, n = n, k = k)
  
  expect_true(is.matrix(x))
  expect_identical(ncol(x), length(n))
  expect_identical(nrow(x), nn)
  
  expect_true(all(x >= 0))
  expect_true(all(t(t(x) <= n)))
  expect_true(all(rowSums(x) == k))
})

test_that("rmvhyper_base handles edge cases deterministically", {
  n <- c(5L, 10L, 3L)
  
  # nn = 0
  x0 <- rmvhyper_base(nn = 0L, n = n, k = 4L)
  expect_true(is.matrix(x0))
  expect_identical(nrow(x0), 0L)
  expect_identical(ncol(x0), length(n))
  
  # k = 0 => all zeros
  xk0 <- rmvhyper_base(nn = 10L, n = n, k = 0L)
  expect_true(all(xk0 == 0L))
  
  # single category => always k
  x1 <- rmvhyper_base(nn = 12L, n = 10L, k = 7L)
  expect_true(all(x1[, 1] == 7L))
  
  # k == sum(n) => always n
  xall <- rmvhyper_base(nn = 8L, n = n, k = sum(n))
  expect_true(all(xall == rep(n, each = 8L)))
})

test_that("rmvhyper_base marginal means are close to expectation", {
  set.seed(42)
  
  n <- c(15L, 25L, 10L, 50L)
  N <- sum(n)
  k <- 30L
  nn <- 40000L
  
  x <- rmvhyper_base(nn = nn, n = n, k = k)
  
  mu <- k * (n / N)
  # Hypergeom marginal variance
  var <- k * (n / N) * (1 - n / N) * ((N - k) / (N - 1))
  se <- sqrt(var / nn)
  
  expect_true(all(abs(colMeans(x) - mu) <= 6 * se + 0.02))
})
