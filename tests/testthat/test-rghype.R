test_that("rghype returns an error if all odds are zeros", {
  adj <- matrix(c(0,4,0,2,0,1,3,0,0),3,3)
  fit <- scm(adj,T,T)
  fit$omega <- fit$omega * 0
  expect_error(rghype(10, fit))
})
