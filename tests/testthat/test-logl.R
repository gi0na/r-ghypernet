test_that("logl returns an error if all odds are zeros", {
  adj <- matrix(c(0,4,0,2,0,1,3,0,0),3,3)
  fit <- scm(adj,T,T)
  
  expect_error(ghype(adj, T,T, xi=fit$xi, omega =  fit$omega * 0))
})
