####################################
####################################
####################################
####################################
####################################
################### pvalues ##########
####################################

test_that("link significance returns correct probability for finding exactly or more edges and mhd case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = TRUE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = FALSE)) + 
           mapply(FUN = stats::dhyper, x = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                  MoreArgs = list(k = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or more edges and wallenius case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = TRUE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = FALSE)) + 
                       mapply(FUN = BiasedUrn::dWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                              odds = fit$omega/omegabar, MoreArgs = list(n = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or more edges and binomial approx case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = TRUE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                           MoreArgs = list(size = fit$m,lower.tail = FALSE)) + 
                       mapply(FUN = stats::dbinom, x = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega),
                              MoreArgs = list(size = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

################

test_that("link significance returns correct probability for finding exactly or less edges and mhd case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = TRUE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or less edges and wallenius case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = TRUE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or less edges and binomial approx case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = TRUE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})


####################################
################### log.p ##########
####################################

test_that("link significance returns correct probability for finding exactly or more edges and mhd case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = TRUE)
  
  xi <- fit$xi
  manual_p <- log(matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = FALSE)) + 
                       mapply(FUN = stats::dhyper, x = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                              MoreArgs = list(k = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or more edges and wallenius case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = TRUE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- log(matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = FALSE)) + 
                       mapply(FUN = BiasedUrn::dWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                              odds = fit$omega/omegabar, MoreArgs = list(n = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or more edges and binomial approx case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = TRUE, binomial.approximation = TRUE)
  
  manual_p <- log(matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = FALSE)) + 
                       mapply(FUN = stats::dbinom, x = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega),
                              MoreArgs = list(size = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

################

test_that("link significance returns correct probability for finding exactly or less edges and mhd case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = TRUE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = TRUE, log = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or less edges and wallenius case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = TRUE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- log(matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = TRUE)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding exactly or less edges and binomial approx case, log value", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = TRUE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = TRUE, log = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

####################################
####################################
####################################
####################################
####################################
################### not pvalues ##########
####################################

test_that("link significance returns correct probability for finding strictly more edges and mhd case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = FALSE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = FALSE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding strictly more edges and wallenius case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = FALSE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = FALSE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding strictly more edges and binomial approx case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = FALSE, give_pvals = FALSE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = FALSE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

################

test_that("link significance returns correct probability for finding strictly less edges and mhd case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = FALSE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = TRUE)) -
                       mapply(FUN = stats::dhyper, x = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                              MoreArgs = list(k = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding stricly less edges and wallenius case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = FALSE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = TRUE)) -
                       mapply(FUN = BiasedUrn::dWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                              odds = fit$omega/omegabar, MoreArgs = list(n = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding stricly less edges and binomial approx case", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = FALSE, give_pvals = FALSE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = TRUE)) -
                       mapply(FUN = stats::dbinom, x = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                              MoreArgs = list(size = fit$m)),3,3)
  
  expect_equal(pvals_over, manual_p)
})


####################################
################### log.p ##########
####################################

test_that("link significance returns correct probability for finding strictly more edges and mhd case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = FALSE)
  
  xi <- fit$xi
  manual_p <- matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = FALSE, log=TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding strictly more edges and wallenius case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = FALSE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- log(matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = FALSE)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding strictly more edges and binomial approx case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = FALSE, log.p = TRUE, give_pvals = FALSE, binomial.approximation = TRUE)
  
  manual_p <- matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = FALSE, log = TRUE)),3,3)
  
  expect_equal(pvals_over, manual_p)
})

################

test_that("link significance returns correct probability for finding strictly less edges and mhd case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- scm(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = FALSE)
  
  xi <- fit$xi
  manual_p <- log(matrix(mapply(FUN = stats::phyper, q = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                            MoreArgs = list(k = fit$m, lower.tail = TRUE)) -
                       mapply(FUN = stats::dhyper, x = as.vector(adj), m = as.vector(xi), n = fit$m^2 - as.vector(xi), 
                              MoreArgs = list(k = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding stricly less edges and wallenius case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = FALSE)
  
  xi <- fit$xi
  xibar <- fit$m^2 - xi
  omegabar <- (sum(fit$xi*fit$omega)-fit$xi*fit$omega)/xibar
  
  manual_p <- log(matrix(mapply(FUN = BiasedUrn::pWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                            odds = fit$omega/omegabar, MoreArgs = list(n = fit$m,lower.tail = TRUE)) -
                       mapply(FUN = BiasedUrn::dWNCHypergeo, x = as.vector(adj), m1 = as.vector(xi), m2 = fit$m^2 - as.vector(xi), 
                              odds = fit$omega/omegabar, MoreArgs = list(n = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})

test_that("link significance returns correct probability for finding stricly less edges and binomial approx case, log", {
  adj <- matrix(c(0,4,0,2,0,0,3,0,0),3,3)
  fit <- ghype(adj,T,T)
  pvals_over <- link_significance(graph = adj, model = fit, under = TRUE, log.p = TRUE, give_pvals = FALSE, binomial.approximation = TRUE)
  
  manual_p <- log(matrix(mapply(FUN = stats::pbinom, q = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                            MoreArgs = list(size = fit$m,lower.tail = TRUE)) -
                       mapply(FUN = stats::dbinom, x = as.vector(adj), prob = fit$xi*fit$omega/sum(fit$xi*fit$omega), 
                              MoreArgs = list(size = fit$m)),3,3))
  
  expect_equal(pvals_over, manual_p)
})
