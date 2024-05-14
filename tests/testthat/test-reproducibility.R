test_that("reversible jump chains are reproducible", {
  fam <- "free"
  maximal <- mrf2d::mrfi(3)
  set.seed(1)
  initial <- sample(c(T,F), size = length(maximal), replace = T)
  C <- 2
  theta <- mrfrj:::zero_array(maximal, fam, C)
  kernel_probs <- c(1,1,1,1,1)
  set.seed(1)
  a <- mrf_reversible_jump(z = mrf2d::Z_potts,
                           initial_theta = theta,
                           maximal_mrfi = maximal,
                           nsamples = 200,
                           initial_included = initial,
                           kernel_probs = kernel_probs,
                           family = fam)
  set.seed(1)
  b <- mrf_reversible_jump(z = mrf2d::Z_potts,
                           initial_theta = theta,
                           maximal_mrfi = maximal,
                           nsamples = 200,
                           initial_included = initial,
                           kernel_probs = kernel_probs,
                           family = fam)
  d <- mrf_reversible_jump(z = mrf2d::Z_potts,
                           initial_theta = theta,
                           maximal_mrfi = maximal,
                           nsamples = 200,
                           initial_included = initial,
                           kernel_probs = kernel_probs,
                           family = fam)
  expect_identical(a$last_theta, b$last_theta)
  expect_false(identical(a$last_theta, d$last_theta))
})

test_that("reversible jump chains in batches are reproducible", {
  fam <- "free"
  z <- mrf2d::Z_potts
  maximal_mrfi <- mrf2d::mrfi(3)
  set.seed(1)
  initial_included <- sample(c(T,F), size = length(maximal_mrfi), replace = T)
  C <- length(unique(as.vector(z))) - 1
  initial_theta <- mrfrj:::zero_array(maximal_mrfi, fam, C)
  a <- mrfrj_batch(z, maximal_mrfi,
                   initial_included, initial_theta,
                   family = fam,
                   logpenalty_prior = 0,
                   id = "abc",
                   batch_size = 100, nbatches = 3,
                   seeds = 1:3)
  b <- mrfrj_batch(z, maximal_mrfi,
                   initial_included, initial_theta,
                   family = fam,
                   logpenalty_prior = 0,
                   id = "def",
                   batch_size = 100, nbatches = 3,
                   seeds = 1:3)
  expect_identical(a$chain, b$chain)
})

