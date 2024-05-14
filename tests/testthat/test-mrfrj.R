test_that("mrfrj works", {
  set.seed(2)
  fam <- "oneeach"
  maximal <- mrf2d::mrfi(2)
  initial <- c(T,F,F,T,F,F)
  C <- 2
  theta <- mrfrj:::zero_array(maximal, fam, C)
  kernel_probs <- c(1,1,1,1,1)
  a <- mrf_reversible_jump(z = mrf2d::Z_potts,
                           initial_theta = theta,
                           maximal_mrfi = maximal,
                           nsamples = 500,
                           initial_included = initial,
                           kernel_probs = kernel_probs,
                           family = fam)
  expect_type(a, "list")
  expect_s3_class(a, "mrfrj")
  expect_true(nrow(a$chain) > 0)
})
