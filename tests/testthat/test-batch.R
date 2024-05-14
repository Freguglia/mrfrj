test_that("batching works", {
  set.seed(2)
  fam <- "oneeach"
  maximal <- mrf2d::mrfi(2)
  initial <- c(T,F,F,T,F,F)
  C <- 2
  theta <- mrfrj:::zero_array(maximal, fam, C)
  kernel_probs <- c(1,1,1,1,1)
  a <- mrfrj_batch(mrf2d::Z_potts, maximal,
                        initial, theta,
                        kernel_probs = kernel_probs,
                        family = fam,
                        verbose = interactive(),
                        id = "teste",
                        batch_size = 200, nbatches = 3)
  expect_type(a, "list")
  expect_s3_class(a, "mrfrj")
  expect_true(nrow(a$chain) > 0)
})
