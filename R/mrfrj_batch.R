#' @importFrom glue glue
#' @rdname mrfrj
#' @param id A string used as an internal identifier to save batch partial
#' results. If an execution of `\code{mrfrj_batch}` is interrupted, files
#' corresponding to batches with the same \code{id} will be reused.
#' @param batch_size Number of steps of the chain to run within each batch.
#' @param nbatches Number of batches to run.
#' @param seeds A vector of integers with the seed value to set with
#' \code{set.seed()} in each batch. Used to ensure reproducibility.
#' @export
mrfrj_batch <- function(z, maximal_mrfi,
                        initial_included, initial_theta,
                        sdprior = 1, sdkernel = 0.005,
                        sdbirth = 0.05,
                        kernel_probs = c(4,1,1,1,1),
                        family = "free",
                        logpenalty_prior = log(prod(dim(z))),
                        verbose = interactive(),
                        id = rid(),
                        batch_size = 1000, nbatches = 3,
                        seeds = seq_len(nbatches)){

  start_time <- Sys.time()
  last_included <- initial_included
  last_theta <- initial_theta


  for(batch in seq_len(nbatches)){

    if(verbose){cat(glue("\r{paste0(rep(' ', getOption('width')-15), collapse = '')}Batch {batch}/{nbatches}"))}

    if(!file.exists(glue("{id}-{batch}.rds"))){
      set.seed(seeds[batch])
      chain <- mrf_reversible_jump(z, maximal_mrfi, nsamples = batch_size,
                                   initial_included = last_included,
                                   initial_theta = last_theta,
                                   sdprior = sdprior, sdkernel = sdkernel,
                                   sdbirth = sdbirth,
                                   kernel_probs = kernel_probs,
                                   family = family,
                                   logpenalty_prior = logpenalty_prior,
                                   verbose = verbose)
      chain$chain$t <- chain$chain$t + (batch-1)* batch_size
      saveRDS(chain, file = glue("{id}-{batch}.rds"))
    } else {
      chain <- readRDS(glue("{id}-{batch}.rds"))
    }
    last_included <- chain$last_included
    last_theta <- chain$last_theta

  }

  result_list <- lapply(glue("{id}-{seq_len(nbatches)}.rds"),
                        readRDS)
  dfchain <- do.call(rbind, lapply(result_list, function(x) x$chain))
  clear <- lapply(glue("{id}-{seq_len(nbatches)}.rds"),
                  file.remove)

  end_time <- Sys.time()
  out <- list(chain = dfchain,
              duration = sapply(result_list, function(x) x$duration))
  list(chain = dfchain,
       batch_duration = lapply(result_list, function(x) x$duration),
       duration = end_time - start_time)

  class(out) <- c("mrfrj", "mrfjbatch")
  if(verbose){cat("\n")}

  return(out)
}
