#' @description
#' Implements the Reversible Jump algorithm based on pseudoposteriors
#' from \insertCite{paper}{mrfrj},
#' used to generate a Markov Chain of subsets of parameters for
#' sparse relative position sets.
#' @references
#' \insertAllCited{}
#' @keywords internal
"_PACKAGE"

#' @title RJMCMC for MRFs
#' @rdname mrfrj
#' @description
#' Generates a Markov Chain of subsets of parameters of a Markov Random Field
#' model corresponding to
#' relative position sets, using the Reversible Jump algorithm
#' based on pseudoposteriors from \insertCite{paper}{mrfrj}.
#'
#' @param z The observed MRF realization.
#'   A `matrix` object containing values in \eqn{\{0, 1, \dots, C\}}.
#' @param maximal_mrfi The maximal RPS
#' (represented by a \code{\link[mrf2d]{mrfi}} object).
#' @param nsamples The number steps for the Markov chain to run.
#' @param initial_included A `logical` vector with length matching the
#' `maximal_mrfi` argument length indicating which positions are included
#' in the initial RPS (specified with `TRUE`).
#' @param initial_theta A `numerical` vector containing the initial values for
#' the interaction parameters \eqn{\theta}. See \code{\link[mrf2d]{smr_array}}
#' for more information on how the vector is expanded to the theta array.
#' @param sdprior Sample deviation considered for the prior distribution
#' of the interaction parameters \eqn{\theta}.
#' @param sdkernel The sample deviation of the proposal distribution for
#' the within-model random walk move proposal move.
#' @param sdbirth The sample deviation of the newly proposed theta values
#' in the birth, split and merge moves.
#' @param kernel_probs A length 5 `numeric` vector containing the probabilities
#' (or weights) of the Reversible Jump moves in the proposal kernel (in order):
#' \itemize{
#' \item{Within-model Random Walk.}
#' \item{Position Swap.}
#' \item{Death.}
#' \item{Birth.}
#' \item{Split and Merge.}
#' }
#' @param family The family of restrictions to be considered for the parameters.
#' For no restriction, use `family = "free"`. See \code{\link[mrf2d]{mrf2d-family}}.
#' @param logpenalty_prior A penalty to be considered in the prior distribution
#' of the RPS (in log-scale) whenever a new position is included.
#' `0.0` represents a uniform distribution over all RPSs that are a
#' subset of `maximal_mrfi`.
#' @param verbose A `logical` value to set whether the function progress
#' is printed in the screen.
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class `\code{mrfrj}`. The `\code{chain}`
#' attribute contains a data frame with all theta values in the chain,
#' identified by interaction, position and iteration. Positions not included
#' in the sampled vector for an iteration will not appear in the data.
#'
#' The `\code{duration}` attribute represents the total duration of the function
#' execution for both `\code{mrf_reversible_jump}` and `\code{mrfrj_batch}`.
#' The duration of a `mrfrj_batch` can be misleading if previous results
#' with the same id are reused. The `\code{batch_duration}` attribute
#' contains individual batch execution duration, which is always the
#' actual duration of the batch, even when it is reused.
#'
#' Results from `\code{mrf_reversible_jump}` also contains
#' `\code{last_theta}` and `\code{last_included}` that can be passed
#' as arguments for subsequent calls of the same functions.
#'
#' @import data.table
#' @importFrom tidyr unite
#' @importFrom dplyr pull
#' @importFrom stats dnorm rgamma rnorm runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
mrf_reversible_jump <- function(z, maximal_mrfi, nsamples,
                                initial_included, initial_theta,
                                sdprior = 1, sdkernel = 0.005,
                                sdbirth = 0.05,
                                kernel_probs = c(4,1,1,1,1),
                                family = "free",
                                logpenalty_prior = 0.0,
                                verbose = interactive()){

  start_time <- Sys.time()

  fdim <- dimension(z, maximal_mrfi, family)
  dim_per_group <- fdim / length(maximal_mrfi)
  C <- length(unique(as.vector(z))) - 1

  if(length(initial_included) != length(maximal_mrfi)){
    stop("initial_included must have the same length as maximal_mrfi.")
  }

  resmat <- matrix(0.0, nrow = nsamples, ncol = fdim)
  colnames(resmat) <- mrf2d::vec_description(maximal_mrfi, family, C) |>
    unite(name, position, interaction) |>
    pull(name)

  included <- initial_included
  current_theta <- initial_theta
  current_pl <- pl_mrf2d(z, maximal_mrfi[which(included)],
                         expand_array(current_theta, family, maximal_mrfi, C)[,,included, drop = FALSE])

  if(verbose & nsamples > 0) {
    pb <- txtProgressBar(0, nsamples, style = 2, char = "%",
                         width = getOption("width") - 15)
  }

  for(t in seq_len(nsamples)){

    if(verbose){setTxtProgressBar(pb, t)}
    move <- sample(moves, prob = kernel_probs, size = 1)

    if(move == "within") {
      move_result <- move_within(current_theta, current_pl, z, fdim, sdkernel,
                                 sdprior, family, maximal_mrfi, C, included)
      current_pl <- move_result$pl
      current_theta <- move_result$theta

    } else if(move == "swap") {
      move_result <- move_swap(current_theta, current_pl, z, included, fdim,
                               dim_per_group, family, maximal_mrfi, C)
      current_pl <- move_result$pl
      current_theta <- move_result$theta
      included <- move_result$included

    } else if(move == "death") {
      move_result <- move_death(current_theta, current_pl, z, dim_per_group, sdbirth,
                                family, maximal_mrfi, C, included, sdprior,
                                logpenalty_prior, kernel_probs)
      current_pl <- move_result$pl
      current_theta <- move_result$theta
      included <- move_result$included
    } else if(move == "birth") {
      move_result <- move_birth(current_theta, current_pl, z, dim_per_group, sdbirth,
                                family, maximal_mrfi, C, included, sdprior,
                                logpenalty_prior, kernel_probs)
      current_pl <- move_result$pl
      current_theta <- move_result$theta
      included <- move_result$included

    } else if(move == "jump") {
      move_result <- move_jump(included, current_theta, current_pl, z,
                               dim_per_group, maximal_mrfi, C, family,
                               sdprior, sdbirth, logpenalty_prior, fdim)
      current_pl <- move_result$pl
      current_theta <- move_result$theta
      included <- move_result$included
    } else stop()

    resmat[t,] <- current_theta

  }

  resdf <- table_iterations(resmat, iters = seq_len(nsamples))
  gc()

  end_time <- Sys.time()

  output <- list(chain = resdf,
                 duration = end_time - start_time,
                 last_theta = current_theta,
                 last_included = included)

  class(output) <- "mrfrj"

  return(output)
}
