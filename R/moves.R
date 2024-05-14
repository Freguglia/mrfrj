#' @importFrom mrf2d pl_mrf2d expand_array
move_within <- function(current_theta, current_pl, z, fdim, sdkernel, sdprior,
                        family, maximal_mrfi, C, included){
  included_theta <- rep(included, each = fdim / length(maximal_mrfi))
  proposed_theta <- current_theta + rnorm(fdim, mean = 0, sd = sdkernel)*(included_theta)
  proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
  proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(included)],
                          proposed_theta_arr[,,included, drop = FALSE])
  logA <- proposed_pl +
    sum(dnorm(proposed_theta[!is.na(proposed_theta)], sd = sdprior, log = TRUE)) -
    current_pl -
    sum(dnorm(current_theta[!is.na(current_theta)], sd = sdprior, log = TRUE))

  if(log(runif(1)) < logA){
    return(list(theta = proposed_theta, pl = proposed_pl,
                logA = logA, accept = TRUE))
  } else {
    return(list(theta = current_theta, pl = current_pl,
                logA = logA, accept = FALSE))
  }
}

move_swap <- function(current_theta, current_pl, z, included, fdim,
                      dim_per_group, family, maximal_mrfi, C){
  if(sum(included) > 1 && sum(included) < length(included)){
    to_remove <- -sample(-which(included), 1)
    to_include <- -sample(-which(!included), 1)
    proposed_theta <- current_theta
    proposed_included <- included
    proposed_included[to_remove] <- FALSE
    proposed_included[to_include] <- TRUE

    idx_include <- ((to_include - 1)*dim_per_group + 1):(to_include*dim_per_group)
    idx_remove <- ((to_remove - 1)*dim_per_group + 1):(to_remove*dim_per_group)
    proposed_theta[idx_include] <- proposed_theta[idx_remove]
    proposed_theta[idx_remove] <- NA

    proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
    proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(proposed_included)],
                            proposed_theta_arr[,,proposed_included, drop = FALSE])


    logA <- proposed_pl -
      current_pl
    if(log(runif(1)) < logA){
      return(list(theta = proposed_theta, pl = proposed_pl,
                  included = proposed_included,
                  logA = logA, accept = TRUE))
    }
    else {
      return(list(theta = current_theta, pl = current_pl,
                  included = included,
                  logA = logA, accept = FALSE))
    }
  } else {
    logA <- NA
    return(list(theta = current_theta, pl = current_pl, logA = logA,
                included = included,
                accept = NA))
  }
}

move_death <- function(current_theta, current_pl, z, dim_per_group, sdbirth,
                       family, maximal_mrfi, C, included, sdprior,
                       logpenalty_prior, kernel_probs){
  if(sum(included) > 0){
    pos_to_remove <- -sample(-which(included), size = 1)
    idx_remove <- ((pos_to_remove-1)*dim_per_group + 1):(pos_to_remove*dim_per_group)
    theta_to_remove <- current_theta[idx_remove]
    proposed_included <- as.logical(included - (seq_len(length(included)) == pos_to_remove))

    proposed_theta <- current_theta
    proposed_theta[idx_remove] <- NA

    if(sum(proposed_included) > 0){
      wts <- rgamma(sum(proposed_included), shape = 1/10)
      wts <- wts/sum(wts)
      mult <- rep(proposed_included, each = dim_per_group)
      mult[mult] <- mult[mult]*rep(wts, each = dim_per_group)
      proposed_theta <- proposed_theta + theta_to_remove*mult
    }

    proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
    proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(proposed_included)],
                            proposed_theta_arr[,,proposed_included, drop = FALSE])

    logA <- proposed_pl +
      sum(dnorm(proposed_theta[rep(proposed_included, each = dim_per_group)], sd = sdprior, log = TRUE)) -
      current_pl -
      sum(dnorm(current_theta[rep(included, each = dim_per_group)], sd = sdprior, log = TRUE)) +
      logpenalty_prior+
      sum(dnorm(theta_to_remove, sd = sdbirth, log = TRUE)) +
      log(kernel_probs[4] * 1/(length(included)-sum(included)+1)) -
      log(kernel_probs[3] * 1/sum(included))

    if(log(runif(1)) < logA) {
      return(list(theta = proposed_theta, pl = proposed_pl,
                  included = proposed_included,
                  logA = logA, accept = TRUE))
    } else {
      return(list(theta = current_theta, pl = current_pl,
                  included = included,
                  logA = logA, accept = FALSE))
    }
  } else {
    return(list(theta = current_theta, pl = current_pl,
                included = included,
                logA = -Inf, accept = NA))
  }
}

move_birth <- function(current_theta, current_pl, z, dim_per_group, sdbirth,
                       family, maximal_mrfi, C, included, sdprior,
                       logpenalty_prior, kernel_probs){
  if(sum(included) < length(included)){
    to_include <- -sample(-which(!included), size = 1)
    idx_new <- ((to_include-1)*dim_per_group + 1):(to_include*dim_per_group)
    theta_to_add <- rnorm(dim_per_group, mean = 0, sd = sdbirth)

    proposed_included <- included
    proposed_included[to_include] <- TRUE

    proposed_theta <- current_theta
    proposed_theta[idx_new] <- theta_to_add

    wts <- rgamma(sum(included), shape = 1/10)
    wts <- wts/sum(wts)
    mult <- rep(included, each = dim_per_group)
    mult[mult] <- mult[mult]*rep(wts, each = dim_per_group)
    proposed_theta <- proposed_theta - theta_to_add*mult

    proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
    proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(proposed_included)],
                            proposed_theta_arr[,,proposed_included, drop = FALSE])

    logA <- proposed_pl +
      sum(dnorm(proposed_theta[rep(proposed_included, each = dim_per_group)], sd = sdprior, log = TRUE)) -
      current_pl -
      sum(dnorm(current_theta[rep(included, each = dim_per_group)], sd = sdprior, log = TRUE)) -
      logpenalty_prior -
      sum(dnorm(theta_to_add, sd = sdbirth, log = TRUE)) +
      log(kernel_probs[3] * 1/(sum(included) + 1)) -
      log(kernel_probs[4] * 1/(sum(!included)))

    if(log(runif(1)) < logA) {
      return(list(theta = proposed_theta, pl = proposed_pl,
                  included = proposed_included,
                  logA = logA, accept = TRUE))
    } else {
      return(list(theta = current_theta, pl = current_pl,
                  included = included,
                  logA = logA, accept = FALSE))
    }
  } else { # Could not add more positions, full model
    return(list(theta = current_theta, pl = current_pl,
                included = included,
                logA = -Inf, accept = NA))
  }
}

move_jump <- function(included, current_theta, current_pl, z,
                      dim_per_group, maximal_mrfi, C, family,
                      sdprior, sdbirth, logpenalty_prior, fdim){
  jump_pos <- sample(seq_along(included), size = 1)
  vec_jump <- seq_along(included) == jump_pos
  proposed_included <- included
  proposed_included[jump_pos] <- !proposed_included[jump_pos]

  if(included[jump_pos]){ # Delete the selected position
    proposed_theta <- ifelse(rep(included*(vec_jump), each = dim_per_group), NA, current_theta)
    proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
    proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(proposed_included)],
                            proposed_theta_arr[,,proposed_included, drop = FALSE])

    logA <- proposed_pl +
      sum(dnorm(proposed_theta[rep(proposed_included, each = dim_per_group)],
                sd = sdprior, log = TRUE)) -
      current_pl -
      sum(dnorm(current_theta[rep(included, each = dim_per_group)],
                sd = sdprior, log = TRUE)) +
      sum(dnorm(current_theta[rep(vec_jump, each = dim_per_group)], sd = sdbirth, log = TRUE)) +
      logpenalty_prior

    if(log(runif(1)) < logA){
      return(list(theta = proposed_theta, pl = proposed_pl,
                  included = proposed_included,
                  logA = logA, accept = TRUE))
    } else {
      return(list(theta = current_theta, pl = current_pl,
                  included = included,
                  logA = logA, accept = FALSE))
    }

  } else { # Add the selected position
    proposed_theta <- ifelse(rep(vec_jump, each = dim_per_group),
                             rnorm(fdim, mean = 0, sd = sdbirth),
                             current_theta)
    proposed_theta_arr <- expand_array(proposed_theta, family, maximal_mrfi, C)
    proposed_pl <- pl_mrf2d(z, maximal_mrfi[which(proposed_included)],
                            proposed_theta_arr[,,proposed_included, drop = FALSE])

    logA <- proposed_pl +
      sum(dnorm(proposed_theta[rep(proposed_included, each = dim_per_group)],
                sd = sdprior, log = TRUE)) -
      current_pl -
      sum(dnorm(current_theta[rep(included, each = dim_per_group)], sd = sdprior, log = TRUE)) -
      sum(dnorm(proposed_theta[rep(vec_jump, each = dim_per_group)], sd = sdbirth, log = TRUE)) -
      logpenalty_prior

    if(log(runif(1)) < logA){
      return(list(theta = proposed_theta, pl = proposed_pl,
                  included = proposed_included,
                  logA = logA, accept = TRUE))
    } else {
      return(list(theta = current_theta, pl = current_pl,
                  included = included,
                  logA = logA, accept = FALSE))
    }
  }

}
