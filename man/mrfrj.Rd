% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrfrj.R, R/mrfrj_batch.R
\name{mrf_reversible_jump}
\alias{mrf_reversible_jump}
\alias{mrfrj_batch}
\title{RJMCMC for MRFs}
\usage{
mrf_reversible_jump(
  z,
  maximal_mrfi,
  nsamples,
  initial_included,
  initial_theta,
  sdprior = 1,
  sdkernel = 0.005,
  sdbirth = 0.05,
  kernel_probs = c(4, 1, 1, 1, 1),
  family = "free",
  logpenalty_prior = 0,
  verbose = interactive()
)

mrfrj_batch(
  z,
  maximal_mrfi,
  initial_included,
  initial_theta,
  sdprior = 1,
  sdkernel = 0.005,
  sdbirth = 0.05,
  kernel_probs = c(4, 1, 1, 1, 1),
  family = "free",
  logpenalty_prior = log(prod(dim(z))),
  verbose = interactive(),
  id = rid(),
  batch_size = 1000,
  nbatches = 3,
  seeds = seq_len(nbatches)
)
}
\arguments{
\item{z}{The observed MRF realization.
A `matrix` object containing values in \eqn{\{0, 1, \dots, C\}}.}

\item{maximal_mrfi}{The maximal RPS
(represented by a \code{\link[mrf2d]{mrfi}} object).}

\item{nsamples}{The number steps for the Markov chain to run.}

\item{initial_included}{A `logical` vector with length matching the
`maximal_mrfi` argument length indicating which positions are included
in the initial RPS (specified with `TRUE`).}

\item{initial_theta}{A `numerical` vector containing the initial values for
the interaction parameters \eqn{\theta}. See \code{\link[mrf2d]{smr_array}}
for more information on how the vector is expanded to the theta array.}

\item{sdprior}{Sample deviation considered for the prior distribution
of the interaction parameters \eqn{\theta}.}

\item{sdkernel}{The sample deviation of the proposal distribution for
the within-model random walk move proposal move.}

\item{sdbirth}{The sample deviation of the newly proposed theta values
in the birth, split and merge moves.}

\item{kernel_probs}{A length 5 `numeric` vector containing the probabilities
(or weights) of the Reversible Jump moves in the proposal kernel (in order):
\itemize{
\item{Within-model Random Walk.}
\item{Position Swap.}
\item{Death.}
\item{Birth.}
\item{Split and Merge.}
}}

\item{family}{The family of restrictions to be considered for the parameters.
For no restriction, use `family = "free"`. See \code{\link[mrf2d]{mrf2d-family}}.}

\item{logpenalty_prior}{A penalty to be considered in the prior distribution
of the RPS (in log-scale) whenever a new position is included.
`0.0` represents a uniform distribution over all RPSs that are a
subset of `maximal_mrfi`.}

\item{verbose}{A `logical` value to set whether the function progress
is printed in the screen.}

\item{id}{A string used as an internal identifier to save batch partial
results. If an execution of `\code{mrfrj_batch}` is interrupted, files
corresponding to batches with the same \code{id} will be reused.}

\item{batch_size}{Number of steps of the chain to run within each batch.}

\item{nbatches}{Number of batches to run.}

\item{seeds}{A vector of integers with the seed value to set with
\code{set.seed()} in each batch. Used to ensure reproducibility.}
}
\value{
An object of class `\code{mrfrj}`. The `\code{chain}`
attribute contains a data frame with all theta values in the chain,
identified by interaction, position and iteration. Positions not included
in the sampled vector for an iteration will not appear in the data.

The `\code{duration}` attribute represents the total duration of the function
execution for both `\code{mrf_reversible_jump}` and `\code{mrfrj_batch}`.
The duration of a `mrfrj_batch` can be misleading if previous results
with the same id are reused. The `\code{batch_duration}` attribute
contains individual batch execution duration, which is always the
actual duration of the batch, even when it is reused.

Results from `\code{mrf_reversible_jump}` also contains
`\code{last_theta}` and `\code{last_included}` that can be passed
as arguments for subsequent calls of the same functions.
}
\description{
Generates a Markov Chain of subsets of parameters of a Markov Random Field
model corresponding to
relative position sets, using the Reversible Jump algorithm
based on pseudoposteriors from \insertCite{paper}{mrfrj}.
}
\references{
\insertAllCited{}
}
