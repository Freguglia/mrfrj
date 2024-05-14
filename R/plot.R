position <- value <- t <- interaction <- NULL
variable <- name <- prob <- x <- y <- NULL

#' @import ggplot2
chain_plot <- function(df, thin = 5){
  df <- df[df$t %% thin == 0, ]
  ggplot(df, aes(x = t, y = value, color = position)) +
    facet_wrap(~interaction) +
    geom_point()
}

#' @rdname plot.mrfrj
#' @title Plotting `mrfrj` Objects.
#' @author Victor Freguglia
#'
#' @description
#' Plots the RJMCMC chain produced by the `mrfrj()` function.
#'
#' @param x A `mrfrj` object.
#' @param thin An integer representing the interval between observations in the
#' chain to be used. If `thin = 1`, all observations are used.
#' @param ... Unused.
#'
#' @details `ggplot2` is used to produce the plot, thus the method returns
#' a `ggplot` object which can be further modified with standard `ggplot2`
#' layers.
#'
#' @export
plot.mrfrj <- function(x, thin = 5, ...){
  chain_plot(x$chain, thin)
}

. <- list
#' @title Relative Positions Summary and Plot
#' @rdname rpos_summary
#' @importFrom dplyr as_tibble
#' @importFrom data.table setDT uniqueN
#'
#' @description
#' Computes the marginal proportions (estimates of the probability)
#' of each relative position in the pseudoposterior sample generated.
#' The results are returned in table format for \code{rpos_summary} and
#' as a plot for \code{rpos_plot}.
#'
#' @details
#' To speed up computations and reduce serial correlations, use
#' the \code{thin} and \code{start_at} parameters to obtain proportions
#' based on a subsample of the chain.
#'
#' @param chain A `\code{mrfrj}` object.
#' @param start_at First observation in the chain to be considered in the
#' subsample.
#' @inheritParams plot.mrfrj
#'
#'
#' @export
rpos_summary <- function(chain, thin = 5, start_at = 0){
  df <- setDT(chain$chain)
  inter <- df$interaction[1]
  df <- df[interaction == inter & df$t >= start_at & df$t %% thin == 0, ]
  n <- as.numeric(df[,.(uniqueN(t))])
  return(as_tibble(df[,.(prob = .N/n), position][order(-prob)]))
}


#' @title Relative Position Plots
#' @rdname rpos_summary
#' @importFrom tidyr separate_wider_delim
#' @importFrom dplyr mutate
#' @export
rpos_plot <- function(chain, thin = 5, start_at = 0){
  smr <- rpos_summary(chain)
  smr <- smr |>
    separate_wider_delim(cols = "position", delim = ",",
                         names = c("x", "y"), cols_remove = FALSE) |>
    mutate(x = as.integer(gsub("\\(", "", x)),
           y = as.integer(gsub("\\)", "", y)))
  smr <- rbind(smr, data.frame(x = 0L, y = 0L, position = "0,0", prob = NA))
  ggplot(smr, aes(x = x, y = y, fill = prob)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = "red", na.value = "black") +
    theme_void()
}
