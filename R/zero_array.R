zero_array <- function(mrfi, family, C){
  n_R <- length(mrfi)
  len <- switch (family,
    "free" = n_R*((C+1)^2 - 1),
    "onepar" = 1,
    "oneeach" = n_R,
    "absdif" = n_R*C,
    "dif" = n_R*2*C
  )
  return(rep(0.0, len))
}
