dimension <- function(z, maximal_mrfi, family){
  length(mrf2d::smr_stat(z, maximal_mrfi, family))
}
