table_iterations <- function(resmat, iters){
  if(nrow(resmat) == 0){
    return(tibble::as_tibble(data.frame(t = numeric(0),
                                        value = numeric(0),
                                        position = character(0),
                                        interaction = character(0))))
  }
  resdf <- as.data.table(resmat)
  gc()
  resdf[,t:=iters]
  resdf <- data.table::melt(as.data.table(resdf), id.vars = "t")[order(t)][
      !is.na(value)][
        , c("position", "interaction") := tstrsplit(variable, "_", fixed=TRUE)][,variable:=NULL]
  resdf <- resdf[,interaction := as.factor(interaction)]
  resdf <- resdf[,position := as.factor(position)]
  gc()
  return(tibble::as_tibble(resdf))
}
