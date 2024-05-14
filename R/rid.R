rid <- function(size = 8){
  paste0(sample(c(letters, 0:9), size = size, replace = TRUE), collapse = "")
}
