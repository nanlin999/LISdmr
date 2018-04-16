

local_smooth <- function(x, window_size = 5){
  res = filter(x, rep(1/window_size,window_size ), sides=2)
  idx = which(is.na(res))
  res[idx] = x[idx]
  return(res)
}
