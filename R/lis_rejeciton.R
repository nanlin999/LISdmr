lis_rejection <- function(lis, alpha = 0.01) {

  m=length(lis)
  sort_lis<-sort(lis)
  if (min(alpha) > alpha)
  {
    res = rep(0, m)
  }
  else
  {
    k = 1
    while(k<m && (1/k)*sum(sort_lis[1:k]) < alpha){
      k = k + 1
    }
    threshold = sort_lis[k - 1]
    rej_id = which(lis < threshold)
    res = rep(0, m)
    res[rej_id] = 1
  }
  return (res)
}
