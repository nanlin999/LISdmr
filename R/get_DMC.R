#' Predict DMC by the LIS-based multiple testing procedure, where LIS values are predicted from the trained CRF model.
#'
#' @param dt Data matrix resulted from load_data(), containing CRF features.
#' @param alpha Nominal FDR level.
#'
#' @return Data matrix containing raw MeDIP-seq, MRE-seq read counts and predicted DMC label.
#'
#' @export
get_DMC <- function(dt, alpha = 0.01) {
  if(!file.exists('/tmp/LISdmr/model')){
    warning('CRF model is not trained.')
  }

  write.table(dt[,8:33], '/tmp/LISdmr/test_dt', quote = F, sep = "\t", col.names = F, row.names = F)

  package_dir = system.file(package = "LISdmr")

  commend = paste(package_dir, '/src/crfpp/crf_test -v2 -m /tmp/LISdmr/model /tmp/LISdmr/test_dt > /tmp/LISdmr/prediction', sep = '')

  system(commend)

  result = read.table('/tmp/LISdmr/prediction')

  prob0 = result[,ncol(result)-1]
  prob0 = as.character(prob0)
  prob0 = as.numeric(substr(prob0,3,10))

  DMC = lis_rejection(prob0, alpha)

  res = dt[,1:7]
  res[,'DMC'] = DMC
  return(res)
}
