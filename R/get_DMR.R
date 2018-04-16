#' Merge nearby DMCs to report DMRs, using heuristic from DSS method.
#'
#' @param dt Data matrix resulted from get_DMC(), containing raw MeDIP-seq, MRE-seq read counts and predicted DMC label.
#' @param minlen Minimum length (in basepairs) required for DMR. Default is 100 bps.
#' @param minlen In all DMRs, the percentage of CG sites with significant p-values (less than p.threshold) must be greater than this threshold. Default is 0.7. This is mainly used for correcting the effects of merging of nearby DMRs.
#' @param minCG Minimum number of CpG sites required for DMR. Default is 10.
#'
#' @return A matrix of reported DMRs, specifying positions, total length and number of CpG sites contained.


get_DMR <- function(dt, minlen = 100, pct.sig = 0.7, minCG = 10) {
  if (!('DMC' %in% colnames(dt))) {
    warning('DMC is not predicted yet.')
  }

  tmp_all <- dt[,c('chr','start')]
  tmp_all[,'mu1'] = 0
  tmp_all[,'mu2'] = 0
  tmp_all[,'diff'] = 0
  tmp_all[,'diff.se'] = 0
  tmp_all[,'stat'] = 0
  tmp_all[,'phi1'] = 0
  tmp_all[,'phi2'] = 0
  tmp_all[,'pval'] = 1 - dt$DMC
  tmp_all[,'fdr'] = 0
  colnames(tmp_all)[2] = 'pos'

  DMR = callDMR(tmp_all, delta=0, p.threshold=0.01, minlen, pct.sig, minCG)

  return(DMR[,1:5])
}
