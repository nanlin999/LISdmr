#' Load the data (MeDIP-seq, MRE-seq) and generate corresponding features.
#'
#' @param ref_gene_path Path of the reference gene file, which is a bed file with 3 columns (chr, start, end) specifying the CpG sites under analysis.
#' @param medip1_path Path of the MeDIP-seq data for the first sample.
#' @param medip2_path Path of the MeDIP-seq data for the second sample.
#' @param mre1_path Path of the MRE-seq data for the first sample.
#' @param mre2_path Path of the MRE-seq data for the second sample.
#'
#' @return Data matrix containing squencing data features and genomic features, ready to be fitted into CRF model.
#'
#' @export
load_data <- function(ref_gene_path, medip1_path, medip2_path, mre1_path, mre2_path) {
  print('Loading data and generate features. This may take a while...')

  n_cluster = 8000

  all = read.table(ref_gene_path, header = F)
  colnames(all) = c('chr', 'start', 'end')
  gr_ref_gene = GenomicRanges::GRanges(all[,1], ranges = IRanges::IRanges(start = all[,2], end = all[,3]))

  data_paths = c(medip1_path, medip2_path, mre1_path, mre2_path)
  names = c('medip1', 'medip2', 'mre1', 'mre2')
  for(i in 1:4) {
    raw_data = read.table(data_paths[i], header = F)
    if (i < 3) {
      gr_raw_data = GenomicRanges::GRanges(raw_data[,1], ranges = IRanges::IRanges(
        start = raw_data[,2], end = raw_data[,3]))
    } else {
      gr_raw_data = GenomicRanges::GRanges(raw_data[,1], ranges = IRanges::IRanges(
        start = raw_data[,2], end = raw_data[,2] + 2))
    }
    count = GenomicRanges::countOverlaps(gr_ref_gene, gr_raw_data)
    all[,names[i]] = count
  }

  for(window_size in c(1,10,50,200,1000)) {
    medip1 = local_smooth(all$medip1, window_size)
    medip2 = local_smooth(all$medip2, window_size)
    mre1 = local_smooth(all$mre1, window_size)
    mre2 = local_smooth(all$mre2, window_size)

    medip1 = medip1 / quantile(medip1[which(medip1>0)], 0.75) * 10
    medip2 = medip2 / quantile(medip2[which(medip2>0)], 0.75) * 10
    mre1 = mre1 / quantile(mre1[which(mre1>0)], 0.75) * 10
    mre2 = mre2 / quantile(mre2[which(mre2>0)], 0.75) * 10

    ratio_medip = abs(medip1 - medip2) / min(medip1 + 1, medip2 + 1)
    ratio_mre = abs(mre1 - mre2) / min(mre1 + 1, mre2 + 1)
    dif_combine = abs((medip1+1)*(mre2+1) - (medip2+1)*(mre1+1))

    ratio_medip = infotheo::discretize(ratio_medip, disc = 'equalfreq', nbins = n_cluster)
    ratio_mre = infotheo::discretize(ratio_mre, disc = 'equalfreq', nbins = n_cluster)
    dif_combine = infotheo::discretize(dif_combine, disc = 'equalfreq', nbins = n_cluster)

    all[,paste('ratio_medip_', window_size, sep='')] = ratio_medip
    all[,paste('ratio_mre_', window_size, sep='')] = ratio_mre
    all[,paste('dif_combine_', window_size, sep='')] = dif_combine
  }

  all[,'dist'] = all$start - c(all$start[1], all$start[1:nrow(all)-1])
  all[,'dist'] = infotheo::discretize(all$dist, disc = 'equalfreq', nbins = n_cluster)

  names = c('cpgi', 'exon', 'intron', 'utr3', 'utr5', 'rmsk_DNA', 'rmsk_LINE', 'rmsk_LowSimple', 'rmsk_SINE', 'rmsk_other')

  for(i in 1:10) {
    gr_genFeature = GenomicRanges::GRanges(genFeatures[[i]][,1], ranges = IRanges::IRanges(
      start = genFeatures[[i]][,2], end = genFeatures[[i]][,3]))
    count = GenomicRanges::countOverlaps(gr_ref_gene, gr_genFeature)
    all[,names[i]] = count
  }
  return(all)
}



