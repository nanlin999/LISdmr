#' Train the CRF model with subset of DMC labels.
#'
#' @param dt Data matrix resulted from load_data(), containing CRF features.
#' @param train_path Path of the training DMC labels, which is a bed file containg a subset of dt where DMC labels are available.
#' @param c Hyper-parameter for the CRF. With larger c value, CRF tends to overfit the training data.
#' @param e Hyper-parameter for the CRF specifiying the termination criterion.
#'


train <- function(dt, train_path, c = 0.1, e = 0.005) {
  train_label = read.table(train_path, header = F)
  colnames(train_label) = c('chr', 'start', 'end', 'DMC')
  train_dt = merge(dt, train_label, by = c('chr', 'start', 'end'), all.x = T)
  train_dt$chr = as.character(train_dt$chr)
  train_dt[which(is.na(train_dt$DMC)),] = ' '
  dir.create('/tmp/LISdmr/', showWarnings = F)
  write.table(train_dt[,8:34], '/tmp/LISdmr/train_dt', quote = F, sep = "\t", col.names = F, row.names = F)

  package_dir = system.file(package = "LISdmr")

  commend = paste(package_dir, '/src/crfpp/crf_learn -c ', c, ' -e ', e, ' ', package_dir, '/src/crfpp/template/template1 /tmp/LISdmr/train_dt /tmp/LISdmr/model', sep = '')
  system(commend)
  print('Trained model is saved as /tmp/LISdmr/model')
}


