# LISdmr R pacakge

R package for DMR detection combining MeDIP-seq and MRE-seq data through based on a conditional random field (CRF) model. The decteions are made from a mulitple testing procedure based on the local index of significance (LIS). 

## Installing

1. Install the [devtools](https://github.com/r-lib/devtools) package.

```
install.packages("devtools")
```

2. Load the devtools library.

```
library(devtools)
```

3. Install the [LISdmr](https://github.com/xiaoyudai/LISdmr) package from github.

```
install_github("xiaoyudai/LISdmr")
```

## Example use

An example to detect DMR between Brain cell and H1ES cell, based on the first 1000 CpGs at Chromesome 18. The example data is stored at [/extdata](https://github.com/xiaoyudai/LISdmr/tree/master/inst/extdata).

1. Load the [LISdmr](https://github.com/xiaoyudai/LISdmr) package.

```
library(LISdmr)
```

2. Specify the paths for the reference gene, MeDIP-seq and MRE-seq data for both samples (1 and 2).

```
filedir = system.file("extdata", package = "LISdmr")
ref_gene_path = paste(filedir, '/ref_gene.bed', sep = '')
medip1_path   = paste(filedir, '/brain_medip.bed', sep = '')
medip2_path   = paste(filedir, '/es_medip.bed', sep = '')
mre1_path     = paste(filedir, '/brain_mre.bed', sep = '')
mre2_path     = paste(filedir, '/es_mre.bed', sep = '')
```

3. Load the data and generate corresponding sequencing data features and genomic features.

```
dt = load_data(ref_gene_path, medip1_path, medip2_path, mre1_path, mre2_path)
```

4. Train the CRF model with a subset of available DMC labels.

```
train_path = paste(filedir, '/train_gene.bed', sep = '')
train(dt, train_path, c = 0.1, e = 0.005)
```

5. Predict the DMCs for the whole reference gene.

```
res_DMC = get_DMC(dt, alpha = 0.01)
```

6. Report DMRs by merging nearby DMCs.

```
res_DMR = get_DMR(res_DMC)
```

## Built With

* [GemoicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) - For converting MeDIP, MRE read counts and annotating genomic features.
* [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html) - For converting MeDIP, MRE read counts and annotating genomic features.
* [infotheo](https://cran.r-project.org/web/packages/infotheo/index.html) - For discretization of contiuous features in CRF.
* [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html) - For merging detected DMCs into DMRs.
* [crfpp](https://taku910.github.io/crfpp/) - For running CRF in modeling LIS values.

## Author

* **Xiaoyu Dai** - xiaoyu.dai@go.wustl.edu - Washington University in St. Louis

 
