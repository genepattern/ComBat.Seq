## About ComBat-seq

ComBat-seq is a batch effect adjustment tool for bulk RNA-seq count data. It is an improved model based on the popular ComBat[1], to address its limitations through novel methods designed specifically for RNA-Seq studies. ComBat-seq takes **untransformed, raw count matrix** as input. Same as ComBat, it requires a known batch variable.

We use a negative binomial regression to model batch effects, then provide adjusted data by mapping the original data to an expected distribution if there were no batch effects. This approach better captures the properties of RNA-Seq count data compared to the Gaussian distribution assumed by ComBat. ComBat-seq specifies different dispersion parameters across batches, allowing for flexible modeling of the variance of gene expression. In addition, ComBat-seq provides adjusted data which preserves the integer nature of counts, so that the adjusted data are compatible with the assumptions of state-of-the-art differential expression software (e.g. edgeR, DESeq2, which specifically request untransformed count data). 


### Overview
ComBat-Seq is a batch effect adjustment tool for bulk RNA-seq count data. Improved model based on ComBat. Takes two different sets of data: 
  1. Combination of GCT + CLS file
<br> <br> or <br> <br>
  3. Combination of two .csv files, a raw count matrix and labels for samples.

### Documentation
  - ComBat-Seq author: Yuqing Zhang. Original GitHub Repo: https://github.com/zhangyuqing/ComBat-seq
  - Docker image used: genepattern/combat_seq:v1
  - Written mostly in R 4.2.1 


### Contact 
  - Edwin Huang
  - edh021@cloud.ucsd.edu

### Citation
Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078
