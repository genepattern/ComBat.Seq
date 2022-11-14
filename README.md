# About ComBat-Seq

## Overview
ComBat-Seq is a batch effect adjustment tool for bulk RNA-seq count data. Improved model based on ComBat. <br>
#### Parameters
To run ComBat Seq, these inputs are used: 
  1. input matrix (required)
      - The input matrix is a counts matrix with dimensions gene x sample. The input counts matrix. Can be a .GCT or a .tsv (tab separated value) file. Rows should be genes, column should be samples, and the value are counts. File format documentation for GCT files: https://www.genepattern.org/file-formats-guide Some sample inputs can be found here: https://github.com/genepattern/ComBat_Seq/tree/develop/data
  2. batch information (required)
      - Batch information contain information on batches. It is a table that looks like the following:
      - The file can be a .CLS file contaning only batch information, or a .tsv (tab separated file) containing batch, group, and any other additional information as long as it follows the format below. <br>

| Samples      | Sample 1 Name | Sample 2 Name | Sample 3 Name | Sample 4 Name | Sample 5 Name | ... |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| Batch      | 1       | 1       | 2       | 2      | 2     |  ... |
| Group   | Group1 info        | Group1 info        | Group1 info        | Group1 info        | Group1 info        |  ... |
| ... | ... | ... | ... | ... | ... | ... |
<br>
<br> 
      
  3. covariates (optional)
      - Row names for covariates to use for this run of ComBat Seq. 
  4. output prefix (required)
      - Prefix for output filenames. 
#### Advanced Parameters
  5. Shrink
      - Whether to apply empirical Bayes estimation on dispersion.
  6. gene subset n
      - Number of genes to use in empirical Bayes estimation, only useful when shrink = Yes
  7. covariance matrix
      - If you wish to specify multiple biological variables. Model matrix for other covariates to include in the linear model besides batch and condition of interest.

## Documentation
  - ComBat-Seq author: Yuqing Zhang. Original GitHub Repo: https://github.com/zhangyuqing/ComBat-seq
  - Docker image used: genepattern/combat_seq:v1



## Citation
Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078


## Contact
  - Edwin Huang: edh021@cloud.ucsd.edu




