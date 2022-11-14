## About ComBat-seq

### Overview
ComBat-Seq is a batch effect adjustment tool for bulk RNA-seq count data. Improved model based on ComBat. There are two different ways to run ComBat-Seq: 
  1. Combination of GCT + CLS file <br> 
    (Guide to GCT and CLS formats here:<br> https://www.genepattern.org/file-formats-guide#_Creating_Input_Files_GCT) <br> 
     Example dataset here:https://github.com/genepattern/ComBat-Seq/tree/develop/data/test_data_set_1
<br> <br> or <br> <br>
  2. A raw count matrix and labels for samples in .TSV (tab-separated values) format. <br> Example dataset here: https://github.com/genepattern/ComBat-Seq/tree/develop/data/test_data_set_2
  

### Documentation
  - ComBat-Seq author: Yuqing Zhang. Original GitHub Repo: https://github.com/zhangyuqing/ComBat-seq
  - Docker image used: genepattern/combat_seq:v1



### Citation
Yuqing Zhang, Giovanni Parmigiani, W Evan Johnson, ComBat-seq: batch effect adjustment for RNA-seq count data, NAR Genomics and Bioinformatics, Volume 2, Issue 3, 1 September 2020, lqaa078, https://doi.org/10.1093/nargab/lqaa078


### Contact
  - Edwin Huang: edh021@cloud.ucsd.edu



| Samples      | Sample 1 Name | Sample 2 Name | Sample 3 Name | Sample 4 Name | Sample 5 Name | ... |
| ----------- | ----------- | ----------- | ----------- | ----------- | ----------- | ----------- |
| Batch      | 1       | 1       | 2       | 2      | 2     |  ... |
| Group   | Group1 info        | Group1 info        | Group1 info        | Group1 info        | Group1 info        |  ... |
| ... | ... | ... | ... | ... | ... | ... |
