########################################################
# Parse Parameters
########################################################
print('================================================')
print("Loading library: optparse")
library("optparse")
source('/opt/genepatt/src/ComBat-seq_wrapper_helper.R')
library(tools)

# Parse input arguments
parser = OptionParser()
# parameter types: 'character', 'integer', 'logical', 'double', or 'complex'
# =================================================
# Input file
parser <- add_option(parser, c("--input_matrix"), help="Input count matrix file to load (containing the raw counts).", default = "NO_FILE")

parser <- add_option(parser, c("--input_class"), help="Input batch information file to load (containing the same phenotypes).", default = "NO_FILE")

parser <- add_option(parser, c("--covariates"), help = "Row names for covariates to use for this run of ComBat Seq. Separate by commas.", default = "None")

parser <- add_option(parser, c("--shrink"), help = "Whether to apply empirical Bayes estimation on parameters.", default = "No")

parser <- add_option(parser, c("--shrink_disp"), help = "Whether to apply empirical Bayes estimation on dispersion.", default = "No")

parser <- add_option(parser, c("--gene_n"),help = "Number of genes to use in empirical Bayes estimation, only useful when shrink = Yes", default = NULL)

parser <- add_option(parser, c("--cov_mat"),help = "If you wish to specify multiple biological variables.
Model matrix for other covariates to include in the linear model besides batch and condition of interest. ", default = "NO_FILE")

parser <- add_option(parser, "--raw_counts", help = "If only using raw count matrix", default = FALSE)
# ===================================================

# ===================================================
# parameter for save_it
parser <- add_option(parser, c("--file_name"), type='character',
                     default='combatseq_normalized', help="Basename of the
                     file to be saved.")
# ====================================================



print('================================================')
args <- parse_args(parser)
print('Parameters used:')
print(args)
print('================================================')

USE_GCT <<- FALSE
## Checking if the user uploaded a GCT file or Not
input_fp <- args$input_matrix
if ("gct" == file_ext(input_fp)){
  USE_GCT <<- TRUE
} else if (("csv" == file_ext(input_fp) || "tsv" == file_ext(input_fp) || "csv" == file_ext(input_fp))){
  ## parse using tsv parser
}

## Checking if the user uploaded a CLS file
USE_CLS <<- FALSE
if (grepl(".cls", args$input_class, fixed = TRUE)){
  USE_CLS <<- TRUE
}

## check for other arguments
SHRINK <<- FALSE
if (args$shrink == "Yes"){
  SHRINK <<- TRUE
  if (!is.null(args$gene_n)){
    GENE_N <<- args$gene_n
  }
}

SHRINK_DISP <- FALSE
if (args$shrink_disp == "Yes"){
  SHRINK_DISP <- TRUE
}

COV_MOD <<- NULL
if (args$covariates != "None"){
  COV_MOD <<- args$covariates
}

FULL_MOD <<- FALSE
if (args$covariates != "None"){
  FULL_MOD <<- TRUE
}



print(paste("GCT FILE IS PROVIDED? ", USE_GCT))
print(paste("GCT FILE IS PROVIDED? ", USE_CLS))


# Setting up the PDF file for the plots
# pdf(file=paste(args$file_name, '.pdf', sep=''))

##########################################################
# Function Definitions
##########################################################


# Run the wrapper
run_combat_seq <- function(args){
  ### Run ComBat_Seq based on what the user inputted.


  if (USE_GCT){
    data <- read_gct(args$input_matrix)
    batch_info <- read_multiline_batch_info_file(data$data, args$input_class, USE_CLS)

    adjusted <- as.data.frame(ComBat_seq(data$data,
      batch=batch_info$batches, group=batch_info$groups, shrink = SHRINK,
      shrink.disp = SHRINK_DISP, gene.subset.n = GENE_N, covar_mod = COV_MOD,
      full_mod = FULL_MOD))

    print('... done')
    save_data(paste(args$file_name, '.gct', sep=''), batch_corrected_matrix=adjusted,
        samples = data$samples, genes = data$genes, descriptions= data$descriptions, use_gct = USE_GCT)
  }else{
    data <- read_delimited_file(args$input_matrix)
    batch_info <- read_multiline_batch_info_file(data$data, args$input_class, USE_CLS)

    adjusted <- as.data.frame(ComBat_seq(data$data,
      batch=batch_info$batches, group=batch_info$groups, shrink = SHRINK,
      shrink.disp = SHRINK_DISP, gene.subset.n = GENE_N, covar_mod = COV_MOD,
      full_mod = FULL_MOD))
    print('... done')
    save_data(fileName = paste(args$file_name, ".tsv", sep = ""),
      batch_corrected_matrix = adjusted,
      samples = data$samples, genes = data$genes,use_gct = FALSE)
  }


}

############################################################
# Begin Running the functions
############################################################

## Loading Libraries
library("BiocManager")
library("edgeR")

source('/opt/genepatt/src/ComBat_seq.R')
source('/opt/genepatt/src/helper_seq.R')
print("ComBat Seq has been loaded")





# Add at each step a display size figure for the Seurat Object
print("**************************************************")
print("************LOAD GCT and CLS**********************")
print("**************************************************")

run_combat_seq(args)


print("**************************************************")
print("***********    Running CombatSeq      ************")
print("**************************************************")

# Run Combat Seq here
# CALL COMBAT SEQ ON raw_data_matrix and phenotype_vector
# adjusted <- as.data.frame(ComBat_seq(count_matrix, batch=batch, group=NULL))



print("**************************************************")
print("****************    SAVE GCT      ****************")
print("**************************************************")
