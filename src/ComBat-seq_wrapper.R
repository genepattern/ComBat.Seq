########################################################
# Parse Parameters
########################################################
print('================================================')
print("Loading library: optparse")
library("optparse")


# Parse input arguments
parser = OptionParser()
# parameter types: 'character', 'integer', 'logical', 'double', or 'complex'
# =================================================
# Input file
parser <- add_option(parser, c("--input_matrix"), help="GCT file to load
                     ( containing the raw counts).", default = "NO_FILE")
parser <- add_option(parser, c("--input_class"), help="CLS file to load
                               (containing the same phenotypes).", default = "NO_FILE")
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
if (grepl(".gct", args$input_matrix, fixed = TRUE)){
  USE_GCT <<- TRUE
}

print(paste("GCT FILE IS PROVIDED?", USE_GCT))

# Setting up the PDF file for the plots
# pdf(file=paste(args$file_name, '.pdf', sep=''))

##########################################################
# Function Definitions
##########################################################

## Load the GCT file
read_gct <- function(input_file){
  write("Reading GCT file...", stdout())
  # Read GCT file as a matrix

  a_df <- read.table(input_file, sep = "\t", skip = 2) # read the gct file as a dataframe
  my_colnames <- a_df[1,]                              # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                        # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                    # dropping the first row now that it's the column names
  a_df <- subset(a_df, select=-c(Description))         # Dropping the "Description" column
  a_df <- subset(a_df, select=-c(Name))                # Dropping the "Name" column
  a_df <- sapply(a_df,as.numeric)                      # turn the dataframe into a matrix
  names(a_df) <- NULL
  return(as.matrix(a_df))
}


# returns the number of columns given a file
cols <- function(input_file){
  a_df <- read.table(input_file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]
  return (ncol(a_df)-2)
}

# returns the number of rows given a file
rows <- function(input_file){
  a_df <- read.table(input_file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]
  return (nrow(a_df))
}


# need to return the name column bc when turning the gct to df to matrix, from df to matrix, the index column disappears
return_name <- function(input_file){
  a_df <- read.table(input_file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                     # dropping the first row now that it's the column names
  name <- a_df$Name
  return (name)
}


# returns the description column
return_description <- function(input_file){
  a_df <- read.table(input_file, sep = "\t", skip = 2)  # read the gct file as a dataframe
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                     # dropping the first row now that it's the column names
  description <- a_df$Description
  return (description)
}

# read the cls file and output an array for the batch labels
read_cls <- function(input_file){
  write("Reading CLS file...", stdout())
  cls_list = scan(file=input_file, sep=" ", skip = 2)  # skip first two lines, capture batch numbers as a list
  cls_vector = array(as.numeric(unlist(cls_list)))                     # turn list into a vector (1-D array)
  return(cls_vector)
}


# save the Combat-seq adjusted data as a GCT file
save_data <- function(fileName, batch_corrected_matrix, input_file){

  a_df <- as.data.frame(batch_corrected_matrix)
  a_df$Description <- return_description(input_file)                              # insert description column
  a_df$Name <- return_name(input_file)                                            # insert name column
  a_df <- a_df[c(cols(input_file)+2, cols(input_file)+1, 1:cols(input_file))]     # move name and description column to the front
  col_names <- names(a_df)
  a_df <- rbind(col_names, a_df)                                                  # add the columns names as the first row in df
  colnames(a_df) <- NULL

  f <- file(fileName)
  first_second = paste("#1.2", paste(as.character(rows(input_file)), as.character(cols(input_file)), sep="\t"), sep="\n")
  write(first_second, file=fileName)

  write.table(a_df, quote=FALSE, row.names = FALSE, file=fileName, sep="\t", append=TRUE) # write data from a_df to the file
  close(f)
}

# Run the wrapper
run_combat_seq <- function(args){
  ### Run ComBat_Seq based on what the user inputted.
  if (USE_GCT){
    raw_data_matrix <- read_gct(args$input_matrix)
    phenotype_vector <- read_cls(args$input_class)
    adjusted <- as.data.frame(ComBat_seq(raw_data_matrix, batch=phenotype_vector))
    print('... done')
    print(adjusted)
    save_data(paste(args$file_name, '.gct', sep=''), batch_corrected_matrix=adjusted, input_file = args$input_matrix)
  }else{
    raw_data_matrix <- as.matrix(sapply(read.table(args$input_matrix),as.numeric))
    phenotype_vector <- array(as.numeric(unlist(scan(file=args$input_class))))
    adjusted <- as.data.frame(ComBat_seq(raw_data_matrix, batch=phenotype_vector))
    print('... done')
    print(adjusted)
    write.table(adjusted, file = paste(args$file_name, "_adjusted.csv", sep=","))
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
