## Load the GCT file
read_gct <- function(input_file){
  write("Reading GCT file...", stdout())
  # Read GCT file as a matrix

  a_df <- read.table(input_file, sep = "\t", skip = 2) # read the gct file as a dataframe
  genes <- return_name(a_df)
  descriptions <- return_description(a_df)
  samples <- a_df[1,]
  my_colnames <- a_df[1,]                              # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                        # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                    # dropping the first row now that it's the column names
  a_df <- subset(a_df, select=-c(Description))         # Dropping the "Description" column
  a_df <- subset(a_df, select=-c(Name))                # Dropping the "Name" column
  a_df <- sapply(a_df,as.numeric)                      # turn the dataframe into a matrix


  names(a_df) <- NULL
  returnlist <- list("data" = as.matrix(a_df), "genes" = genes,
          "samples" = my_colnames[3:length(my_colnames)], "descriptions" = descriptions)
  return(returnlist)
}

read_delimited_file <- function(input_file){
  #' This function takes in the input file path and returns a matrix of the count values
  #' input_file --> filepath to the input counts matrix

  write("Reading Tab Delimited file ... ", stdout())
  df <- read.table(input_file, sep = '\t')
  raw_counts <- all(sapply(df, is.numeric))

  if(raw_counts){
    ## if only raw counts
    num_rows <- nrow(df)
    num_cols <- ncol(df)
    df <- sapply(df, as.numeric)

    returnlist <- list("data" = as.matrix(df),
        "num_rows" = num_rows,
        "num_cols" = num_cols,
        "sample_names"=NULL,
        "gene_names"=NULL,
        "raw_counts"= TRUE)
    return(returnlist)
  }else{
    # Reading the tab delimited file as a matrix
    # if reading in delimited file with samples and gene names
    num_rows <- nrow(df)
    num_cols <- ncol(df)
    samples <- df[1,]
    genes <- df[,1]

    #resize df, so only the raw data is used
    df <- df[-1,]
    df <- df[, -1]
    df <- sapply(df,as.numeric)
    returnlist <- list("data" = as.matrix(df),
        "num_rows" = num_rows,
        "num_cols" = num_cols,
        "sample_names" = samples,
        "gene_names" = genes,
        "raw_counts"= FALSE)
    return(returnlist)
  }

}

read_multiline_batch_info_file <- function(data_matrix, batch_file, use_cls){
  #' This function reads in the batch information, parses it, and
  #' returns batch, groups, and any other information

  print("Reading batch file...")

  if (use_cls){
    df <- read.table(batch_file, skip = 2)
    print(df)
  }else{
    df <- read.table(batch_file, sep = '\t')
    print(df)
  }


  if (nrow(df) > 1){
    if ((ncol(df)) == ncol(data_matrix)){
      for (row_name in row.names(df)){
        if ((tolower(row_name)) == 'batch'){
          ## make a list called batch, these are the batches
          batches <- df[row_name,]
        } else if (tolower(row_name) == 'group'){
          groups <- df[row_name,]
        }
      }
    }else{
      print(paste("data columns and batches shape mismatch! columns of count data and number of batches should match. # batches: ", (ncol(df) - 1), " # data columns: ", (ncol(data_matrix))))
    }
  }else{
    batches <- df[1,]
    groups <- NULL
  }

  returnlist <- list("batches" = batches, "groups" = groups)
  return(returnlist)

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
return_name <- function(a_df){
  my_colnames <- a_df[1,]                               # saving the data in the first row to variable my_colnames
  colnames(a_df) <- my_colnames                         # setting my_colnames to be the column name in the df
  a_df <- a_df[-1,]                                     # dropping the first row now that it's the column names
  name <- a_df$Name
  return (name)
}


# returns the description column
return_description <- function(a_df){
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

# printing dimensions
print_dimensions <- function(arr, arr_name){
  ## Prints the dimensions of *arr*, with *arr_name*
  ## returns the dimensions as well
  print_str = paste(
    "Dimension of ", arr_name, " is: ",
    nrow(arr), " rows, and ",
    ncol(arr), " columns.",
    sep = ""

  )
  print(print_str)
  return(list(
    "rows" = nrow(arr),
    "cols" = ncol(arr)
  ))


}


# save the Combat-seq adjusted data as a GCT file
save_data <- function(fileName, batch_corrected_matrix,
    samples = NULL, genes = NULL, descriptions = NULL, use_gct){

  ## First, print dimensions to make sure everything lines up
  print_dimensions(batch_corrected_matrix, "corrected matrix")
  print_dimensions(samples, "samples vector")
  print_dimensions(genes, "genes vector")
  print_dimensions(descriptions, "descriptions vector")


  if (use_gct && !is.null(samples) && !is.null(genes)){
    ## append columns, genes, descriptions to build GCT file

    colnames(batch_corrected_matrix) <- samples
    batch_corrected_matrix$Description <- descriptions
    batch_corrected_matrix$Name <- genes


    f <- file(fileName)
    ## create first two rows for GCT file
    row1 = paste("#1.2", paste(as.character(nrow(batch_corrected_matrix)), as.character(ncol(batch_corrected_matrix)), sep="\t"), sep="\n")
    write(row1, file=fileName)
    write.table(batch_corrected_matrix, quote=FALSE, row.names = FALSE, file=fileName, sep="\t", append=TRUE) # write data from a_df to the file
    close(f)
  }else if (!use_gct){
    if (!is.null(samples) && !is.null(genes)){
      colnames(batch_corrected_matrix) = samples
      rownames(batch_corrected_matrix) = genes
      write.table(batch_corrected_matrix, file = fileName, sep = '\t')
      print(paste('Table written as: ', fileName, sep = ''))
    }else{
      ## raw counts were provided, can just return
      write.table(batch_corrected_matrix, file = fileName, sep = '\t', row.names = FALSE, col.names = FALSE)
      print(paste('Table written as: ', fileName, sep = ''))
    }

  }


}
