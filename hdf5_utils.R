library("rhdf5")

load.from.taiga <- function(dataset_id=NULL, name=NULL, version=NULL) {
  if(!is.null(dataset_id)) {
    source <- paste("http://datasci-dev:8999/rest/v0/datasets/",dataset_id,"?format=rdata", sep='');
  } else {
    stopifnot(!is.null(name))
    source <- paste("http://datasci-dev:8999/rest/v0/namedDataset?fetch=content&format=rdata&name=",name,sep='');
    if(!is.null(version)) {
      source <- paste(source, "&version=", version, sep='')
    }
  }
  
  cat("Loading", source, "from Taiga\n")
  
  e <- new.env()
  load(url(source), envir=e);
  
  e$data
}

##
# Loads .Rdata from file. The R object must have been saved as 'data'
# @name: file path
#
load.from.file <- function(name){
  cat("Loading", name, "from file\n")
  
  e <- new.env()
  load(name, envir=e);
  
  e$data
}


##
# Creates a group for dataset.vals where dataset is the dataset name. Each column of the matrix is stored as dataset.val/colname,
# creates a group for dataset.headers and stores the rownames, adds the dataset label at labels/dataset
# @name: dataset filename
# @label: human readable label for dataset
# @data_matrix: matrix or data.frame to insert column-wise to hdf5
# @hdf5_file: file path
#
convert_to_hdf5 <- function(name,label,data_matrix,version,out_dir){
  
  #Create emtpy hdf5 file
  hdf5_file <- file.path(out_dir,paste0(version,".h5"))
  h5createFile(hdf5_file)
  
  #Write each column to an index in the values group
  h5createGroup(hdf5_file,"values")
  for (g in colnames(data_matrix)){
    h5write(data_matrix[,g],file = hdf5_file,paste("values",g,sep="/"))
  }
  
  #Write the column headers, row headers, and dataset label to the annotaitons group
  h5createGroup(hdf5_file, "annotations")
  h5write(rownames(data_matrix),file = hdf5_file,paste("annotations","rownames",sep="/"),write.attributes=F)
  h5write(colnames(data_matrix),file = hdf5_file,paste("annotations","colnames",sep="/"),write.attributes=F)
  h5write(label,file = hdf5_file,paste("annotations","label",sep="/"),write.attributes=F)
  
  H5close()
}

##
# Get a column vector from the dataset. 
# For use with HDF5 files created using convert_to_hdf5 function
# @dataset_id: unique identifier for dataset version, for datasets on taiga it is the taiga id
# @column_name: name of the column from original matrix or df that will be returned from hdf5
# @hdf5_dir: path to directory with hdf5 files
# @return: named vector where values are the column's values from given dataset and names are rownames
#
get_values <- function(dataset_id,column_name,hdf5_dir){
  
  hdf5_file <- paste(hdf5_dir,dataset_id,".h5",sep="")
  vec <- h5read(hdf5_file,paste("values",column_name,sep="/"),read.attributes=F)
  names(vec) <- h5read(hdf5_file,paste("annotations","rownames",sep="/"),read.attributes=F)
  H5close()
  return(vec)
}

##
# Gets human readable label for the dataset contained in the hdf5 file for use in labelling the axes
# For use with HDF5 files created using convert_to_hdf5 function
# @dataset_id: unique identifier for dataset version, for datasets on taiga it is the taiga id
# @hdf5_dir: path to directory with hdf5 files
# @return: string 
#
get_label <- function(dataset_id,hdf5_dir){
  hdf5_file <- paste(hdf5_dir,dataset_id,".h5",sep="")
  label <- h5read(hdf5_file,paste("annotations","label",sep="/"),read.attributes=F)
  H5close()
  return(label)
}

##
# Gets column names for the specified dataset
# For use with HDF5 files created using convert_to_hdf5 function
# @dataset_id: unique identifier for dataset version, for datasets on taiga it is the taiga id
# @hdf5_dir: path to directory with hdf5 files
# @return: character vector of column names 
#
get_colnames <- function(dataset_id,hdf5_dir){
  hdf5_file <- paste(hdf5_dir,dataset_id,".h5",sep="")
  cnames <- h5read(hdf5_file,paste("annotations","colnames",sep="/"),read.attributes=F)
  H5close()
  return(cnames)
}

##
# Adds taiga datasets to hdf5
# Currently adds Demeter and lineage info from file since it isn't on tiaga yet
# @out_file: ABSOLUTE file path of hdf5 file to create
# @taiga_config_file: tab delimited text file with column headers -- "taiga_id", "name", "label", "transpose"
#
compile_hdf5s_from_taiga <- function(out_dir,taiga_config_file){
  
  taiga_config <- read.delim(taiga_config_file, stringsAsFactors=FALSE)
  
  for (i in 1:nrow(taiga_config)){
    
    if (length(list.files(out_dir,pattern=paste(taiga_config[i,"taiga_id"],".h5",sep=""))) == 0){
      
      data <- load.from.taiga(dataset_id=taiga_config[i,"taiga_id"])
      
      #Must transpose if necessary before doing dataset specific modifications
      if (taiga_config[i,"transpose"]){
        data <- t(data)
      }
      
      #Dataset specific modifications
      if (taiga_config[i,"name"] == "demeter"){
        data <- data[,colSums(is.na(data)) != nrow(data)]
      }
      if (taiga_config[i,"name"] == "gene.rpkm"){
        colnames(data)<-gsub(" .+$", "", colnames(data))
      }
      
      convert_to_hdf5(taiga_config[i,"name"],taiga_config[i,"label"],data,taiga_config[i,"taiga_id"],out_dir)
      
    }
    
  }
  
}