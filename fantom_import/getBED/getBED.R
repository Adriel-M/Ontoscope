# Purpose:   Access the Fantom Database and etrieve the Relevant .bed files
# Version:   1.2.0
# Date:      2016-03-26
# Author(s): Dmitry Horodetsky
#  
#
# Input:     List/Character Vector
# Output:    .bed files downloaded
#
#
# V 1.0.0    Initial Commit
#
# V 1.2.0    Fixed the 1 fantomID, 2 matches bug, new naming convention 

#Check whether the BED_DB file is present

.check_BED_DB <- function(){
  if (file.exists("BED_DB.RData")){
    load("BED_DB.RData",envir = globalenv())
    message ('BED_DB Loaded!')
  } else { stop("BED_DB not found. Please put it in your working directory")
  }
}
.check_BED_DB()

#Main Function
getBED <- function(IDs){
  for (i in IDs){
    if (substr(i,start = 1, stop = 3) == "FF:"){
      fixed_ID <- gsub("FF:","",i)
      dl_index <- grep(fixed_ID,BED_DB[,1])
      
      for (j in dl_index){
        message(paste("Downloading...",i))
        download.file(as.character(BED_DB[j,1]),paste0(j,"_",fixed_ID,".bed.gz"))
        message(paste0(j,"_",fixed_ID,".bed.gz"," Saved!"))
      }
    } else {
      stop("Ontology IDs must be in a FF:XXXXX format")
    }
  }
}




