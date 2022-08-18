#
# This function takes as inputs:
#   -a path to a folder where are stored grouped categories of genes
#   -a name of a tissue, and 
#   -a threshold level.
# The name of the tissue is referred to the field 'name' of the 'gtexTissueV8' table
# of the 'hgFixed' database (e.g. "adiposeSubcut")
#

geneExpressionPerTissue <- function(fileInput, tissue, threshold){

  #loading packages
  lapply(c("DBI","RMySQL","RMariaDB","tibble","dplyr","sys","tidyr","magrittr","scriptName"), library, character.only = TRUE)
  
  #creating the output folder (checking if exists first)
  if(!dir.exists("./Output"))
    dir.create("./Output")
  
  #initializing errors list
  errors <- list()
  
  #getting the name of the script that called the function
  filename <- current_filename()
  if (is.null(filename)){
    scriptName <- "Unknown"
    errors[["Filename"]] <- "Can't get the script filename"
  }
  else {
    filename <- strsplit(filename, "/")
    scriptName <- last(filename[[1]])
  }
  
  #writing name of the script, date and time in the txt output file
  scriptLine <- paste("Script:\t\t\t\t\t\t\t\t", scriptName)
  timeLine <- paste("\nDate and Time of execution:\t\t\t\t\t", Sys.time())
  write(c(scriptLine,timeLine), "./Output/output.txt", sep = "\n")
  
  #initializing the strings of input, output and error messages
  inputs <- paste("\nInputs\n\t-File containing categories folder path:\t\t",fileInput,"\n\t-Tissue selected:\t\t\t\t\t\t",tissue,"\n\t-Threshold expScores:\t\t\t\t\t",threshold, "\n\t-Input files:\n")
  outputs <- "Outputs\n\t-tibbles generated from:\n"
  errorLine <-  "\nErrors\n"

  
  #check of the inputs
  if(!is.character(fileInput)){
    errors[["fileInputType"]] <- "fileInput must be character"
    write(paste(errorLine,errors[["fileInputType"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["fileInputType"]])
  }
  else if (!file.exists(fileInput)){
    errors[["fileInputExists"]] <- "fileInput do not exists"
    write(paste(errorLine,errors[["fileInputExists"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["fileInputExists"]])
  }
  else if(!is.character(tissue)){
    errors[["tissue"]] <- "tissue must be character"
    write(paste(errorLine,errors[["tissue"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["tissue"]])
  }
  else if(!is.numeric(threshold)){
    errors[["threshold"]] <- "threshold must be numeric"
    write(paste(errorLine,errors[["threshold"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["threshold"]])
  }
  
  #connecting to UCSC hgFixed DB 
  con <- dbConnect(drv = RMariaDB::MariaDB(), username = "genomep", password = "password", db = "hgFixed", host = "genome-mysql.soe.ucsc.edu", port = 3306)
  
  #retrieving tissues query
  tissues <- dbGetQuery(con, "SELECT name FROM gtexTissueV8")
  
  #getting the position of the selected tissue
  pos <- match(tissue, tissues[,1])
  
  #checking if tissue exists
  if(is.na(pos)){
    errors[["position"]] <- "The selected tissue could not be found"
    write(paste(errorLine,errors[["position"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["position"]])
  }
  
  #closing connection
  dbDisconnect(con)
  
  #reading the text file with the path to the directory which contains txt file constituted by grouped categories
  path <- readLines(fileInput)
  
  #connecting to UCSC hg38 DB  
  con <- dbConnect(drv = RMariaDB::MariaDB(), username = "genomep", password = "password", db = "hg38", host = "genome-mysql.soe.ucsc.edu", port = 3306)
  
  #creating the list of tibbles used to store results
  tibbles <- lst()
  
  #opening and processing text files 
  for(file in list.files(path)){
    inputs <- paste(inputs,"\t\t\t\t\t\t\t\t\t", file, "\n")
    
    #retrieving categories inside the file
    table <- read.table(paste(path, file, sep = "\\"), header = TRUE, skip = 1, sep="\t", quote = "\"")
    genes <- toString(sprintf("'%s'", table[,1]))
    resultsTibble <- tibble(category = strsplit(file, "\\.")[[1]][1], symbol = table[,1])
    
    #retrieving genes and expScores in relation to the selected categories
    query <- dbSendQuery(con, paste("SELECT name, CONVERT(expScores USING utf8) AS expScores FROM gtexGeneV8 WHERE name in (",genes,")") )
    results <- dbFetch(query)
    dbClearResult(query)
    
    #processing expScores to get a list of separated values
    results$expScores <- strsplit(results$expScores, ",")
    
    #retrieving only the value for expScores of the selected position
    results$expScores <- as.numeric(sapply(results$expScores,"[[",pos))
    
    #adding a column to the tibble indicating if the expScore for the selected tissue is higher than the threshold
    results <- mutate(results, isOver = results$expScores > threshold)
    
    #creating a tibble from the results of the query
    resultsTibble <- left_join(resultsTibble, select(results, name, isOver), by = c("symbol" = "name"))
    
    #generating a wide format tibble 
    wide <- resultsTibble %>% pivot_wider(names_from = symbol, values_from = isOver)
    
    #saving the results in the list of tibbles
    tibbles[[strsplit(file, "\\.")[[1]][1]]] <- wide
  }
  
  #adding tibbles to outputs section of the txt file
  outputTibbles <- paste("\t\t\t\t\t\t\t\t\t", names(tibbles), "\n", collapse = " ")
  outputs <- paste(outputs, outputTibbles)
  
  #closing connection
  dbDisconnect(con)
  
  #adding errors to errorLine
  errorLine <- paste(errorLine, "\t-", errors, "\n")
  
  #writing the input, ouputs and errors in the output text file
  fileOutput <- c(inputs, outputs, errorLine)
  write(fileOutput, "./Output/output.txt", append = TRUE, sep = "\n")
  
  #return the list of tibbles generated
  return(tibbles)
}

