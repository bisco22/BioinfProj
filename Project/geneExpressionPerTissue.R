#
# This function takes as inputs:
#   -the path to a folder where are stored grouped categories of genes
#   -the name of a tissue, and 
#   -the threshold level.
# The name of the tissue is referred to the field 'name' of the 'gtexTissueV8' table
# of the 'hgFixed' database (e.g. "adiposeSubcut")
#

geneExpressionPerTissue <- function(fileInput, tissue, threshold){
  
  #write name of the script, date and time and the inputs of the function in the txt output file
  scriptLine <- paste("Script:\t\t\t\t\t\t\t\t\t\t geneExpressionPerTissue.R")
  timeLine <- paste("\nDate and Time of execution:\t\t\t\t\t\t\t", Sys.time())
  inputsSect <- paste("\nInputs\n\t-File containing the path to the categories folder:\t\t",fileInput,"\n\t-Tissue selected:\t\t\t\t\t\t\t\t",tissue,"\n\t-Threshold expScores:\t\t\t\t\t\t\t",threshold)
  write(c(scriptLine,timeLine, inputsSect), "./Output/output.txt", sep = "\n")
  
  #initialize warnings and errors lists
  warnings <- list()
  errors <- list()
  
  #creating list of packages needed
  packages <- c("DBI","RMySQL","RMariaDB","tibble","dplyr","sys","tidyr","magrittr","scriptName")
  
  #check packages existence
  for (package in packages){
    if(is.na(match(package,installed.packages()))){
      errors[['Package']] <- paste(package,"Package not installed")
      write(paste("\nError\n", errors[["Package"]], "Output/output.txt"), append = TRUE)
      stop(errors[["Package"]])
    }
  }  
  
  #load needed packages
  lapply(packages, library, character.only = TRUE)
  
  #create the output folder (checking if exists first)
  if(!dir.exists("./Output"))
    dir.create("./Output")
  
  #get the name of the script that called the function
  filename <- current_filename()
  if (is.null(filename)){
    scriptName <- "Unknown"
    warnings[["Filename"]] <- "Can't get the script filename"
  }
  else {
    filename <- strsplit(filename, "/")
    scriptName <- last(filename[[1]])
  }
  
  
  
  #check of the inputs
  if(!is.character(fileInput)){
    errors[["fileInputType"]] <- "\t-Input file must be character"
    write(paste("\nError\n",errors[["fileInputType"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["fileInputType"]])
  }
  else if (!file.exists(fileInput)){
    errors[["fileInputExists"]] <- "\t-Input file do not exists"
    write(paste("\nError\n",errors[["fileInputExists"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["fileInputExists"]])
  }
  else if(!is.character(tissue)){
    errors[["tissue"]] <- "\t-Tissue must be character"
    write(paste("\nError\n",errors[["tissue"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["tissue"]])
  }
  else if(!is.numeric(threshold)){
    errors[["thresholdType"]] <- "\t-Threshold must be numeric"
    write(paste("\nError\n",errors[["thresholdType"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["thresholdType"]])
  }
  else if(threshold < 0){
    warnings[["thresholdValue"]] <- "Threshold value should be positive"
  }
  
  #connect to UCSC hgFixed DB 
  con <- dbConnect(drv = RMariaDB::MariaDB(), username = "genomep", password = "password", db = "hgFixed", host = "genome-mysql.soe.ucsc.edu", port = 3306)
  
  #retrieve tissues query
  tissues <- dbGetQuery(con, "SELECT name FROM gtexTissueV8")
  
  #get the position of the selected tissue
  pos <- match(tissue, tissues[,1])
  
  #check if tissue exists
  if(is.na(pos)){
    errors[["position"]] <- "\t-The selected tissue could not be found"
    write(paste("\nError\n",errors[["position"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["position"]])
  }
  
  #close DB connection
  dbDisconnect(con)
  
  #read the text file with the path to the directory which contains txt file constituted by grouped categories
  path <- readLines(fileInput)
  
  #connect to UCSC hg38 DB  
  con <- dbConnect(drv = RMariaDB::MariaDB(), username = "genomep", password = "password", db = "hg38", host = "genome-mysql.soe.ucsc.edu", port = 3306)
  
  #create the list of tibbles to store results
  tibbles <- list()
  
  #adding input files caption to input section
  inputsSect <- "\t-Input files:\n"
  
  #open and process text files 
  for(file in list.files(path)){
    #add name of the input file to the output file
    inputsSect <- paste(inputsSect,"\t\t\t\t\t\t\t\t\t\t\t", file, "\n")
    
    #retrieve genes inside the category file
    #skip = 1 because examples files always start with a copyright caption on line 1
    #the header is on the 2nd line
    table <- read.table(paste(path, file, sep = "\\"), header = TRUE, skip = 1, sep="\t", quote = "\"")
    genes <- toString(sprintf("'%s'", table[,1]))
    
    #build a tibble formed by 2 columns (category and symbol)
    resultsTibble <- tibble(category = strsplit(file, "\\.")[[1]][1], symbol = table[,1])
    
    #get the expScores for the genes that are both in the file and in the DB
    results <- dbGetQuery(con, paste("SELECT name, CONVERT(expScores USING utf8) AS expScores FROM gtexGeneV8 WHERE name in (",genes,")") )
    
    #process expScores to get a list of separated values
    results$expScores <- strsplit(results$expScores, ",")
    
    #retrieve only the value for expScores of the selected position
    results$expScores <- as.numeric(sapply(results$expScores,"[[",pos))
    
    #add a column to the tibble indicating if the expScore for the selected tissue is higher than the threshold
    results <- mutate(results, isOver = results$expScores > threshold)
    
    #build a tibble from the results of the query
    resultsTibble <- left_join(resultsTibble, select(results, name, isOver), by = c("symbol" = "name"))
    
    #generate a wide format tibble 
    wide <- resultsTibble %>% pivot_wider(names_from = symbol, values_from = isOver)
    
    #save the results in the list of tibbles
    tibbles[[strsplit(file, "\\.")[[1]][1]]] <- wide
  }
  
  #initialize the outputs section
  outputsSect <- "Outputs\n\t-tibbles generated from the following grouped categories:\n"
  
  #add tibbles to outputs section of the txt file
  outputsText <- paste("\t\t\t\t\t\t\t\t\t\t\t", names(tibbles), "\n", collapse = " ")
  outputsSect <- paste(outputsSect, outputsText)
  
  #close connection
  dbDisconnect(con)
  
  #add warnings to warningsSect
  if (length(warnings) != 0){
    warningsSect <- "Warnings\n"
    warningsText <- paste("\t-", warnings, "\n", collapse = "", sep = "")
    warningsSect <- paste(warningsSect, warningsText)
  } else 
    warningsSect <- ""

  #write the inputs, ouputs and warnings sections in the output text file
  fileOutput <- c(inputsSect, outputsSect, warningsSect)
  write(fileOutput, "./Output/output.txt", append = TRUE, sep = "\n")
  
  #return the list of tibbles generated
  return(tibbles)
}

