geneExpressionPerTissue <- function(fileInput, tissue, threshold){
  
  #loading packages
  lapply(c("DBI","RMySQL","RMariaDB","tibble","dplyr","sys","magrittr","scriptName"), library, character.only = TRUE)
  
  #getting the name of the script that called the function
  filename <- current_filename()
  filename <- strsplit(filename, "/")
  scriptName <- last(filename[[1]])
  
  #writing name of the script, date and time in the txt output file
  scriptLine <- paste("Script:\t\t\t\t\t", scriptName)
  timeLine <- paste("\nDate and Time:\t\t\t\t", Sys.time())
  write(c(scriptLine,timeLine), "./Output/output.txt", sep = "\n")
  
  #initializing the strings of input, output and error messages
  inputs <- paste("\nInputs\n\t-Categories folder's path:\t",fileInput,"\n\t-Tissue selected:\t\t\t",tissue,"\n\t-Threshold expScores:\t\t",threshold, "\n\t-Input files:\n")
  outputs <- "Outputs\n\t-tibble for:\n"
  errorLine <-  "\nErrors\n\t"
  errors <- list()
  
  #check of the inputs
  if(!is.character(fileInput)){
    errors[["fileInput"]] <- "fileInput must be character"
    write(paste(errorLine,errors[["fileInput"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["fileInput"]])
  }else if(!is.character(tissue)){
    errors[["tissue"]] <- "tissue must be character"
    write(paste(errorLine,errors[["tissue"]]), "./Output/output.txt", append = TRUE)
    stop(errors[["tissue"]])
  }else if(!is.numeric(threshold)){
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
  if(pos < 1 || pos > 54 || is.na(pos)){
    cat(pos)
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
    inputs <- paste(inputs,"\t\t\t\t\t\t", file, "\n")
    
    #retrieving categories inside the file
    categories <- 0
    categories <- read.table(paste(getwd(), "/Input/Categories/", file, sep=""))
    categories <- toString(sprintf("'%s'", categories))
    
    #retrieving genes and expScores in relation to the selected categories
    query <- dbSendQuery(con, paste("SELECT geneType, name, CONVERT(expScores USING utf8) AS expScores FROM gtexGeneV8 WHERE geneType in (",categories,")") )
    results <- dbFetch(query)
    dbClearResult(query)
    
    #creating a tibble from the 
    resultsTibble <- as_tibble(results)
    
    #processing expScores to get a list of separated values
    resultsTibble$expScores <- strsplit(resultsTibble$expScores, ",")
    resultsTibble$expScores <- as.numeric(sapply(resultsTibble$expScores,"[[",pos))
    
    #adding a column to the tibble indicating if the expScore for the selected tissue is higher than the threshold
    resultsTibble <- mutate(resultsTibble, isOver = resultsTibble$expScores > threshold)
    
    #group the lines of the tibble by category and rearrange the format of the tibble to wide saving only one gene (the first) per category
    group <- resultsTibble %>%
      group_by(geneType) %>% 
      summarise(geneType = first(geneType), name = first(name), isOver = first(isOver))
    
    #saving the results of the process in the list of tibbles
    tibbles[[file]] <- as_tibble(group)
    outputs <- paste(outputs,"\t\t\t\t\t\t", file, "\n")
  }
  
  #closing connection
  dbDisconnect(con)
  
  #writing the input, ouputs and errors in the output text file
  fileOutput <- c(inputs, outputs, errorLine)
  write(fileOutput, "./Output/output.txt", append = TRUE, sep = "\n")
  
  #return the list of tibbles generated
  return(tibbles)
}

