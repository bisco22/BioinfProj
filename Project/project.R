#load of the function
source("geneExpressionPerTissue.R")

#file where is stored the path to the folder containing the text files
fileInput <- paste(getwd(), "/Input/path.txt", sep = "")

#tissue selected
tissue <- "adiposeSubcut"

#threshold level for expScores
threshold <- "CIAO"

#apply the function to 
tibbles <- geneExpressionPerTissue(fileInput, tissue, threshold)