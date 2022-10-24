#load of the function
source("geneExpressionPerTissue.R")

#file where is stored the path to the folder containing the text files
fileInput <- paste(dirname(filename), "/Input/path.txt", sep = "")

#tissue selected
tissue <- "adiposeSubcut"

#threshold level for expScores
threshold <- 10

#apply the function to 
tibbles <- geneExpressionPerTissue(fileInput, tissue, threshold)
