##
# DNAphi analysis
##
library(DNAshapeR)
library(caret)
library(ggplot2)
library(grid)

## Set parameters
groupName <- "test"
folderPath <- "C:\\Users\\tsupeich\\Documents\\GitHub\\DNAphi_analysis\\example_datasets\\"
faFileNames <- list.files(folderPath, pattern = "txt.fa$")
outputFileName <- paste0(format(Sys.time(), format = "%Y-%m-%d-%H%M%S"), ".txt")

featureNames <- list(
  c("1-MGW"),
  c("1-EP"),
  c("1-mer"),
  c("1-mer", "1-MGW"),
  c("1-mer", "1-EP"),
  c("1-mer", "1-shape"),
  c("1-mer", "1-ProT", "1-Roll", "1-HelT", "1-EP"),
  c("1-mer", "1-shape", "1-EP")
)

## Run prediction
i <- 0
for(faFileName in faFileNames){
  i <- i + 1
  print ( paste("Running progress:", i, "/", length(faFileNames)) )
  for(featureName in featureNames){
    fn <- paste0(folderPath,"\\", faFileName)
    print(featureName)
    print(faFileName)
    
    print("shape prediction...")
    pred <- getShape(fn)
    
    print("encoding...")
    featureVector <- encodeSeqShape(fn, pred, featureName, normalize = TRUE)
    
    print("machine learning...")
    experimentalData <- read.table( sub("\\.fa", replacement = "", fn) )
    df <- data.frame(affinity=experimentalData$V2, featureVector)
    
    set.seed(168)
    trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
    model <- train(affinity~ ., data = df, trControl=trainControl, method="glmnet",
                   tuneGrid=data.frame(alpha = 0, lambda=c(2^c(-15:15))))
    
    rs <- model$results$Rsquared[1]
    rsd <- model$results$RsquaredSD[1]
    
    count <- nrow(featureVector)
    result <- paste(faFileName, groupName, count, paste(featureName, collapse="|"), rs, rsd)
    
    write(result, file = paste0(folderPath,"\\", outputFileName), append = TRUE)
    
  }
}

# Write output into file
# Output format:
# fileName | groupName | count | MGW | EP | seq | seq+MGW | seq+EP | seq+shape | seq+3shape+ep | seq+shape+ep
fileName <- paste0(folderPath, "\\", outputFileName)
qt <- read.table(fileName)

for( i in 1:floor(nrow(qt)/8) ){
  result <- unname( cbind(qt[8*i-7, 1:3], t(qt[(8*i-7):(8*i), 5])) )
  write.table(result, file = paste0(fileName, ".formated.txt"), append = TRUE, quote = FALSE)
}

