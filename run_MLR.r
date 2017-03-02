##
# DNAphi analysis
##
library(DNAphiR)
library(caret)
library(ggplot2)
library(grid)

## Set parameters
groupName <- "group"
folderPath <- paste0( "file_path", groupName )
faFileNames <- list.files(folderPath, pattern = "txt.fa$")
flagPermute <- 0

featureNames <- list(
    c("1-MGW"),
    c("1-EP"),
    c("1-mer"),
    c("1-mer", "1-MGW"),
    c("1-mer", "1-EP"),
    c("1-mer", "1-shape"),
    c("1-mer", "1-ProT", "1-Roll", "1-HelT", "1-EP"),
    c("1-mer", "1-shape", "1-EP"),
    c("1-mer", "1-EP", "1-MGW")
)

outputFileName <- paste0(format(Sys.time(), format = "%Y-%m-%d-%H%M%S"), ".txt")
weightFileName <- paste0(folderPath, "\\", outputFileName, ".weighted.txt")


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
    write.table(featureVector, file = paste0(folderPath,"\\", outputFileName, ".", faFileName, ".", paste(featureName, collapse="_"), ".feature.txt"))

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

    # Write weight file
    weights <- coef(model$finalModel, model$bestTune$lambda)
    weights <- as.matrix(weights)
    write(faFileName, file = weightFileName, append = TRUE)
    write.table(weights, file = weightFileName, append = TRUE)
  }
}


## Format output file
# fileName | groupName | count | MGW | EP | seq | seq+MGW | seq+EP | seq+shape | seq+3shape+ep | seq+shape+ep | 1-mer+1-EP+1-MGW

formatedFileName <- paste0(folderPath, "\\", outputFileName, ".formated.txt")
qt <- read.table( paste0(folderPath,"\\", outputFileName) )

sizeFeature <- length(featureNames)

#####
