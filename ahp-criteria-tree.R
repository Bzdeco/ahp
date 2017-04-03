# install.packages("XML")
library(XML)

setwd(dir = "/home/bzdeco/Documents/agh/badania/r/")

# parses "comparison" xml node into pairwise comparison matrix
parsePairCompMatrix <- function(comparisonsNode) {
  
  n <- length(xmlChildren(comparisonsNode)) # number of comparisons
  ratios <- as.numeric(xmlToList(comparisonsNode))
  
  # given the number of elements above the diagonal in the pairwise comparison matrix
  # we can calculate its dimensions using formula for sum of arithmetic series, where x+1
  # is seeked dimension - x^2 + x - 2n = 0
  zeros <- ceiling(Re(polyroot(c(-2*n, 1, 1)))) # need to round up (inprecise numeric result)
  dim <- zeros[zeros > 0] + 1
  
  # fill matix with pair comparison values
  pairCompMatrix <- matrix(1, nrow = dim, ncol = dim)
  r <- 1
  for(i in seq(1, dim)) {
    for(j in seq(1, dim)) {
      if(j > i) {
        pairCompMatrix[i,j] <- ratios[r]
        pairCompMatrix[j,i] <- 1 / ratios[r]
        r <- r + 1
      }
    }
  }
  
  pairCompMatrix
}

# converts pairwise comparison matrix into the weight vector using eigenvector method
getEigenWeightVector <- function(pairCompMatrix) {
  
  # the weight vector is the eigenvector corresponding the largest eigenvalue (which is reals)
  weightVector <- Re(eigen(pairCompMatrix)$vectors[,1])
  
  # normalize vector so its elements sum up to 1
  weightVector <- weightVector / sum(weightVector)
}

# converts pairwise comparison matrix into the weight vector using geometric mean method
getGeoMeanWeightVector <- function(pairCompMatrix) {
  
  # calculate geometric mean of all matrix entries in each row
  weightVector <- apply(pairCompMatrix, 1, prod) ** (1 / (dim(pairCompMatrix)[1]))
  
  # normalize vector so its elements sum up to 1
  weightVector <- weightVector / sum(weightVector)
}

# recursive function rating alternatives from the bottom of the tree
rateAlternatives <- function(rootNode) {
  
    # weight vector with rating subcriterion importance
    lvlPairCompMatrix <- parsePairCompMatrix(xmlElementsByTagName(rootNode, "comparisons")[[1]])
    lvlWeightVector <- getEigenWeightVector(lvlPairCompMatrix)
    
    subcriterions <- xmlElementsByTagName(rootNode, "criterion")
    #print(lvlWeightVector)
    
    # if there are no subcriterions it is the deepest level
    if(length(subcriterions) == 0)
      return(lvlWeightVector)
    # go deeper into the tree
    else {
      
      # list of vectors rating alternatives respectively for all subcriterions
      subVectors <- vector("list", length(subcriterions))
      i <- 1
      for(subcriterion in subcriterions) {
        # recursive call
        subVectors[[i]] <- rateAlternatives(subcriterion)
        i <- i + 1
      }
      #print(subVectors)
      
      # linear combination of subVectors with critWeights as coefficients
      weightVector <- rep(0, length(subVectors[[1]]))
      i <- 1
      for(critWeight in lvlWeightVector) {
        weightVector <- weightVector + critWeight * subVectors[[i]]
        i <- i + 1
      }
      
      return(weightVector)
    }
}

### testing ###
# print(parsePairCompMatrix(xmlChildren(rootNode)[[1]]))
# print(class(xmlElementsByTagName(rootNode, "criterion")[[1]]))
# print(getEigenWeightVector(matrix(c(1, 2, 8, 0.5, 1, 0.25, 0.125, 4, 1), nrow=3, ncol = 3, byrow = TRUE)))
# print(getGeoMeanWeightVector(matrix(c(1, 2, 8, 0.5, 1, 0.25, 0.125, 4, 1), nrow=3, ncol = 3, byrow = TRUE)))
# print(class(xmlChildren(rootNode)))

tree <- xmlTreeParse("data.xml")
rootNode <- xmlRoot(tree)
goalNode <- xmlElementsByTagName(rootNode, "goal")[[1]]

rateAlternatives(goalNode)