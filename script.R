### ======================================================================== ###
### === A Monte Carlo simulation study of the Gaussian graphical model ===== ###
### === using graphical lasso based on the extended Bayes information ====== ###
### === criterion. ========================================================= ### 
### ======================================================================== ###


# 1. Necessary libraries ------------------------------------------------------- 

library(BDgraph)
library(corpcor)
library(qgraph)
library(bootnet)
library(tidyverse)
library(ggplot2)


# 2. Simulation conditions -----------------------------------------------------

# network parameters
numberOfNodes <- c(5, 10, 15, 20, 25)
ratioOfEdgesPresence <- c(.40, .60, .80)
ratioOfPositiveEdges <- 1.0

# data parameters
skewness <- c("symmetric", "skewed", "highly skewed")
numberOfObservations <- c(100, 200, 500, 1000, 2000)

# number of samples for each condition
numberOfSamples <- 1000


# 3. Function creating a partial correlation matrix which has: -----------------
# -- a specific number of variables and number of non-zero coefficients --------

generatePcorMatrix <- function(nodes, edgeRatio)
{
  numberOfEdges <- edgeRatio * (nodes * (nodes - 1) / 2)
  adjacencyMatrix <- bdgraph.sim(p = nodes, size = numberOfEdges,
                                 graph = "cluster", class = 1)$G
  constant <- 1.5
  parRange <- c(0.5, 1)
  precMatrix <- adjacencyMatrix
  class(precMatrix) <- "matrix"
  
  # Make edges negative and add weights:
  precMatrix[upper.tri(precMatrix)] <-
    precMatrix[upper.tri(precMatrix)] *
    sample(c(-1, 1), sum(upper.tri(precMatrix)), TRUE,
           prob = c(ratioOfPositiveEdges, 1-ratioOfPositiveEdges)) *
    runif(sum(upper.tri(precMatrix)), min(parRange), max(parRange))
  
  # Symmetrize:
  precMatrix[lower.tri(precMatrix)] <- t(precMatrix)[lower.tri(precMatrix)]  
  
  # Make positive definite:
  diag(precMatrix) <- constant * rowSums(abs(precMatrix))
  diag(precMatrix) <- ifelse(diag(precMatrix)==0,1,diag(precMatrix))
  precMatrix <- precMatrix/diag(precMatrix)[row(precMatrix)]
  precMatrix <- (precMatrix + t(precMatrix)) / 2
  
  pcorMatrix <- as.matrix(wi2net(precMatrix))
  
  return(pcorMatrix)
}


# 4. Creation of populational networks for all conditions ----------------------

set.seed(1)

# this list will contain matrices of populational edge weights for all conditions
# the coordinates are: [[numberOfNodes]][[ratioOfEdgesPresence]]
matrixPopulation <- vector(mode = "list", length = length(numberOfNodes))

for (i in 1:length(numberOfNodes))
{
  matrixPopulation[[i]] <- vector(mode = "list", length = length(ratioOfEdgesPresence))
  
  for (j in 1:length(ratioOfEdgesPresence))
  {
    matrixPopulation[[i]][[j]] <- generatePcorMatrix(numberOfNodes[i],
                                           ratioOfEdgesPresence[j])
  }
}


# 5. Functions calculating mean and SD of edges --------------------------------

calculateMean <- function(matrix)
{
  return(round(mean(matrix[lower.tri(matrix)]), digits = 3))
}

calculateSD <- function(matrix)
{
  return(round(sd(matrix[lower.tri(matrix)]), digits = 3))
}

# 6. Visualization of populational networks ------------------------------------

# in order to visualize all networks using the same scale for their edges
# find the maximum value of edge weights in all of the populational matrices
maxEdge <- 0

for (i in 1:length(numberOfNodes))
{
  for (j in 1:length(ratioOfEdgesPresence))
  {
    if (max(matrixPopulation[[i]][[j]]) > maxEdge)
    {
      maxEdge <- max(matrixPopulation[[i]][[j]])
    }
  }
}

# plot all populational networks on a figure
ratioDescription <- "Ratio = "
meanDescription <- "; Mean = "
sdDescription <- "; SD = "
maxDescription <- "; Max = "

jpeg("figures/Fig1.jpeg", units = "in", width=9, height=15, quality = 100, res=300)
par(mfrow = c(length(numberOfNodes), length(ratioOfEdgesPresence)))

for (i in 1:length(numberOfNodes))
{
  for (j in 1:length(ratioOfEdgesPresence))
  {
    # creation of reader-friendly description of percentage of present edges
    ratio <- paste0(ratioOfEdgesPresence[j] * 100, "%")
    
    # calculation of mean, sd, and max edge with 3 decimal places accuracy
    edgeMean <- format(calculateMean(matrixPopulation[[i]][[j]]), nsmall = 3)
    edgeSD <- format(calculateSD(matrixPopulation[[i]][[j]]), nsmall = 3)
    edgeMax <- format(round(max(matrixPopulation[[i]][[j]]), digits = 3), nsmall = 3)
    
    qgraph(input = matrixPopulation[[i]][[j]], layout = "spring",
           maximum = maxEdge, theme = "colorblind",
           vsize = 12, border.width = 2,
           title = paste0(ratioDescription, ratio,
                          meanDescription, edgeMean,
                          sdDescription, edgeSD,
                          maxDescription, edgeMax),
           title.cex = 1)
  }
}
dev.off()
  
  
# 7. Function generating one sample of data ------------------------------------

# this is a simplified version of the function obtained from running
# the ggmGenerator(ordinal = TRUE, nLevels = 5,
#                  skewFactor = 1, type = "uniform", missing = 0)
# the modifications included:
# - removal of unnecessary parts of code (unused conditions)
# - defining the particular threshold values
generateData <- function(input, skewness = "symmetric", n = 100)
{
  Sigma <- cov2cor(solve(diag(ncol(input)) - input))
  Data <- mvtnorm::rmvnorm(n, sigma = Sigma)
  
  if (skewness == "symmetric")
  {
    thresholds <- c(qnorm(0.1), qnorm(0.3), qnorm(0.7), qnorm(0.9))
  }
  else if (skewness == "skewed")
  {
    thresholds <- c(qnorm(0.3), qnorm(0.7), qnorm(0.9), qnorm(0.975))
  }
  else if (skewness == "highly skewed")
  {
    thresholds <- c(qnorm(0.50), qnorm(0.75), qnorm(0.90), qnorm(0.975))
  }
  
  
  for (i in 1:ncol(Data))
  {
    Data[,i] <- as.numeric(cut(Data[,i],sort(c(-Inf,thresholds,Inf))))
  }
  
  return(Data)
}

# 8. Function generating all sample of data for a single condition -------------

identifyMatrix <- function(matrix)
{
  for (i in 1:length(numberOfNodes))
  {
    for (j in 1:length(ratioOfEdgesPresence))
    {
      if (identical(matrixPopulation[[i]][[j]], matrix))
      {
        return(c(i, j))
      }
    }
  }
}

generateCondition <- function(matrix)
{
  iter <- 1
  conditions <- length(skewness) * length(numberOfObservations) * numberOfSamples
  matrixList <- vector(mode = "list", length = length(skewness))
  
  for (i in 1:length(skewness))
  {
    matrixList[[i]] <- vector(mode = "list", length = length(numberOfObservations))
    for (j in 1:length(numberOfObservations))
    {
      matrixList[[i]][[j]] <- vector(mode = "list", length = numberOfSamples)
      for (k in 1:numberOfSamples)
      {
        cat("Started", iter, "/", conditions, "conditions.\n")
        
        data <- generateData(input = matrix,
                             skewness = skewness[i],
                             n = numberOfObservations[j])
        if(is.positive.definite(qgraph::cor_auto(data, verbose = FALSE)))
        {
          matrixList[[i]][[j]][[k]] <- estimateNetwork(data,
                                                       default = "EBICglasso",
                                                       corMethod = "cor_auto",
                                                       tuning = 0.5,
                                                       refit = TRUE,
                                                       verbose = FALSE)$graph
        }
        else
        {
          errorsCounter <<- errorsCounter + 1
          matrixInfo <- identifyMatrix(matrix)
          matrixInfo <- append(matrixInfo, c(i, j, k, iter))
          errorsNonPositive[[errorsCounter]] <<- matrixInfo
          matrixList[[i]][[j]][[k]] <- NA
        }
        
        iter <- iter + 1
      }
    }
  }
  return(matrixList)
}


# 9. Estimate network based on the sample data ---------------------------------

errorsNonPositive <- list()
errorsCounter <- 0

set.seed(1)

matrixSample11 <- generateCondition(matrixPopulation[[1]][[1]])
save(matrixSample11, file = "matrixSample11.RData")
matrixSample12 <- generateCondition(matrixPopulation[[1]][[2]])
save(matrixSample12, file = "matrixSample12.RData")
matrixSample13 <- generateCondition(matrixPopulation[[1]][[3]])
save(matrixSample13, file = "matrixSample13.RData")

matrixSample1 <- list(matrixSample11, matrixSample12, matrixSample13)
save(matrixSample1, file = "matrixSample1.RData")

gc()

set.seed(2)

matrixSample21 <- generateCondition(matrixPopulation[[2]][[1]])
save(matrixSample21, file = "matrixSample21.RData")
matrixSample22 <- generateCondition(matrixPopulation[[2]][[2]])
save(matrixSample22, file = "matrixSample22.RData")
matrixSample23 <- generateCondition(matrixPopulation[[2]][[3]])
save(matrixSample23, file = "matrixSample23.RData")

matrixSample2 <- list(matrixSample21, matrixSample22, matrixSample23)
save(matrixSample2, file = "matrixSample2.RData")

gc()

set.seed(3)

matrixSample31 <- generateCondition(matrixPopulation[[3]][[1]])
save(matrixSample31, file = "matrixSample31.RData")
matrixSample32 <- generateCondition(matrixPopulation[[3]][[2]])
save(matrixSample32, file = "matrixSample32.RData")
matrixSample33 <- generateCondition(matrixPopulation[[3]][[3]])
save(matrixSample33, file = "matrixSample33.RData")

matrixSample3 <- list(matrixSample31, matrixSample32, matrixSample33)
save(matrixSample3, file = "matrixSample3.RData")

gc()

set.seed(4)

matrixSample41 <- generateCondition(matrixPopulation[[4]][[1]])
save(matrixSample41, file = "matrixSample41.RData")
matrixSample42 <- generateCondition(matrixPopulation[[4]][[2]])
save(matrixSample42, file = "matrixSample42.RData")
matrixSample43 <- generateCondition(matrixPopulation[[4]][[3]])
save(matrixSample43, file = "matrixSample43.RData")

matrixSample4 <- list(matrixSample41, matrixSample42, matrixSample43)
save(matrixSample4, file = "matrixSample4.RData")

gc()

set.seed(5)

matrixSample51 <- generateCondition(matrixPopulation[[5]][[1]])
save(matrixSample51, file = "matrixSample51.RData")
matrixSample52 <- generateCondition(matrixPopulation[[5]][[2]])
save(matrixSample52, file = "matrixSample52.RData")
matrixSample53 <- generateCondition(matrixPopulation[[5]][[3]])
save(matrixSample53, file = "matrixSample53.RData")

matrixSample5 <- list(matrixSample51, matrixSample52, matrixSample53)
save(matrixSample5, file = "matrixSample5.RData")

gc()

matrixSample <- list(matrixSample1, matrixSample2,
                     matrixSample3, matrixSample4,
                     matrixSample5)

save(matrixSample, file = "matrixSample.RData")


# a copy from the bootnet package
cor0 <- function(x,y){
  if (all(is.na(x)) || all(is.na(y))){
    return(NA)
  }
  
  if (sd(x,na.rm=TRUE)==0 | sd(y,na.rm=TRUE) == 0){
    return(0)
  }
  
  return(cor(x,y,use="pairwise.complete.obs"))
}


simulationResults <- data.frame(ID = numeric(),
                                numberOfNodes = factor(),
                                ratioOfEdgesPresence = factor(),
                                skewness = factor(),
                                numberOfObservations = factor(),
                                sampleNumber = factor(),
                                isEmpty = logical(),
                                correctModel = logical(),
                                sensitivity = numeric(),
                                specificity = numeric(),
                                correlation = numeric(),
                                meanBias = numeric(),
                                sdBias = numeric(),
                                maxBias = numeric(),
                                meanEdgeWeightDiff = numeric(),
                                maxEdgeWeightDiff = numeric(),
                                maxFalseEdgeWeight = numeric())

iter <- 1
conditions <- length(numberOfNodes) * length(ratioOfEdgesPresence) *
  length(skewness) * length(numberOfObservations) * numberOfSamples

for (nNodes in 1:length(numberOfNodes))
{
  for (ratio in 1:length(ratioOfEdgesPresence))
  {
    for (skew in 1:length(skewness))
    {
      for (nObs in 1:length(numberOfObservations))
      {
        for (nSample in 1:numberOfSamples)
        {
          cat("Calculates", iter, " / ", conditions, "\n")
          
          estNet <- matrixSample[[nNodes]][[ratio]][[skew]][[nObs]][[nSample]]
          trueNet <- matrixPopulation[[nNodes]][[ratio]]
          
          # Estimated edges:
          est <- estNet[upper.tri(estNet)]
          # Real edges:
          real <- trueNet[upper.tri(trueNet)]
          
          # True positives:
          TruePos <- sum(est != 0 &  real != 0)
          
          # False pos:
          FalsePos <- sum(est != 0 & real == 0)
          
          # True Neg:
          TrueNeg <- sum(est == 0 & real == 0)
          
          # False Neg:
          FalseNeg <- sum(est == 0 & real != 0)


          simulationResults <- simulationResults %>%
            add_row(ID = iter,
                    numberOfNodes = as.factor(paste0(numberOfNodes[nNodes], " nodes")),
                    ratioOfEdgesPresence = as.factor(paste0(ratioOfEdgesPresence[ratio] * 100, "% of edges")),
                    skewness = as.factor(skewness[skew]),
                    numberOfObservations = as.factor(numberOfObservations[nObs]),
                    sampleNumber = as.factor(nSample),
                    isEmpty = all(matrixSample[[nNodes]][[ratio]][[skew]][[nObs]][[nSample]]==0),
                    correctModel = all((est == 0) == (real == 0)),
                    sensitivity = TruePos / (TruePos + FalseNeg),
                    specificity = TrueNeg / (TrueNeg + FalsePos),
                    correlation = cor0(est, real),
                    meanBias = mean(abs(est - real)),
                    sdBias = sd(est - real),
                    maxBias = max(abs(est - real)),
                    meanEdgeWeightDiff = (mean(est)-mean(real))/mean(real),
                    maxEdgeWeightDiff = ((est[real == max(real)])-max(real))/max(real),
                    maxFalseEdgeWeight = max(abs(est[real==0 & est!=0])))
          iter <- iter + 1
        }
      }
    }
  }
}

save(simulationResults, file = "simulationResults.RData")

limits01 <- c(0, 1)
calculatedEmpty <- simulationResults %>% group_by(numberOfNodes,
                                      ratioOfEdgesPresence,
                                      skewness,
                                      numberOfObservations) %>%
  summarise(emptyPercentage = mean(isEmpty))


jpeg("figures/Fig2.jpeg", units = "in", width=9, height=6, quality = 100, res=300)
ggplot(data = distDF,
       aes(x = obs)) +
  geom_bar((aes(y = (..count..)/1000)),
           color = "black",
           fill = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  facet_wrap(~factor(type, levels = c("symmetric", "skewed", "highly skewed"))) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x = "Response category", y = "Percentage of responses")
dev.off()

jpeg("figures/Fig3.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = calculatedEmpty,
       aes(x = numberOfObservations,
           y = emptyPercentage,
           fill = skewness)) +
  geom_bar(position="dodge",
           stat="identity",
           color = "black") +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L), breaks = seq(limits01[1], limits01[2], by = 0.1)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Percentage of empty networks", color = "Distribution")
dev.off()

jpeg("figures/Fig4.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = sensitivity,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(limits = limits01, breaks = seq(limits01[1], limits01[2], by = 0.1)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Sensitivity", color = "Distribution")
dev.off()

jpeg("figures/Fig5.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = specificity,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(limits = limits01, breaks = seq(limits01[1], limits01[2], by = 0.1)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Specificity", color = "Distribution")
dev.off()

jpeg("figures/Fig6.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = correlation,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(limits = limits01, breaks = seq(limits01[1], limits01[2], by = 0.1)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Correlation", color = "Distribution")
dev.off()

jpeg("figures/Fig7.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = meanBias,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  scale_y_continuous(limits = c(0, 0.35)) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Mean bias", color = "Distribution")
dev.off()

jpeg("figures/Fig8.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = sdBias,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  scale_y_continuous(limits = c(0, 0.35)) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Standard deviation of bias", color = "Distribution")
dev.off()

jpeg("figures/Fig9.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = maxBias,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Max bias", color = "Distribution")
dev.off()

jpeg("figures/Fig10.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = meanEdgeWeightDiff,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Difference between mean edge weight", color = "Distribution")
dev.off()

jpeg("figures/Fig11.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = maxEdgeWeightDiff,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Difference between max edge weight", color = "Distribution")
dev.off()

jpeg("figures/Fig12.jpeg", units = "in", width=10, height=10, quality = 100, res=300)
ggplot(data = simulationResults,
       aes(x = numberOfObservations,
           y = maxFalseEdgeWeight,
           fill = skewness)) +
  geom_boxplot() +
  scale_fill_grey(start = 0.3, end = 1) +
  facet_grid(numberOfNodes ~ ratioOfEdgesPresence) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top") +
  labs(x = "Number of observations", y = "Max false edge weight", color = "Distribution")
dev.off()

