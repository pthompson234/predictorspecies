##Predictor species: Improving assessments of rare species occurrence by modelling environmental co-responses
##Sample R Code
##Peter Thompson, Phillip Staniczenko, William Fagan
##January 16, 2020

##This file will run one iteration of the Bayesian Network analysis conducted in Thompson et al. (2020), Ecology and Evolution.
##Here four different models are compared on simulated species presence/absence data along with simulated environmental data
##While this is simulated data, real-life data of similar nature can be obtained from https://datadryad.org/stash/dataset/doi:10.5061/dryad.g1pk3
##This data is also randomized, so do not expect the eGLM+BN to pull out more accurate results, but this can easily be applied to any other data

## setwd(wherever this file is located)

print("Beginning model development")

####Load packages####
#may have to install BiocManager and install some packages through that if not done already
require(gRain)
require(MASS) #for multivariate normal random variables
require(pROC) #for AUC evaluations

####Functions for Bayesian network analysis####
GetBinaryORFacInhResponse <- function(netcontributingunits, maxunits, p){
  # OR response relating number of faciliatory and/or inhibitory species present to probability of focal species presence
  if(netcontributingunits > 0){
    stepfactor <- 1
  }
  else if(netcontributingunits < 0){
    stepfactor <- -1
  }
  else{
    stepfactor <- 0
  }
  if(p > 0.5){
    return(p + stepfactor * (1 - p))
  }
  else {
    return(p + stepfactor * p)
  }
}
getNames <- function(prey){
  # get rownames of a state table defined by which species are present and which are not
  n <- length(prey)
  ch <- matrix("", 2^n, n)
  status <- rep("", 2^n)
  for(i in 0:((2^n) - 1)){
    x <- i
    for(j in 1:n){
      if(x%%2 == 0){
        ch[i+1, j] <- as.character(prey[j])
        x <- floor(x/2)
      } else {
        ch[i+1, j] <- "-"
        x <- floor(x/2)
      }
      status[i + 1] <- paste(status[i + 1], " ", ch[i + 1, j], sep='')
    }
  }
  return(status)
}
buildTable <- function(NET, Pi, model, alpha=1, beta=1) {
  # form input for getMarginals_gRain()  including state tables for each species
  use_sf <- 1  # USER-DEFINED OPTION
  if (use_sf==1){significant_figures <- 7}
  PiVector <- Pi
  S <- length(PiVector)
  Table <- list()
  Table$PiVector <- PiVector
  Table$M <- NET
  for (i in 1:S){
    variable_name <- paste("V", i, sep='')
    prey <- which(NET[, i] != 0)
    n <- length(prey)
    if (n==0){
      # species is basal
      X <- matrix(0, 1, 1)
      if (use_sf==1) {
        X[1, 1] <- signif(PiVector[i], digits=significant_figures)
      } else {
        X[1, 1] <- PiVector[i]
      }
      Table[[variable_name]] <- X
    } else {
      # species is a consumer
      X <- matrix(0, 2^n, 1)
      rownames(X) <- getNames(prey)
      if(model == "binaryplusminus_OR"){
        posEdges <- sum(NET[,i] > 0)  # Total number of facilitatory dependencies for i
        negEdges <- sum(NET[,i] < 0)  # Total number of inhibitory dependencies for i
        maxunits <- max(posEdges, negEdges)
        for(j in 1:(2^n)) {
          tmp <- j - 1
          v <- numeric(n)  # create a vector of length n
          ## Determine which resources to consider for case j
          for(k in 1:n){
            if((tmp %% 2) == 0){
              v[k] = 1
              tmp <- floor(tmp/2)
            } else {
              v[k] = 0
              tmp <- floor(tmp/2)
            }
          }
          netcontributingunits <- NET[NET[, i] != 0, i]%*%v  # Total number of facilitatory minus inhibitory resources considered in case j
          if(use_sf==1){
            X[j, 1] <- signif(GetBinaryORFacInhResponse(netcontributingunits, maxunits, PiVector[i]), digits=significant_figures)
          }
          else{
            X[j, 1] <- GetBinaryORFacInhResponse(netcontributingunits, maxunits, PiVector[i])
          } 
        }
      } else {
        stop("You must specify a response form: binaryplusminus_linear")
      }	
      Table[[variable_name]] <- X
    }
  }
  Table$model <- model
  return(Table)
}
getPreys <- function(v){
  # get string of names of preys of a consumer species with "+" in between
  n <- length(v)
  p <- sum(v > 0)
  prey_str <- ''
  for(i in 1:n){
    if(p > 1){
      if(v[i] > 0){
        prey_str <- paste(prey_str, 'V', i, "+", sep='')
        p <- p - 1
      }
    } else{
      break
    }	
  }
  for(j in i:n){
    if(v[j] > 0){
      prey_str <- paste(prey_str, 'V', j, sep='')
    }
  }
  return(prey_str)
}
getMarginals <- function(Table, file="test.R"){
  # automatically produce an R script file for running gRain
  write("require(gRain)", file=file)
  write("temp_fn <- function(){", file, append=TRUE)
  sp_levels <- "levels=c('absent', 'present'))"
  S <- length(Table$PiVector)
  for(i in 1:S){
    sp_name <- paste("V", i, sep='')
    if(dim(Table[[sp_name]])[1]==1){  # basal species
      X <- Table[[sp_name]]
      command <- paste("V", i, "<-cptable(~V", i, ", values=c(", 1 - X[1,1], ",", X[1,1], "), ", sp_levels, sep='')
      write(command, file=file, append=TRUE)
    } else {  # consumers 
      X <- Table[[sp_name]]
      nX <- dim(X)[1]
      command <- paste("V", i, "<-cptable(~V", i, "|", sep='')
      prey_str <- getPreys(abs(Table$M[, i]))
      command <- paste(command, prey_str, ", values=c(", sep='')
      for(j in 1:(nX - 1)){
        command <- paste(command, 1 - X[nX - j + 1, 1], ",", X[nX - j + 1, 1], ",", sep='')  # l-to-r: all prey present to all extinct
      }
      command <- paste(command, 1 - X[1, 1], ",", X[1, 1], "), ", sp_levels, sep='')
      write(command, file, append=TRUE)
    }
  }
  command2 <- paste("plist<-compileCPT(list(", sep='')
  for(i in 1:(S-1)){
    command2 <- paste(command2, "V", i, ",", sep='')
  }
  command2 <- paste(command2, "V", S, "))", sep='')
  write(command2, file, append=TRUE)
  write("BN<-grain(plist)", file, append=TRUE)
  write("Results_tmp<-querygrain(BN)", file, append=TRUE)
  write("return(Results_tmp)", file, append=TRUE)
  write("}", file, append=TRUE)
  
  # Run the script
  source(file)
  r <- temp_fn()
  results <- list()
  # Re-orgainse the sequence in the output
  for(i in 1:S){
    sp_name <- paste("V", i, sep='')
    results[[sp_name]] <- r[[sp_name]]
  }
  return(results)	
}

####Load in external data####

YFull = read.csv("presenceabsencedata.csv")
XFull = read.csv("envdata.csv")

m = length(YFull[1,]) # number of species
nFull = length(YFull[,1]) # number of sample locations
k = length(XFull[1,]) + 1 # number of environmental covariates, add 1 for intercept

print("Packages and files have been loaded")

####Divide the random data into training and test partitions####
partitionSize = 0.5 # size of training partition for data; varies from 0.25 to 0.75 in empirical study
n = floor(partitionSize * nFull) # "floor" in case partitionSize > 1
randomSort = sample(1:nFull) # randomly sorts the location numbers to select a training partition
trainPartition = sort(randomSort[1:n])
testPartition = sort(randomSort[(n + 1):nFull])
##Store data as matrices
Y = matrix(nrow = n, ncol = m, 0)
YT = matrix(nrow = (nFull - n), ncol = m, 0)
X = matrix(nrow = n, ncol = k, 0)
XT = matrix(nrow = (nFull - n), ncol = k, 0)
##represents Y-intercept (for prediction purposes)
X[,1] = 1
XT[,1] = 1
##Fill matrices with species presence/absence data as well as environmental data
for (i in 1:n) {
  for (j in 1:m) {
    Y[i, j] = as.numeric(as.character(YFull[trainPartition[i],j]))
  }
}
for (i in 1:(nFull - n)) {
  for (j in 1:m) {
    YT[i, j] = as.numeric(as.character(YFull[testPartition[i],j]))
  }
}
for (i in 1:n) {
  for (j in 2:k) {
    X[i,j] = as.numeric(as.character(XFull[trainPartition[i],(j-1)]))
  }
}
for (i in 1:(nFull - n)) {
  for (j in 2:k) {
    XT[i,j] = as.numeric(as.character(XFull[testPartition[i],(j-1)]))
  }
}
##convert species data to binary presence/absence
Y = sign(Y)
YT = sign(YT)

print("Data has been randomly partitioned into training and test sets")

####The eGLM: A Baseline Model Incorporating Only Abiotic Conditions####

print("Beginning model #1: eGLM")

##Get B by running logistic regression on each species
B = matrix(nrow = k, ncol = m)
X = as.data.frame(X) # convert to data.frame for model fitting
for (i in 1:m) {
  # fit a logistic GLM for each individual species
  model = glm(Y[,i] ~ . , data = X[,2:k], family = binomial(link = "logit"))
  B[,i] = model$coefficients
}
X = as.matrix(X) # back to a matrix so we can apply matrix multiplication
M = XT %*% B # predictions of each regression coefficient for test data
OUTPUT_1_NOLOGIT = M
for (i in 1:(length(M[,1]))) {
  for (j in 1:(length(M[1,]))) {
    if (abs(M[i,j]) > 700) {
      # if any values are excessively large, correct (won't really matter for logit)
      M[i,j] = sign(M[i,j]) * 700
    }
  }
}
M = apply(M, 2, function(x) {exp(x) / (1 + exp(x))}) # get inverse logit link for binomial predictions
OUTPUT_1 = M

####Obtaining a matrix of Pearson correlations between species####

CORR = matrix(nrow = m, ncol = m) # correlations for all species from the entire data set
##We can basically think of this as known interactions between species, and if the study system is well-understood, we don't need a correlation matrix
for (i in 1:m) {
  for (j in 1:m) {
    CORR[i, j] = cor(sign(as.numeric(as.character(YFull[,i]))), sign(as.numeric(as.character(YFull[,j]))))
  }
}

#this is random data so none of the species are actually correlated; therefore we use a low threshold
CORR_THRESHOLD = 0.25 # may vary depending on the nature of the data (i.e., this threshold may overpopulate the BN for some study systems)
##Read in a hierarchy of influence for the BN in order to direct edges
#In the BN "rows influence columns" so to find the number of incoming edges a species has look at its column sum
hierarchy = as.numeric(read.csv("hierarchy.csv"))
####Constructing the BN based on the correlation matrix####
BN = matrix(nrow = m, ncol = m)
for (i in 1:m) {
  for (j in 1:m) {
    if (i == j || abs(CORR[i, j]) < CORR_THRESHOLD || hierarchy[i] < hierarchy[j]) {
      #if the edge is on the wrong side of hierarchy, or relationship isn't strong enough
      BN[i, j] = 0
    } else {
      #either 1 or -1
      BN[i, j] = sign(CORR[i, j])
    }
  }
}

####The sGLM: A Basic Implementation of a JSDM####

print("Beginning model #2: sGLM")

##Get B by running logistic regression on each species
M = matrix(nrow = length(testPartition), ncol = m)
for (i in 1:m) {
  influences = as.vector(which(BN[,i] != 0)) # finding any species that influence species i
  #Manually buildling the GLM formula based on what species are included as covariates
  form1 = paste("Y[,",i,"] ~",sep = "")
  form2 = paste("X[,",2:k,"]",sep = "", collapse = " + ")
  form3 = paste("Y[,",influences,"]",sep = "",collapse = " + ")
  if (length(influences) > 0) {
    formTOT = paste(form1, form2, "+", form3, sep = " ")
  } else {
    formTOT = paste(form1, form2, sep = " ")
  }
  #notice that we don't make any assumptions about "sign" of the influence, it's just included and fit like other covariates
  model = glm(formula = as.formula(formTOT), family = binomial(link = "logit"))
  BI = model$coefficients
  BI[is.na(BI)] = 0
  modXT = matrix(nrow = length(testPartition), ncol = k + length(influences))
  modXT[,1:k] = XT
  if (length(influences) > 0) {
    modXT[,(k+1):(length(influences)+k)] = YT[,influences]
  }
  M[,i] = modXT %*% BI
}
for (i in 1:(length(M[,1]))) {
  for (j in 1:(length(M[1,]))) {
    if (abs(M[i,j]) > 700) {
      M[i,j] = sign(M[i,j]) * 700
    }
  }
}
M = apply(M, 2, function(x) {exp(x) / (1 + exp(x))}) # get inverse logit link for binomial predictions
OUTPUT_2 = M

####The eGLM+BN: Our Novel Method Incorporating Bayesian Networks####

print("Beginning model #3: eGLM+BN")

M = matrix(nrow = length(testPartition), ncol = m)
for (loc in 1:length(testPartition)) {
  for (i in 1:m) {
    if (sum(abs(BN[,i])) == 0) {
      # if there are no incoming BN edges to i our BN result is the same as the eGLM result
      M[loc, i] = OUTPUT_1[loc, i]
    } else {
      prior = YT[loc,]
      prior[i] = OUTPUT_1[loc, i]
      #assume known presences in all species except i (as we do in sGLM), improving our prior guess through the BN
      table = buildTable(BN, prior, model = "binaryplusminus_OR")
      fileName = paste0("out_temp_file.R")
      marginals = unlist(getMarginals(table, fileName))[seq(2, 2 * m, 2)]
      M[loc,i] = marginals[i]
    }
  }
}
OUTPUT_3 = M

####Evaluating Model Accuracy Using AUC####

print("Model fitting complete. Beginning AUC analysis")

AUC_1 = numeric(m)
AUC_2 = numeric(m)
AUC_3 = numeric(m)

AUC_ALL = matrix(0, nrow = 3, ncol = m)
for (i in 1:m) {
  if (sum(YT[,i]) == 0 || sum(YT[,i]) == length(YT[,i])) {
    AUC_ALL[1, i] = NA
    AUC_ALL[2, i] = NA
    AUC_ALL[3, i] = NA
  } else {
    AUC_ALL[1, i] = auc(roc(YT[,i], OUTPUT_1[,i]))
    AUC_ALL[2, i]= auc(roc(YT[,i], OUTPUT_2[,i]))
    AUC_ALL[3, i] = auc(roc(YT[,i], OUTPUT_3[,i]))
  }
}

write.csv(AUC_ALL, "AUC_OUTPUTS.csv")
