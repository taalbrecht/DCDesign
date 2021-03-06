#' Linear Model D-Optimality. Review for deprecation and removal.
#'
#' @param CurrentMatrix matrix - A model matrix with appropriate factor numerical coding as created by the model.matrix function. Data.frame formats not accepted
#' @param returncov logical - Whether the variance-covariance matrix should be returned (requires additional calculation time)
#'
#' @return Either the D-Optimality of the model if returncov is FALSE or a list with elements d_eff (D-Optimality) and vcov (covariance matrix) if returncov is TRUE
#' @description Calculates determinant of supplied model and raises it to the power of 1 / ncol(ModelMatrix) and divides the result by the number of rows in ModelMatrix.
#' @export
#'
#' @examples
d_efficiencysimple <- function(CurrentMatrix, returncov = FALSE){

  #Calculate information matrix
  infomat <- t(CurrentMatrix)%*%CurrentMatrix
  #Calculate determinant of information matrix
  det_calc <- det(infomat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Calculate d-efficiency
  d_eff <- ((det_calc)^(1/(ncol(CurrentMatrix))))/nrow(CurrentMatrix)

  if(returncov == TRUE){

    #Calculate variance-covariance matrix
    vocvmat <- tryCatch(solve(infomat), error = function(x) diag(x = Inf, nrow = 2, ncol = 2))
    output <- list(d_eff = d_eff, vcov = vocvmat)

  }else{output <- d_eff}

  #Return objective function and determinants for both current models
  return(output)}


#' Linear model D-Optimality update function. Review for deprecation and removal.
#'
#' @param CurrentMatrix matrix - A model matrix with appropriate factor numerical coding as created by the model.matrix function. Data.frame formats not accepted
#'
#' @return The D-Optimality of the model
#' @description Calculates D-Optimality of the supplied model matrix by calculating the determinant and raising it to the power of 1 / ncol(ModelMatrix) then dividing the result by the number of rows in ModelMatrix.
#' @export
#'
#' @examples
d_efflinearupdate <- function(CurrentMatrix){

  #Calculate information matrix
  infomat <- t(CurrentMatrix)%*%CurrentMatrix

  #Calculate determinant of information matrix
  det_calc <- det(infomat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Return d-efficiency
  return(((det_calc)^(1/(ncol(CurrentMatrix))))/nrow(CurrentMatrix))}


#' D-efficiency calculation
#'
#' @param CurrentMatrix either a matrix or data.frame of base variables only (ex: Col for A, B, C) that will be converted to a model matrix using input_formula
#' @param det_ref reference determinant for calculation of d_efficiency. Should be max attainable determinant of info matrix
#' @param input_formula one-sided formula (inputs only) for model to calculate d efficiency from. Format = ~ A + B + A:B + I(A/B^2), etc
#' required if CurrentMatrix is a data frame of base variables that needs to be converted to a model matrix
#' @param Input_range matrix or data frame listing the range of columns used for standardization. Column names must match input_formula term names.
#' format is matrix or data frame with column in Input_range matching each column name in CurrentMatrix and minimum value in first row and maximum value in second row.
#'
#' @return vector(d_eff,det_calc) where the first element is the d-efficiency and the second element is the determinanat of info matrix created using CurrentMatrix and input_formula
#' @description Calculates determinant ratio based on supplied det_ref. Will convert a base data.frame/matrix to a model matrix and standardize based on supplied input_formula and Input_range.
#' @export
#'
#' @examples
d_efficiency <- function(CurrentMatrix, det_ref, input_formula, Input_range){

  #Convert input data to data.frame in case it was passed as matrix or array
  CurrentMatrix <- data.frame(CurrentMatrix)

  #Get model matrix
  modelmat <- stats::model.matrix(input_formula,data = CurrentMatrix)

  #Standardize model matrix
  modelmat <- standardize_cols(modelmat, colnames(modelmat[,2:ncol(modelmat)]), Input_range = Input_range)

  #Calculate determinant of information matrix
  det_calc <- det(t(modelmat)%*%modelmat)

  #Set determinants equal to zero if less than zero due to R not being able to precisely calculate determinant
  if(det_calc < 0){det_calc <- 0}

  #Calculate ratio of current determinant to optimal determinanat for additive and mechanistic model
  d_eff <- ((det_calc/det_ref)^(1/(ncol(modelmat))))

  #Construct return vector with named elements
  returnvect <- c(d_eff, det_calc)
  names(returnvect) <- c("D efficiency", "Info Matrix Determinant")

  #Return objective function and determinants for both current models
  return(returnvect)}


#' D-Error for multinomial logit model
#'
#' @param CurrentMatrix model matrix as created by model.matrix function with all continuous parameters centered and standardized
#'   and appropriate factor numerical coding (generated using contr.sum for factors) as created by the model.matrix function.
#'   Data.frame formats not accepted
#'   Opt out choices, if included, should be coded as a row of all zeroes for every parameter
#' @param altvect vector with integer corresponding to each row that indicates which choice set it is a member of
#' @param paramestimates estimates for each effect (column) of model matrix sized corresponding to standardized model matrix.
#'   If not supplied, parameter estimates will be assumed equal to zero for all parameters
#' @param returninfomat whether to return the covariance matrix as well. Default is FALSE
#'
#' @return list containing:
#'   \item{d_eff}{D-Optimality of the design}
#'   \item{info_mat}{Information matrix of the design}
#'   \item{p_var}{The sum of the product of the probability of selection of each choice set in the design.}
#'   \item{probvect}{A vector of the probability of selection of each profile in the design.}
#' @description Calculates d-error of supplied model matrix. Will also calculate probability centered d-error if vector of parameter estimates is supplied
#' @export
#'
#' @examples
d_effchoice <- function(CurrentMatrix, altvect, paramestimates = NULL, returninfomat = FALSE){

  #get all unique alternate names
  altnames <- unique(altvect)

  #check supplied parameters and set = 0 if not supplied
  if(is.null(paramestimates)){
    paramestimates <- rep(0, ncol(CurrentMatrix))
  }

  #Get position of intercept column
  iceptcol <- grep("(Intercept)", colnames(CurrentMatrix), fixed = TRUE)

  #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
  if(length(iceptcol) == 1){
    CurrentMatrix[,iceptcol] <- rep(c(0, rep(1, length(altvect)/length(unique(altvect)) - 1)), times = length(unique(altvect)))
  }

  info_mat=matrix(rep(0,ncol(CurrentMatrix)*ncol(CurrentMatrix)), ncol(CurrentMatrix), ncol(CurrentMatrix))
  # compute exp(design matrix times initial parameter values)
  exputilities=exp(CurrentMatrix%*%paramestimates)

  #Initialize vector that will be sum of product of each set probability (variance)
  p_var <- 0

  #Initialize vector that holds calculated probability of each alternative being selected
  probvect <- rep(0,nrow(CurrentMatrix))

  #Loop over each choice set
  for (k_set in 1:length(altnames)) {
    # select row numbers corresponding to current loop alternatives in the choice set
    alternatives= which(altvect == altnames[k_set])
    # obtain vector of choice shares within the choice set
    p_set=exputilities[alternatives]/sum(exputilities[alternatives])
    # also put these probabilities on the diagonal of a matrix that only contains zeros
    p_diag=diag(p_set)
    # compute middle term P-pp'
    middle_term<-p_diag-p_set%o%p_set
    # pre- and postmultiply with the Xs from the design matrix for the alternatives in this choice set
    full_term<-t(CurrentMatrix[alternatives,])%*%middle_term%*%CurrentMatrix[alternatives,]
    # Add contribution of this choice set to the information matrix
    info_mat<-info_mat+full_term

    #Calculate product of all probabilities in set and add to p_var
    p_var <- p_var + prod(p_set)

    #Enter all calculated probabilities for this set into the probability vector
    probvect[alternatives] <- p_set

  }
  #get the inverse of the information matrix (i.e., gets the variance-covariance matrix)
  #Use "try" wrapper to prevent unsolvable matrices from crashing. Return 2x2 diagonal infinite matrix on crash
  #sigma_beta<- tryCatch(solve(info_mat,diag(ncol(CurrentMatrix))), error = function(x) diag(x = Inf, nrow = ncol(CurrentMatrix), ncol = ncol(CurrentMatrix)))


  #Calculate determinant
  det_calc <- det(info_mat)
  #If determinant is negative (should not be possible but sometimes happens), return zero to prevent an error
  if(det_calc < 0){det_calc <- 0}

  if(returninfomat == TRUE){

    output <- list(d_eff = det_calc^(1/ncol(CurrentMatrix)), info_mat = info_mat, p_var = p_var, probvect = probvect)

  }else{output <- det_calc^(1/ncol(CurrentMatrix))}

  #Return objective function and determinants for both current models
  #return(list(d_eff = det(sigma_beta)^(-1/ncol(CurrentMatrix)), vcov = sigma_beta))}
  return(output)}


#' D-efficiency update function for single question for multinomial logit model
#'
#' @param CurrentMatrix subsection of model matrix for one question only. Should include all alternatives for that questionas created by the model.matrix function with all continuous parameters centered and standardized
#'   and appropriate factor numerical coding (generated using contr.sum for factors) as created by the model.matrix function.
#'   Data.frame formats not accepted
#'   Opt out choices, if included, should be coded as a row of all zeroes for every parameter
#' @param paramestimates standardized effect estimates for each effect (column) of model matrix corresponding to standardized model matrix.
#' @param info_mat the information matrix for all other questions in the model matrix
#'
#' @return
#'   \item{d_eff}{D-Optimality of the design}
#'   \item{p_var}{The product of all selection probabilities for this choice set}
#'   \item{info_mat}{Information matrix of the choice set.}
#'   \item{p_set}{A vector of the probability of selection of each profile in the choice set.}
#' @description Fast d-efficiency update function for use during search algorithms to reduce search time. Requires supplied information matrix for all questions not included
#' @export
#'
#' @examples
d_effchoiceupdate <- function(CurrentMatrix, paramestimates, info_mat = 0){

  #Get position of intercept column
  iceptcol <- grep("(Intercept)", colnames(CurrentMatrix), fixed = TRUE)

  #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
  if(length(iceptcol) == 1){
    CurrentMatrix[,iceptcol] <- c(0, rep(1, nrow(CurrentMatrix) - 1))
  }


  # compute exp(design matrix times initial parameter values)
  exputilities <- c(exp(CurrentMatrix%*%paramestimates))

  # obtain vector of choice shares within the choice set
  p_set <- exputilities/sum(exputilities)

  #Calculate product of all probabilities in set for utility balance use
  p_var <- prod(p_set)

  # calculate information matrix of this choice set and add it to the info_matrix
  info_mat<-info_mat + t(CurrentMatrix)%*%(diag(p_set)-p_set%o%p_set)%*%CurrentMatrix

  #Calculate determinant
  det_calc <- det(info_mat)

  #If determinant is negative (should not be possible but sometimes happens due to numerical processing), return zero to prevent an error
  if(det_calc < 0){det_calc <- 0}

  #Return Determinant
  return(list(d_eff = det_calc^(1/ncol(CurrentMatrix)), p_var = p_var, info_mat = info_mat, p_set = p_set))}


#' Matchup probability for choice tournament
#'
#' @param matchupframe data frame with the following columns:
#'   matrixrowid - unique number corresponding to the alternative ("team") as it progresses through the tournament.
#'   questionid - number of the matchup ("game").
#'   uniquesetid - unique number assigned to each potential matchup. There will be multiple potential matchups for each possible questionid
#'   sourcequestion - the previous questionid ("game") that the alternative ("team") must have won to arrive in this questionid ("game"). Should be 0 for the first questionid ("game") for each alternative ("team")
#' @param winprobs vector of win probabilities for each alternative in matchupframe. Must have length equal to number of rows in matchupframe
#' @param occurprobs (optional) vector of pre-existing chance of each matchup occurring. Should be 1 for any questionid that does not have a sourcequestion in matchupframe. If not passed, the entire occurrence probability vector will be generated from scratch.
#' @param updatevect (optional) vector of questionids in matchupframe that should be updated based on the supplied winprobs vector and previous occurprobs vector. If not passed, all occurrence probabilities will be updated.
#'
#' @return
#'   \item{occurprobs}{A vector of the probability of each choice set occurring.}
#' @description Calculates the probability of occurrence for every possible matchup in a choice tournament bracket.
#' @export
#'
#' @examples
matchupprobs <- function(matchupframe, winprobs, occurprobs = NA, updatevect = NA){

  #Initialize occurrence probability vector and update vector if occurrence probability vector was not provided
  if(is.na(occurprobs)[1]){
    occurprobs <- rep(1, times = nrow(matchupframe))
    updatevect <- unique(matchupframe$questionid[matchupframe$sourcequestion > 0])
  }
  #browser()
  #Initialize update vector if it was not provided
  if(is.na(updatevect)[1]){
    updatevect <- unique(matchupframe$questionid[matchupframe$sourcequestion > 0])
  }

  #Loop through all tournament questions
  for(i in updatevect){
    #Loop through all rows for this tournament set in matchup frame to calculate the probability of each alternative making it to this round of the tournament
    for(j in which(matchupframe$questionid == i)){

      temprows <- which((matchupframe$matrixrowid == matchupframe$matrixrowid[j]) &
                          (matchupframe$questionid == matchupframe$sourcequestion[j]))

      occurprobs[j] <- sum(occurprobs[temprows]*winprobs[temprows])

    }
    #Multiply chance of each alternative making it to this question to get probability of specific matchup occuring
    for(j in unique(matchupframe$uniquesetid[matchupframe$questionid == i])){

      occurprobs[matchupframe$uniquesetid == j] <- prod(occurprobs[matchupframe$uniquesetid == j])


    }
  }
  return(list(occurprobs = occurprobs))}


#' Objective function calculator for linear and multinomial logit models for different objective functions
#'
#' @param modmat model matrix as produced by using the \code{\link[stats]{model.matrix}} function on the whole formula to be optimized - needed for all objective functions
#' @param modparams named list of all other elements needed for objective function below. Named elements are as follows:
#' \describe{
#'   \item{objname}{name of objective function to use - needed for all objective functions}
#'   \item{momentmatrix}{moment matrix. Must be a p x p matrix where p is the number of columns in modmat - only needed for I optimality}
#'   \item{referencemat}{reference matrix used to calculate squared error between current design and reference matrix - only needed for SSE objective function}
#'   \item{priorvals}{vector of prior values used for choice based (multinomial logit) optimality criteria}
#'   \item{altvect}{vector of question set IDs used for choice based (multinomial logit) optimality criteria}
#'   \item{matchupframe}{choice tournament matrix generated by function \code{\link{gentourneybracket}}. Used for calculating optimality for tournament based choice optimality criteria}
#'   }
#'
#' @return The value of the objective function
#' @export
#'
#' @description Universal function used to calculate a wide variety of objective functions for linear and multinomial logit models.
#'   Typically called during the model search process to calculate the objective function on every iteration.
#'
#' @examples
objfunc <- function(modmat, modparams){

  #Linear model objective functions

  #Calculate D-Efficiency for specified models (On error, return Inf since goal is to minimize)
  if("D-Opt" == modparams$objname){
    output <- -tryCatch(det(t(modmat)%*%modmat)^(1/ncol(modmat))/nrow(modmat), error = function(x) -Inf)
  }


  #Calculate I-Optimality for specified models (On error, return Inf since goal is to minimize)
  if("I-Opt" == modparams$objname){
    output <- tryCatch(sum(diag(modparams$momentmatrix%*%solve(t(modmat)%*%modmat))), error = function(x) Inf)
  }

  ##Multinomial logit objective functions for choice designs - integrate first two options into getchoicemodinfoprobs if possible

  #Calculate multinomial choice model D-Efficiency for specified models (On error, return Inf since goal is to minimize)
  if("D-OptChoice" == modparams$objname){

    #get all unique alternate names
    altnames <- unique(modparams$altvect)

    #check supplied parameters and set = 0 if not supplied
    if(is.null(modparams$paramestimates)){
      modparams$paramestimates <- rep(0, ncol(modmat))
    }

    #Get position of intercept column
    iceptcol <- grep("(Intercept)", colnames(modmat), fixed = TRUE)
    if(length(iceptcol) == 1){
      #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
      modmat[,iceptcol] <- rep(c(0, rep(1, length(modparams$altvect)/length(unique(modparams$altvect)) - 1)), times = length(unique(modparams$altvect)))
    }

    info_mat=matrix(rep(0,ncol(modmat)*ncol(modmat)), ncol(modmat), ncol(modmat))
    # compute exp(design matrix times initial parameter values)
    exputilities=exp(modmat%*%modparams$paramestimates)

    #Initialize vector that will be sum of product of each set probability (variance)
    p_var <- 0

    #Initialize vector that holds calculated probability of each alternative being selected
    probvect <- rep(0,nrow(modmat))

    #Loop over each choice set
    for (k_set in 1:length(altnames)) {
      # select row numbers corresponding to current loop alternatives in the choice set
      alternatives= which(modparams$altvect == altnames[k_set])
      # obtain vector of choice shares within the choice set
      p_set=exputilities[alternatives]/sum(exputilities[alternatives])
      # also put these probabilities on the diagonal of a matrix that only contains zeros
      p_diag=diag(p_set)
      # compute middle term P-pp'
      middle_term<-p_diag-p_set%o%p_set
      # pre- and postmultiply with the Xs from the design matrix for the alternatives in this choice set
      full_term<-t(modmat[alternatives,])%*%middle_term%*%modmat[alternatives,]
      # Add contribution of this choice set to the information matrix
      info_mat<-info_mat+full_term

      #Calculate product of all probabilities in set and add to p_var
      p_var <- p_var + prod(p_set)

      #Enter all calculated probabilities for this set into the probability vector
      probvect[alternatives] <- p_set

    }

    #Calculate determinant
    output <- -tryCatch(det(info_mat)^(1/ncol(modmat))/length(unique(modparams$altvect)), error = function(x) -Inf)

  }

  #Calculate multinomial choice model I-Optimality for specified models (On error, return Inf since goal is to minimize)
  if("I-OptChoice" == modparams$objname){


    #get all unique alternate names
    altnames <- unique(modparams$altvect)

    #check supplied parameters and set = 0 if not supplied
    if(is.null(modparams$paramestimates)){
      modparams$paramestimates <- rep(0, ncol(modmat))
    }

    #Get position of intercept column
    iceptcol <- grep("(Intercept)", colnames(modmat), fixed = TRUE)
    if(length(iceptcol) == 1){
      #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
      modmat[,iceptcol] <- rep(c(0, rep(1, length(modparams$altvect)/length(unique(modparams$altvect)) - 1)), times = length(unique(modparams$altvect)))
    }


    info_mat=matrix(rep(0,ncol(modmat)*ncol(modmat)), ncol(modmat), ncol(modmat))
    # compute exp(design matrix times initial parameter values)
    exputilities=exp(modmat%*%modparams$paramestimates)

    #Initialize vector that will be sum of product of each set probability (variance)
    p_var <- 0

    #Initialize vector that holds calculated probability of each alternative being selected
    probvect <- rep(0,nrow(modmat))

    #Loop over each choice set
    for (k_set in 1:length(altnames)) {
      # select row numbers corresponding to current loop alternatives in the choice set
      alternatives= which(modparams$altvect == altnames[k_set])
      # obtain vector of choice shares within the choice set
      p_set=exputilities[alternatives]/sum(exputilities[alternatives])
      # also put these probabilities on the diagonal of a matrix that only contains zeros
      p_diag=diag(p_set)
      # compute middle term P-pp'
      middle_term<-p_diag-p_set%o%p_set
      # pre- and postmultiply with the Xs from the design matrix for the alternatives in this choice set
      full_term<-t(modmat[alternatives,])%*%middle_term%*%modmat[alternatives,]
      # Add contribution of this choice set to the information matrix
      info_mat<-info_mat+full_term

      #Calculate product of all probabilities in set and add to p_var
      p_var <- p_var + prod(p_set)

      #Enter all calculated probabilities for this set into the probability vector
      probvect[alternatives] <- p_set

    }

    #Calculate I-optimality criterion
    output <- tryCatch(sum(diag(modparams$momentmatrix%*%solve(info_mat))), error = function(x) Inf)

  }

  #Calculate multinomial choice model Bayesian D-Efficiency for specified models (On error, return Inf since goal is to minimize)
  if("D-OptChoiceBayes" == modparams$objname){

    #Get bayes info matrix
    info_mat <- bayesinfochoice(bayesmat = modparams$bayesmat, modmat = modmat, altvect = modparams$altvect)

    #Calculate determinant
    output <- -tryCatch(det(info_mat)^(1/ncol(modmat))/length(unique(modparams$altvect)), error = function(x) -Inf)

  }

  #Calculate multinomial choice model Bayesian I-Efficiency for specified models (On error, return Inf since goal is to minimize)
  if("I-OptChoiceBayes" == modparams$objname){

    #Get bayes info matrix
    info_mat <- bayesinfochoice(bayesmat = modparams$bayesmat, modmat = modmat, altvect = modparams$altvect)

    #Calculate I-optimality criterion
    output <- tryCatch(sum(diag(modparams$momentmatrix%*%solve(info_mat))), error = function(x) Inf)

  }

  #Calculate multinomial tournament choice model D-Efficiency for specified models (On error, return Inf since goal is to minimize)
  if("D-OptChoiceTournament" == modparams$objname){

    #Check to see if usewinprobs was specified in modparams. If so, use it, otherwise use default of TRUE
    if(!is.null(modparams$usepriorwinprobs)){
      usepriorwinprobs <- modparams$usepriorwinprobs
    }else{
      usepriorwinprobs <- TRUE
    }

    #Get information matrix weighted by likelihood of each tournament matchup occurring based on value of modparams$
    info_mat <- getchoicemodinfoprobs(modmat = modmat, paramestimates = modparams$paramestimates, matchupframe = modparams$matchupframe, usepriorwinprobs = usepriorwinprobs)$info_mat

    #Calculate determinant of info matrix
    output <- -tryCatch(det(info_mat)^(1/ncol(modmat))/length(unique(modparams$matchupframe$questionid)), error = function(x) -Inf)

  }

  #Calculate multinomial tournament choice model I-optimality for specified models (On error, return Inf since goal is to minimize)
  if("I-OptChoiceTournament" == modparams$objname){


    #Check to see if usewinprobs was specified in modparams. If so, use it, otherwise use default of TRUE
    if(!is.null(modparams$usepriorwinprobs)){
      usepriorwinprobs <- modparams$usepriorwinprobs
    }else{
      usepriorwinprobs <- TRUE
    }

    #Get information matrix weighted by likelihood of each tournament matchup occurring based on value of modparams$
    info_mat <- getchoicemodinfoprobs(modmat = modmat, paramestimates = modparams$paramestimates, matchupframe = modparams$matchupframe, usepriorwinprobs = usepriorwinprobs)$info_mat

    #Calculate I-optimality criterion
    output <- tryCatch(sum(diag(modparams$momentmatrix%*%solve(info_mat))), error = function(x) Inf)

  }

  ##Other objective functions

  #Calculate SSE for specified models
  if("SSE" == modparams$objname){

    #Calculate sum of squared error between current design and reference matrix
    output <- sum((modparams$referencemat - modmat)^2)

  }

  return(output)

}


#' Calculate Bayesian information matrix for discrete choice designs
#'
#' @param bayesmat matrix - matrix of points at which function should be returned. Each row is the evaluation point for that respective variable. Each column will be interpreted as a separate point.
#' @param modmat matrix - model matrix of single choice set to evaluate (produced using model.matrix(formula, design matrix))
#' @param altvect vector with integer corresponding to each row that indicates which choice set it is a member of
#'
#' @return The Bayesian information matrix of the design
#' @export
#'
#' @description Function to get bayesian information matrix for multinomial logit model based on provision of a random sample of points drawn from the prior parameter distribution
#'
#' @examples
bayesinfochoice <- function(bayesmat, modmat, altvect){

  #get all unique alternate names
  altnames <- unique(altvect)

  #Get position of intercept column
  iceptcol <- grep("(Intercept)", colnames(modmat), fixed = TRUE)

  #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
  if(length(iceptcol) == 1){
    modmat[,iceptcol] <- rep(c(0, rep(1, length(altvect)/length(unique(altvect)) - 1)), times = length(unique(altvect)))
  }

  ##Option 1

  #Calculate the exponentiation of the utility for each alternative under entire matrix of parameter values in x
  #(take transpose to make it easier to work with subsequent functions)
  #exputilities <- t(exp(modmat%*%x))
  exputilities <- exp(bayesmat%*%t(modmat))

  #Translate the exputilities the probability for each entry in the exponetiated utility matrix based on the sum of all utilities for that choice set
  #probmat <- exputilities

  #Create large flattened information matrix to store resulting information matrix
  infoflat <- matrix(nrow = ncol(modmat)^2, ncol = length(bayesmat)/ncol(modmat), data = 0)

  #Loop through each choice set
  for(i in altnames){

    #Select correct rows for current set and calculate probabilities
    altselect <- which(altvect == i)

    #Translate the exponentiated utilities to probabilities based on sum of exponentiated utilities in the current choice set
    exputilities[,altselect] <- exputilities[,altselect]/rowSums(exputilities[,altselect])

    #Translate to information matrix list
    #infoflat <- infoflat + apply(exputilities[,altselect], MARGIN = 1, function(x, Xdes = modmat[altselect,]) t(Xdes)%*%(diag(x) - x%*%t(x))%*%Xdes) #Replaced with faster crossprod function (though only makes sense to replace %*% with crossprod if a transpose is involved) below but accomplishes same calculation
    infoflat <- infoflat + apply(exputilities[,altselect], MARGIN = 1, function(x, Xdes = modmat[altselect,]) crossprod(Xdes,(diag(x) - tcrossprod(x))%*%Xdes))


  }

  # #Option 2:
  #
  # #Calculate the exponentiation of the utility for each alternative under entire matrix of parameter values in x
  # #(take transpose to make it easier to work with subsequent functions)
  # exputilities <- exp(modmat%*%x)
  #
  # #All in one loop
  # for(i in altnames){
  #
  #   #Select correct rows for current set and calculate probabilities
  #   altselect <- which(altvect == i)
  #
  #   #Create info matrix for each
  #   apply(exputilities[altselect,], MARGIN = 1, function(x) )
  #
  # }

  # #Calculate mahalanobis distance between x and mean for each column of x
  # distval <- stats::mahalanobis(t(x), center = parammean, cov = paramsigma)
  #
  # #Calculate the probability of each column of x occuring (based on Gaussian integrand formula)
  # probdens <- exp(-(length(parammean) * log(2 * pi) + logdet + distval)/2)
  # #
  # # #Note, lines below would be used to calculate if a stratified sample of x values is provided
  # # #Calculate weighted information matrix by multiplying info mat for each possible matchup with chance of that matchup happening and convert to matrix format
  # # info_mat <- matrix(infoflat%*%t(probdens), nrow = length(parammean))
  #
  # #Return the information matrix for each value of x multiplied by the density at x - if using cubature package but is very slow so no need to use anymore
  # #return(infoflat%*%diag(probdens))
  # return(matrix(data = apply(scale(infoflat, scale = 1/probdens, center = FALSE), MARGIN = 2, function(x) det(matrix(data = x, nrow = sqrt(length(x))))), nrow = 1))

  ##Note, lines below would be used to calculate if a stratified sample of x values is provided
  ##Calculate weighted information matrix by multiplying info mat for each possible matchup with chance of that matchup happening and convert to matrix format
  #info_mat <- matrix(infoflat%*%t(probdens), nrow = length(parammean))

  ##Note, line below would be used to calculate if a random normal vector of x values is drawn from an appropriate multivariate normal distribution
  ##Calculate weighted information matrix by multiplying info mat for each possible matchup with chance of that matchup happening and convert to matrix format
  info_mat <- matrix(rowSums(infoflat), nrow = ncol(modmat))/(length(bayesmat)/ncol(modmat))

}
