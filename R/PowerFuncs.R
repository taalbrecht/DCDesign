########################Power Calculation Method using variance covariance matrix######################################
#Calculates power to detect each supplied coefficient difference based on supplied alpha and beta error.

#Output: data frame with the following:
##first column containing parameter name
##second column containing difference to detect
##third column containing minimum sample size required to reach supplied power value
##fourth column containing alpha value used
##fifth column containing beta value used

#Inputs:
#vcovmat - variance-covariance matrix for model to estimate
#detectdiff - vector containing minmum detectable difference for each effect (column) of model matrix sized corresponding to standardized model matrix.
#test_alpha - numeric value or vector specifying allowable alpha error. If different alpha values are to be used for each estimated parameter, vector order must match
#test_beta - numeric value or vector specifying allowable beta error. If different beta values are to be used for each estimated parameter

vcovpower <- function(vcovmat, detectdiff, test_alpha = 0.05, test_beta = 0.2){

  #If a single value is supplied for test_alpha, copy that to a vector equal to the number of effects
  if(length(test_alpha) == 1){
    test_alpha <- rep(test_alpha, ncol(vcovmat))
  }

  #If a single value is supplied for beta, copy that to a vector equal to the number of effects
  if(length(test_beta) == 1){
    test_beta <- rep(test_beta, ncol(vcovmat))
  }

  #Transform alpha and beta power values to quantile functions
  z_one_minus_alpha<-qnorm(1-test_alpha)
  z_one_minus_beta<-qnorm(1-test_beta)

  #Calculate minimum sample size for each coefficient

#   # Use the parameter values as effect size. Other values can be used here.
#   effectsize<-paramestimates
#
#   # formula for sample size calculation is n>[(z_(beta)+z_(1-alpha))*sqrt(??????)/delta]^2
#   N<-((z_one_minus_beta + z_one_minus_alpha)*sqrt(diag(sigma_beta))/abs(effectsize))^2

  # formula for sample size calculation is n>[(z_(beta)+z_(1-alpha))*sqrt(??????)/delta]^2
  minsamplesize<-((z_one_minus_beta + z_one_minus_alpha)*sqrt(diag(vcovmat))/abs(detectdiff))^2

  output <- data.frame(rownames(vcovmat), detectdiff, minsamplesize, test_alpha, test_beta)

  colnames(output) <- c("Effect Name", "Difference to Detect", "Minimum Sample Size", "Alpha", "Beta")

  return(output)}


###########################################################################################
