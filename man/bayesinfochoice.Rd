% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Defficiencyfuncs.R
\name{bayesinfochoice}
\alias{bayesinfochoice}
\title{Calculate Bayesian information matrix for discrete choice designs}
\usage{
bayesinfochoice(bayesmat, modmat, altvect)
}
\arguments{
\item{bayesmat}{matrix - matrix of points at which function should be returned. Each row is the evaluation point for that respective variable. Each column will be interpreted as a separate point.}

\item{modmat}{matrix - model matrix of single choice set to evaluate (produced using model.matrix(formula, design matrix))}

\item{altvect}{vector with integer corresponding to each row that indicates which choice set it is a member of}
}
\value{
The Bayesian information matrix of the design
}
\description{
Function to get bayesian information matrix for multinomial logit model based on provision of a random sample of points drawn from the prior parameter distribution
}
