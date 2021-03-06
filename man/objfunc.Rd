% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Defficiencyfuncs.R
\name{objfunc}
\alias{objfunc}
\title{Objective function calculator for linear and multinomial logit models for different objective functions}
\usage{
objfunc(modmat, modparams)
}
\arguments{
\item{modmat}{model matrix as produced by using the \code{\link[stats]{model.matrix}} function on the whole formula to be optimized - needed for all objective functions}

\item{modparams}{named list of all other elements needed for objective function below. Named elements are as follows:
\describe{
  \item{objname}{name of objective function to use - needed for all objective functions}
  \item{momentmatrix}{moment matrix. Must be a p x p matrix where p is the number of columns in modmat - only needed for I optimality}
  \item{referencemat}{reference matrix used to calculate squared error between current design and reference matrix - only needed for SSE objective function}
  \item{priorvals}{vector of prior values used for choice based (multinomial logit) optimality criteria}
  \item{altvect}{vector of question set IDs used for choice based (multinomial logit) optimality criteria}
  \item{matchupframe}{choice tournament matrix generated by function \code{\link{gentourneybracket}}. Used for calculating optimality for tournament based choice optimality criteria}
  }}
}
\value{
The value of the objective function
}
\description{
Universal function used to calculate a wide variety of objective functions for linear and multinomial logit models.
  Typically called during the model search process to calculate the objective function on every iteration.
}
