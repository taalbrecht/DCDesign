% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModelFrameProcessing.R
\name{algebraic_range}
\alias{algebraic_range}
\title{Find minimum and maximum for functions of base variables}
\usage{
algebraic_range(base_var_range, algebraic_formula)
}
\arguments{
\item{base_var_range}{matrix or data frame listing the range of base input variables used for algebraic model. Names must match those used in algebraic_formula
format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.}

\item{algebraic_formula}{a formula (one or two sided) that contains all of the algebraic combos to be evaluated. format = x~ A + B + A:B + I(A^2/B)+ I(ln(A)), etc}
}
\value{
matrix(AlgebraicRange), where Algebraic_Range is a matrix with column names matching algebraic combos on input side of formula. Row 1 = min and row 2 = max.
}
\description{
Function to return min and max range for functions of base variables. Calculations are based on provided range of base variables.
}
