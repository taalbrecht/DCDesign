#' Optimal design constructor for multinomial logit designs
#'
#' @param base_input_range 2 row matrix - Matrix with variable minimums and maximums for base variables. First row is minimum, second row is maximum. Names of colums should equal variable names used in formulas
#' @param formulalist list of formulas - list of formulas to simultaneously optimize design for
#' @param objectivelist list of characters - objective functions to use when optimizing each model in formulalist. Options are:
#' \describe{
#'   \item{D-Opt}{Maximizes X'X of centered and standardized model matrix}
#'   \item{I-Opt}{Minimizes the integral of the prediction variance over the design space}
#'   \item{SSE}{Minimizes sum of squared error between optimized design and values in startingmat. Uses centered and standardized model matrix for calculation.
#'              Useful for situations where a target starting design must be kept static with reference to a specific formula in formulalist that exists in a subspace of the entire design space.}
#'          }
#' @param startingmat matrix with named columns - Matrix to use as a design starting point. Also used as the reference for SSE calculations when that loss function is used
#' @param referencemat matrix with named columns - matrix used as the reference for SSE calculations when that loss function is used
#' @param npoints integer - number of design points. Only needed if startingmat is not provided
#' @param searchdirection character - Method used to perform optimization search. Options are:
#' \describe{
#'   \item{row}{Optimize by modifying one entire row at a time}
#'   \item{column}{Optimize by modifying one entire column at a time}
#'   \item{coordinate}{Optimize by modifying a single coordinate at a time}
#'   }
#' @param weight vector of numerics - weight for each formula's loss function to use to calculate overall loss from each loss function
#' @param randomstarts TRUE/FALSE - whether a random start should be used for each row/column iteration through the optimizer
#' @param searchstrat character - the type of search strategy to use. Options are:
#' \describe{
#'   \item{numoptimize}{numeric optimization. Does not require candidate points.}
#'   \item{fedorov}{Row exchange. Requires candmat argument to be provided}
#'   \item{cex}{Coordinate exchange. Requires cexpoints argument to be provided.}
#'   }
#' @param candmat matrix with column names corresponding to base variable names in base_input_range - matrix of candidate points used for fedorov optimization.
#' @param cexpoints integer stating how many evenly spaced points should be used across the allowed variable range for the coordinate exchange search strategy
#' @param constrainnodeslist list of data frames or matrices - list of constraining nodes for each formula if a rectangular space based on the base_input_range shouldn't be used. This is important for I-Optimality
#' @param extraparams list of extra parameters required for non-linear model optimization. May contain
#' \describe{
#'   \item{altvect}{vector of integers - vector listing the choice set that each row belongs to for discrete choice optimizers}
#'   \item{paramestimates}{vector of numerics - vector containing the prior estimates for each coefficient in the choice model. MUST BE IN THE SAME ORDER AS WOULD BE RETURNED BY THE MODEL MATRIX FUNCTION}
#'   }
#' @param scalemat logical, default TRUE - whether the base variable matrix for each model should be scaled from -1 to 1 or not
#' @param verbose logical - whether to print progress statements from optimizer each time an improvement occurs
#' @param prevoptimobject result of previous run of this function - carries fixed constants over from previous execution (formula processing, moment matrix, etc) to help speed execution. ONLY USE IF YOU ARE PASSING THE RESULT FROM A PREVIOUS RUN WITH EXACTLY THE SAME MODEL STRUCTURES AND RANGES
#'
#' @return
#' \describe{
#'   \item{DesignMatrix}{The resulting design matrix}
#'   \item{OptimizerOutput}{The objective function result for the resulting design matrix}
#'   \item{FixedObjects}{A number of properties of the designs that do not change as the design is optimized.
#'         Used for subsequent re-optimization when the the output of this function is recycled and re-optimized.}
#'  }
#' @export
#'
#' @description Creates an experimental design simultaneously optimized for one or more models.
#'   Can choose between D or I optimal objective functions.
#'   Can be used to create linear or discrete choice designs.
#'   May also be used to create forward-looking discrete choice tournament designs.
#'
#' @examples
design_mlogit <- function(profiles_per_set, rounds, base_input_range, formulalist,  first_round_choice_sets = profiles_per_set ^ rounds, tolerance = .00001, ...){

  # choiceform <- formula(~A+B+I(A^2) + I(B^2))
  # paramestimates <- c(0.5, 0.5, 0.5, -1, -1) # Informative prior
  choiceform <- formula(~A+B+A*B+I(A^2) + I(B^2))
  paramestimates <- c(1,1,1,-2, -1, 1) # Informative prior
  paramsigma <- diag(length(paramestimates))
  extraparams <- list(list(altvect = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,
                                       9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16),
                           matchupframe = gentourneybracket(altvect = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,
                                                                        9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16), ChoiceTournamentAlts = 2)$matchupframe,
                           paramestimates = paramestimates,
                           bayesmat = mvtnorm::rmvnorm(n = 1000, mean = paramestimates, sigma = paramsigma)))


  #Varying parameters (fixed parameters will simply be hard-coded into pi groups)
  rangemat <- matrix(c(-1, 1,
                       -1, 1),
                     nrow = 2)
  colnames(rangemat) <- c("A", "B")


  #####################Test design matrix to demonstrate use####################
  #Randomly initialize matrix to start
  startingmat <- apply(rangemat, MARGIN = 2, function(x) runif(32, min = min(x), max = max(x)))
  #Create simple candidate matrix to test row exchange
  candmatsimple <- expand.grid(data.frame(apply(rangemat, MARGIN = 2, function(x) seq(from = min(x), to = max(x), by = (max(x) - min(x))/2))))


  #First optimizer run to calculate fixed parameters
  optimfixedparams <- optimizemodellist(base_input_range = rangemat,
                                        formulalist = list(choiceform),
                                        objectivelist = list("D-OptChoice"),
                                        referencemat = startingmat, startingmat = startingmat,
                                        searchdirection = "coordinate", weight = c(1), randomstarts = FALSE,
                                        searchstrat = "cex", candmat = candmatsimple, extraparams = extraparams,
                                        npoints = 32, cexpoints = 11)

  #Randomize design matrix to start all methods on level playing field
  optimfixedparams$DesignMatrix <- apply(rangemat, MARGIN = 2, function(x) runif(nrow(optimfixedparams$DesignMatrix), min = min(x), max = max(x)))

  ################Designs with forward looking choice tournament and informative prior

  #I-Opt with random start numeric optimization
  currentopt <- Inf
  continuevar <- TRUE
  Icoordnumdesignrandtourney <- list()
  Icoordnumdesignrandtourney[[1]] <- optimfixedparams
  Icoordnumdesignrandtourney[[1]]$cumulativeruntime <- 0
  #Run coordinate exchange selection to convergence
  t1 <- Sys.time()
  while(continuevar){

    prevopt <- currentopt
    #Run with random starting design with cex optimier
    Icoordnumdesignrandtourney[[length(Icoordnumdesignrandtourney) + 1]] <- optimizemodellist(base_input_range = rangemat,
                                                                                              formulalist = list(choiceform),
                                                                                              objectivelist = list("I-OptChoiceTournament"), weight = c(1),
                                                                                              searchdirection = "coordinate", randomstarts = TRUE,
                                                                                              searchstrat = "numoptimize", candmat = candmatsimple, extraparams = extraparams,
                                                                                              startingmat = Icoordnumdesignrandtourney[[length(Icoordnumdesignrandtourney)]]$DesignMatrix,
                                                                                              npoints = 100, cexpoints = 41, prevoptimobject = optimfixedparams)

    #Store cumulative runtime
    Icoordnumdesignrandtourney[[length(Icoordnumdesignrandtourney)]]$cumulativeruntime <- Sys.time() - t1

    #Get current optimality
    currentopt <- Icoordnumdesignrandtourney[[length(Icoordnumdesignrandtourney)]]$OptimizerOutput$value

    #Once improvements get close to convergence, stop the loop
    if((abs((prevopt - currentopt)/currentopt)) < 1/10000){continuevar <- FALSE}

    prevopt <- currentopt

    print(paste("I-Opt =", currentopt, "Time = ", Sys.time() - t1))

  }


}
