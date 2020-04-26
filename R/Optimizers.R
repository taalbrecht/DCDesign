#' Optimizer update function to exchange point and recalculate objetive function
#'
#' @param x vector of numerics - new values for the row, column, or coordinate of design matrix being modified.
#'   Length must match destination length of target row, column, or coordinate of \code{designmat}
#' @param index integer - position of row or column being modified
#' @param roworcol character - whether to replace row, column, or coordinate.
#'   Uses "row", "column", or "coordinate" as argument choices
#' @param designmat whole design matrix in raw variable terms
#' @param base_input_range input type - 2 row matrix - matrix with variable minimums and maximums for base variables (variables common to all formulas in formulacollection). First row is minimum, second row is maximum. Names of colums should equal variable names used in formulas
#' @param formulacollection collection of formula objects used to construct scaled model matrices. Should include the following at a minimum for linear models. Additional arguments can be included for more complex operations
#' \describe{
#'   \item{inputranges}{matrix - A matrix with one column per variable in inputmat and minimum and maximum values in row 1 and 2, respectively}
#'   \item{baseformula}{formula - no-intercept formula of base algebraic combination variables}
#'   \item{range}{matrix - matrix with one column per variable in baseformula and minimum and maximum values in row 1 and 2, respectively}
#'   \item{termsalt}{vector, optional - names of alternate terms to replace terms in baseformula with so that the model matrix function will work. If not provided, returns scaled model matrix for baseformula}
#'   \item{retermedformula}{formula, optional - full model formula with terms replaced by altterms elements. If not provided, returns scaled model matrix for baseformula}
#'   }
#' @param scalemat logical, default TRUE - whether the base variable matrix for each model should be scaled from -1 to 1 or not
#' @param weight vector - weights to apply to each formula's objective function for overall objective function
#'
#' @return The objective function value for the updated design with the new row, column, or coordinate value(s)
#' @export
#'
#' @description Recalculates an objective function using the new row, column, or coordinate. Used for design optimization.
#'
#' @examples
optimizefn <- function(x, index, roworcol, designmat, base_input_range, formulacollection, scalemat = TRUE, weight = rep(1/length(formulacollection), times = length(formulacollection))){

  #Initialize with
  X <- designmat

  #Replace row based on searchdirection
  if(roworcol == "row"){

    X[index,] <- x

  }

  #Replace column based on searchdirection
  if(roworcol == "col"){

    X[,index] <- x

  }

  #Replace point based on searchdirection
  if(roworcol == "coordinate"){

    X[index] <- x

  }

  if(scalemat){
    #Loop through all formulas to get scaled design model matricies
    Xlist <- lapply(formulacollection, function(x) Xmat(inputmat = X, inputranges = base_input_range,
                                                        baseformula = x$baseformula, baseformrange = x$range,
                                                        altterms = x$termsalt, fullformulareterm = x$retermedformula))
  }else{
    #Loop through all formulas to get non-scaled design model matrices
    Xlist <- lapply(formulacollection, function(x) model.matrix(x$fullformula, data.frame(X)))

  }

  #Calculate objective function
  output <- sum(weight*mapply(function(x, y) objfunc(modmat = x, modparams = y), x = Xlist, y = formulacollection))

  #Return weighted penalty
  return(output)

}



#' Optimal design constructor for linear and multinomial logit discrete choice designs
#'
#' @param base_input_range 2 row matrix - Matrix with variable minimums and maximums for base variables. First row is minimum, second row is maximum. Names of colums should equal variable names used in formulas
#' @param formulalist list of formulas - list of formulas to simultaneously optimize design for
#' @param objectivelist list of characters - objective functions to use when optimizing each model in formulalist. Options are:
#' \describe{
#'   \item{D-Opt}{Maximizes X'X of centered and standardized model matrix}
#'   \item{I-Opt}{Minimizes the integral of the prediction variance over the design space}
#'   \item{SSE}{Minimizes sum of squared error between optimized design and values in startingmat. Uses centered and standardized model matrix for calculation.
#'              Useful for situations where a target starting design must be kept static with reference to a specific formula in formulalist that exists in a subspace of the entire design space.}
#'   }
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
#'   \item{DesignMatrix}{The resulting design matrix}
#'   \item{OptimizerOutput}{The objective function result for the resulting design matrix}
#'   \item{FixedObjects}{A number of properties of the designs that do not change as the design is optimized.
#'         Used for subsequent re-optimization when the the output of this function is recycled and re-optimized.}
#' @export
#'
#' @description Creates an experimental design simultaneously optimized for one or more models.
#'   Can choose between D or I optimal objective functions.
#'   Can be used to create linear or discrete choice designs.
#'   May also be used to create forward-looking discrete choice tournament designs.
#'
#' @examples
optimizemodellist <- function(base_input_range, formulalist, objectivelist = rep("D-Opt", times = length(formulalist)), startingmat = NULL, referencemat = NULL, npoints = NULL, searchdirection = "row", weight = rep(1/length(formulalist), times = length(formulalist)), randomstarts = FALSE, searchstrat = "numoptimize", candmat = NULL, cexpoints = NULL, constrainnodeslist = rep(NA, times = length(formulalist)), extraparams = NULL, scalemat = TRUE, verbose = FALSE, prevoptimobject = NULL){

  ##### Fixed object creation/loading

  #If previous model object was not provided, create "fixed" objects like formulas and moment matrices
  if(is.null(prevoptimobject)){
    #Extract terms only from each variable (not including interactions, etc. Used for centering and standardizing from -1 to 1 before building model matrix)
    baseformulalist <- lapply(formulalist, function(x) as.formula(paste("~", paste(paste(attributes(terms(x))$variables)[-1], collapse = "+"), "-1")))

    #Reassign all terms to different variable name to make sure geometric combinations are not deconstructed
    termslist <- lapply(formulalist, function(x) paste(attributes(terms(x))$variables)[-1])
    retermedformulalist <- lapply(formulalist, function(x) paste(x)[2])
    retermedbaseformulalist <- lapply(formulalist, function(x) paste("~", paste(paste(attributes(terms(x))$variables)[-1], collapse = "+"), "-1"))

    termsaltlist <- list()
    for(j in 1:length(formulalist)){

      #Define alternate term names
      termsaltlist[[j]] <- paste0("x", c(1:length(termslist[[j]])))

      #Replace original names with new term names
      for(i in 1:length(termslist[[j]])){

        retermedformulalist[[j]] <- gsub(pattern = termslist[[j]][i], replacement = paste0("x", i), x = retermedformulalist[[j]], fixed = TRUE)
        retermedbaseformulalist[[j]] <- gsub(pattern = termslist[[j]][i], replacement = paste0("x", i), x = retermedbaseformulalist[[j]], fixed = TRUE)

      }
      #Save retermed formulas in formula format
      retermedformulalist[[j]] <- as.formula(paste("~", retermedformulalist[[j]]))
      retermedbaseformulalist[[j]] <- as.formula(paste("~", retermedbaseformulalist[[j]]))

    }

    baseformulaiceptlist <- list()
    for(i in 1:length(baseformulalist)){

      #Reinitialize base variable formulas using re-termed formulas so they are correct and properly eliminate unwanted geometric combinations.
      baseformulalist[[i]] <- as.formula(paste("~", paste(termslist[[i]][termsaltlist[[i]] %in% paste(attributes(terms(retermedformulalist[[i]]))$variables)[-1]], collapse = "+"), "-1"))
      #Create base formulas with intercepts for algebraic rangefinding
      baseformulaiceptlist[[i]] <- as.formula(paste("~", paste(termslist[[i]][termsaltlist[[i]] %in% paste(attributes(terms(retermedformulalist[[i]]))$variables)[-1]], collapse = "+")))
    }

    #Get ranges of all base terms of each formula
    rangelist <- lapply(baseformulaiceptlist, function(x) algebraic_range(base_var_range = data.frame(base_input_range), algebraic_formula = x))

    #Get moment matrix for each input formula for use in I-Optimality
    mommatlist <- list()
    for(i in 1:length(baseformulalist)){

      #Store rangelist element i as temporary looprange to rename variables per retermed formula
      looprange <- rangelist[[i]]

      #Rename range columns to match retermed formula
      colnames(looprange) <- termsaltlist[[i]][1:ncol(looprange)]

      if(is.na(constrainnodeslist[[i]])){

        #Approximate constraining node set by sampling 1000 points from base variable space and including extreme points from base variable space if no constraining nodeset was provided to the function
        constrainnodesloop <- rbind(model.matrix(baseformulalist[[i]], data.frame(expand.grid(data.frame(base_input_range)))),
                                    model.matrix(baseformulalist[[i]], data.frame(apply(base_input_range, MARGIN = 2, function(x) runif(1000, min = min(x), max = max(x))))))

      }else{

        #Use provided constraining node set
        constrainnodesloop <- constrainnodeslist[[i]]

      }

      if(scalemat){
        #Scale constraining nodeset from -1 to 1 if scaling is indicated by scalemat argument
        constrainnodesloop <- scale(constrainnodesloop,
                                    center = apply(rangelist[[i]], MARGIN = 2, mean),
                                    scale = apply(rangelist[[i]], MARGIN = 2, function(x) x[2] - x[1])/2)

        #Scale the range matrix to -1 to 1 range
        looprange <- scale(looprange, center = apply(looprange, MARGIN = 2, mean),
                           scale = apply(looprange, MARGIN = 2, function(x) x[2] - x[1])/2)

      }

      #Get moment matrix for formula i subjected to the provided constraining nodes (if provided)
      mommatlist[[i]] <- designmommat(inputranges = looprange, modelformula = retermedformulalist[[i]],
                                      constrainnodes = constrainnodesloop)$MomentMatrix

    }

  }else{

    #If previous model object was provided, load "fixed" objects calculated above from previous model object to the function environment for use
    base_input_range <- prevoptimobject$FixedObjects$base_input_range
    formulalist <- prevoptimobject$FixedObjects$formulalist
    constrainnodeslist <- prevoptimobject$FixedObjects$constrainnodeslist
    baseformulaiceptlist <- prevoptimobject$FixedObjects$baseformulaiceptlist
    baseformulalist <- prevoptimobject$FixedObjects$baseformulalist
    mommatlist <- prevoptimobject$FixedObjects$mommatlist
    rangelist <- prevoptimobject$FixedObjects$rangelist
    retermedbaseformulalist <- prevoptimobject$FixedObjects$retermedbaseformulalist
    retermedformulalist <- prevoptimobject$FixedObjects$retermedformulalist
    termsaltlist <- prevoptimobject$FixedObjects$termsaltlist
    termslist <- prevoptimobject$FixedObjects$termslist

    ##Doesn't seem to work as one would think but maybe worth exploring later
    #list2env(prevoptimobject$FixedObjects)

  }

  ######End of fixed object creation/loading

  #Initialize staring design if none was provided
  if(is.null(startingmat)){
    startingmat <- apply(base_input_range, MARGIN = 2, function(x) runif(npoints, min = min(x), max = max(x)))
  }


  #Get model matrix in base variables from each input formula for SSE loss based on variable ranges calculated above
  if(!is.null(referencemat)){

    matlist <- list()
    for(i in 1:length(baseformulalist)){

      matlist[[i]] <- model.matrix(baseformulalist[[i]], data.frame(referencemat))

      #Scale from -1 to 1 if indicated by scalemat
      if(scalemat){
        matlist[[i]] <- scale(matlist[[i]],
                              center = apply(rangelist[[i]], MARGIN = 2, mean),
                              scale = apply(rangelist[[i]], MARGIN = 2, function(x) x[2] - x[1])/2)
      }
    }

    #If no reference matrix supplied, return NA for matlist elements
  }else{matlist <- rep(NA, times = length(formulalist))}

  #Initialize design matrix for looping for all variables
  testmat <- as.matrix(startingmat)

  #Scale from -1 to 1 if indicated by scalemat
  if(scalemat){

    testmat <- scale(testmat, center = apply(base_input_range, MARGIN = 2, mean),
                     scale = apply(base_input_range, MARGIN = 2, function(x) x[2] - x[1])/2)

  }

  #Initialize list of formula specific objects
  formulacollection <- list()
  for(i in 1:length(formulalist)){

    formulacollection[[i]] <- list("fullformula" = formulalist[[i]],
                                   "baseformula" = baseformulalist[[i]],
                                   "range" = rangelist[[i]],
                                   "termsalt" = termsaltlist[[i]],
                                   "retermedformula" = retermedformulalist[[i]],
                                   "momentmatrix" = mommatlist[[i]],
                                   "referencemat" = matlist[[i]],
                                   "objname" = objectivelist[[i]])

    #Replace reterming elements with NA for model matrix construction if SSE is to be used as objective function
    if(objectivelist[[i]] == "SSE"){

      formulacollection[[i]]$termsalt <- NULL
      formulacollection[[i]]$retermedformula <- NULL

    }

    #Add "extra elements" needed for non-linear model optimization as required if an extraparams object is provided
    if(length(extraparams[[i]]) > 0){

      formulacollection[[i]][names(extraparams[[i]])] <- extraparams[[i]]

    }

  }

  #Get current loss value from objective function
  bestout <- list()
  bestout$value <- optimizefn(x = testmat[1,], index = 1, roworcol = "row", designmat = testmat, base_input_range = base_input_range, formulacollection = formulacollection, weight = weight, scalemat = scalemat)

  #Numeric optimizer search strategy
  if(searchstrat == "numoptimize"){

    #Calculate minimum and maximum value matrices for optimizer
    minmat <- matrix(rep(apply(base_input_range, MARGIN = 2, min), nrow(testmat)),
                     ncol = ncol(base_input_range), byrow = TRUE)
    maxmat <- matrix(rep(apply(base_input_range, MARGIN = 2, max), nrow(testmat)),
                     ncol = ncol(base_input_range), byrow = TRUE)

    #Scale all minimum and maximum values from -1 to 1 if scalemat parameter is true
    if(scalemat){
      minmat[,] <- -1
      maxmat[,] <- 1

    }

    #Initialize loop count and reshape minimum/maximum boundary matrices appropriately
    if(searchdirection == "row"){loopcount <- nrow(testmat)}
    if(searchdirection == "column"){
      loopcount <- ncol(testmat)
      minmat <- t(minmat)
      maxmat <- t(maxmat)
    }
    if(searchdirection == "coordinate"){
      loopcount <- length(testmat)
      minmat <- matrix(minmat, ncol = 1)
      maxmat <- matrix(maxmat, ncol = 1)
    }

    for(j in 1:loopcount){

      if(randomstarts){

        #optimize with randomized candidate row/column with error handling returning Inf on error for solver errors
        loopout <- tryCatch(optim(apply(rbind(minmat[j,], maxmat[j,]), MARGIN = 2, function(x) runif(1, min = x[1], max = x[2])), # Starting values for optimizer
                                  #runif(length(testmat)/loopcount, min = -1, max = 1),
                                  #lower = rep(-1, length(testmat)/loopcount), upper = rep(1, length(testmat)/loopcount), #Lower and upper bounds
                                  fn = optimizefn, method = "L-BFGS-B", lower = minmat[j,], upper = maxmat[j,], #Optimizer arguments
                                  index = j, roworcol = searchdirection, designmat = testmat, base_input_range = base_input_range, formulacollection = formulacollection, weight = weight, scalemat = scalemat), #Additional static arguments passed to loss function
                            error = function(x) list(value = Inf))  #Error handling to prevent solver errors from crashing looping by returning Inf so results won't be used

      }else{

        #optimize without randomizing current point/row/column
        if(searchdirection == "row"){prevvals <- testmat[j,]}
        if(searchdirection == "column"){prevvals <- testmat[,j]}
        if(searchdirection == "coordinate"){prevvals <- testmat[j]}

        loopout <- tryCatch(optim(prevvals, #Starting values for optimizer
                                  fn = optimizefn, method = "L-BFGS-B", lower = minmat[j,], upper = maxmat[j,], #Optimizer arguments
                                  index = j, roworcol = searchdirection, designmat = testmat, base_input_range = base_input_range, formulacollection = formulacollection, weight = weight, scalemat = scalemat), #Additional static arguments passed to loss function
                            error = function(x) list(value = Inf))  #Error handling to prevent solver errors from crashing looping by returning Inf so results won't be used


      }

      #If loss function has been improved, keep changes. otherwise discard them
      if(loopout$value < bestout$value){
        if(searchdirection == "row"){testmat[j,] <- loopout$par}
        if(searchdirection == "column"){testmat[,j] <- loopout$par}
        if(searchdirection == "coordinate"){testmat[j] <- loopout$par}
        bestout <- loopout
        if(verbose){print(loopout$value)}
      }

    }

  }

  #Row exchange search strategy
  if(searchstrat == "fedorov"){

    #Get parameters used in loops
    candpoints <- nrow(candmat)


    #Pull starting design from starting matrix if one was provided and add it to the candidate set (though the points will not be retried in the candidate point search.
    #This is by design. Add the starting matrix to the candidate set before calling this function if those points should be searched.)
    if(!is.null(startingmat)){

      candmat <- rbind.data.frame(candmat, startingmat)
      bestout$par <- c((candpoints + 1):(candpoints + nrow(startingmat)))

    }else{

      #Randomly sample candidate matrix without replacement
      bestout$par <- sample(c(1:candpoints), size = npoints)

    }


    #Loop through all formulas to get design model matrices for entire candidate set
    if(scalemat){
      #Get scaled design model matrices if indicated by scalemat
      Xlist <- lapply(formulacollection, function(x) Xmat(inputmat =  scale(candmat, center = apply(base_input_range, MARGIN = 2, mean),
                                                                            scale = apply(base_input_range, MARGIN = 2, function(x) x[2] - x[1])/2),
                                                          inputranges = base_input_range,
                                                          baseformula = x$baseformula, baseformrange = x$range,
                                                          altterms = x$termsalt, fullformulareterm = x$retermedformula))
    }else{
      #Get unscaled desigm model matrices if indicated by scalemat
      Xlist <- lapply(formulacollection, function(x) model.matrix(x$fullformula, data.frame(candmat)))

    }


    #Initialize loop variables and calculate starting objective function
    loopout <- bestout
    bestout$value <- sum(weight*mapply(function(x, y) objfunc(modmat = x, modparams = y), x = lapply(Xlist, function(x) x[loopout$par,]), y = formulacollection))

    #Loop through each possible row of design
    for(i in 1:length(bestout$par)){

      #Reinitialize loopout to clear last point tried in previous index
      loopout <- bestout

      #Try replacing current row with every row in candidate matrix and see if there is an improvement
      for(j in 1:candpoints){

        loopout$par[i] <- j

        loopout$value <- sum(weight*mapply(function(x, y) objfunc(modmat = x, modparams = y), x = lapply(Xlist, function(x) x[loopout$par,]), y = formulacollection))

        #If design has been improved, save new point
        if(loopout$value < bestout$value){
          bestout <- loopout
          if(verbose){print(paste(loopout$value, "row", i, "candpoint", j))}
        }

      }

    }

  }


  #Coordinate exchange search strategy
  if(searchstrat == "cex"){

    #Initialize loopout with current design
    loopout <- bestout

    #Step through each point in the design matrix
    for(i in 1:length(testmat)){

      #Initialize step list to try
      if(scalemat){
        #Step from -1 to 1 if -1 to 1 scaling is indicated by scalemat
        pointvec <- seq(from = -1, to = 1, by = 2/(cexpoints-1))
      }else{

        #Step from minimum to maximum for current variable if no scaling should be applied
        pointvec <- seq(from = min(base_input_range[,ceiling(i/nrow(testmat))]),
                        to = max(base_input_range[,ceiling(i/nrow(testmat))]),
                        by = (max(base_input_range[,ceiling(i/nrow(testmat))]) - min(base_input_range[,ceiling(i/nrow(testmat))]))/(cexpoints-1))

      }

      #Try each step size as a replacement
      for(j in pointvec){
        #Get loss function
        loopout$value <- optimizefn(x = j, index = i, roworcol = "coordinate", designmat = testmat, base_input_range = base_input_range, formulacollection = formulacollection, weight = weight, scalemat = scalemat)

        #If design has been improved, save new point
        if(loopout$value < bestout$value){
          testmat[i] <- j
          bestout <- loopout
          if(verbose){print(paste(loopout$value, "element", i, "candpoint", j))}
        }


      }

    }

  }

  #Return final design matrix
  if(searchstrat %in% c("numoptimize", "cex")){

    #Rescale testmat to original variable scale if scaling was used for optimization
    if(scalemat){
      testmat <- data.frame(scale(
        scale(testmat, center = FALSE, scale = 2/apply(base_input_range, MARGIN = 2, function(x) x[2] - x[1])),
        center = -apply(base_input_range, MARGIN = 2, mean), scale = FALSE))
    }

    #Convert model matrix to a data frame
    testmat <- data.frame(testmat)

  }
  if(searchstrat == "fedorov"){

    testmat <- candmat[bestout$par,]

  }

  #Return list of outputs
  list("DesignMatrix" = testmat,
       "OptimizerOutput" = bestout,
       "FixedObjects" = list("base_input_range" = base_input_range, #From previous func call
                             "formulalist" = formulalist, #From previous func call
                             "constrainnodeslist" = constrainnodeslist, #From previous func call
                             "extraparams" = extraparams, #From previous func call
                             "baseformulaiceptlist" = baseformulaiceptlist, #Internal fixed object
                             "baseformulalist" = baseformulalist, #Internal fixed object
                             "mommatlist" = mommatlist, #Internal fixed object
                             "rangelist" = rangelist, #Internal fixed object
                             "retermedbaseformulalist" = retermedbaseformulalist, #Internal fixed object
                             "retermedformulalist" = retermedformulalist, #Internal fixed object
                             "termsaltlist" = termsaltlist, #Internal fixed object
                             "termslist" = termslist #Internal fixed object
       )
  )

}
