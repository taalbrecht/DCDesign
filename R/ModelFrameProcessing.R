#' Change factor coding to contr.sum
#'
#' @param factorin vector of factor values
#' @param speclevels vector of factor level names in desired order
#'
#' @return vector as a factor with contrast sum encoding that respects the order in speclevels
#' @export
#'
#' @description Converts a factor contrast to contr.sum while preserving factor order and applying appropriate column names.
#'   Changes factor codings to contr.sum without changing global factor settings for R.
#'
#' @examples
customcontrsum <- function(factorin, speclevels){

      #Correct factor levels to proper order
  factorin <- factor(factorin, levels = speclevels)

      #Apply contrast sum encoding
      contrasts(factorin) <- contr.sum(nlevels(factorin))

      #Apply contrast column names
      colnames(contrasts(factorin)) <- rownames(contrasts(factorin))[1:ncol(contrasts(factorin))]

      return(factorin)

}


#################################Standardize Columns in Matrix Function##########################

#Calculates: Standardizes (from {-1 to 1}) or reverses standardization of list of column names in matrix

#Output: matrix(StandardMat), where StandardMat matches StartingMat but each column name in column_names has been standardized from -1 to 1 based on min and max values in StartingMat or Input_range if provided

#Inputs, required:
#StartingMat - a raw input matrix of continuous variables with column names matching variable names in input_formula. May also include outputs or variables to not be analyzed
#column_names - a list that contains all of the column names to be standardized. format = list("Response", "X1", ...)

#Inputs, optional:
#Input_range - matrix or data frame listing the range of columns used for standardization. Names must match those used in column_names
#			format is matrix or data frame with column name for each column to be standardized and minimum value in first row and maximum value in second row.
#reverse_standard - logical(TRUE, FALSE) indicating whether standardization should be reversed using Input_range. Default is FALSE

#NOTE: Replace this function with built in scale function after first draft of package is done
standardize_cols <- function(StartingMat, column_names, Input_range, reverse_standard = FALSE){

  #Initialize StandardMat by setting equal to StartingMat
  StandardMat <- StartingMat

  #For loop to standardize each variable column used in the formula from -1 to 1
  for(i in 1:length(column_names)){

    if(reverse_standard == FALSE){
      #If Input_range was not provided, pull min and max values from StartingMat. Otherwise, use min and max from Input_range
      if(missing(Input_range)){
        currentmax <- max(StartingMat[,column_names[i]])
        currentmin <- min(StartingMat[,column_names[i]])
        currentavg <- (currentmax + currentmin)/2
        currentrange <- currentmax-currentmin
      }
      else{currentmax <- max(Input_range[,column_names[i]])
      currentmin <- min(Input_range[,column_names[i]])
      currentavg <- (currentmax + currentmin)/2
      currentrange <- currentmax-currentmin
      }
      if(currentrange == 0){
        StandardMat[,column_names[i]] <- 0
      }
      else{
        StandardMat[,column_names[i]] <- 2*(StartingMat[,column_names[i]]-currentavg)/currentrange
      }
    }
    else{
      #Reverse standardization. If Input_range was not provided, return same matrix as StartingMat
      if(missing(Input_range)){}
      else{currentmax <- max(Input_range[,column_names[i]])
      currentmin <- min(Input_range[,column_names[i]])
      currentavg <- (currentmax + currentmin)/2
      currentrange <- currentmax-currentmin

      StandardMat[,column_names[i]] <- (currentrange/2)*(StartingMat[,column_names[i]])+currentavg
      }
    }
  }
  #Return Standardized Matrix
  return(StandardMat)}

##################################################################################################################

###############################Pull input variables from formula function#######################################

#Calculates: list of all unique fundamental input variables from formula passed to function

#Output: list(InputVarList), where:
#InputVarList - list of all input variables in function, format = list("x", "y", "z") where input formula may include combinations of x, y, and z

#Inputs:
#input_formula - a formula (one or two sided) that contains all of the input variables to be modeled. format = x~ A + B + A:B + I(A^2/B), etc

#NOTE: Replace this function with built in all.vars function after first draft of package is done
list_input_variables <- function(input_formula){

  #Pull inputs only from input_formula and count
  Input_Func <- as.formula(paste(nlme::splitFormula(input_formula)))
  InputVarList <- all.vars(Input_Func)

  # #Pull inputs only from input_formula and count - work on this
  # InputVarList <- attributes(terms(input_formula))$term.labels

  #Return list of input variables
  return(InputVarList)
}

#########################################################################################################################


#' Find minimum and maximum for functions of base variables
#'
#' @param base_var_range matrix or data frame listing the range of base input variables used for algebraic model. Names must match those used in algebraic_formula
#'   format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.
#' @param algebraic_formula a formula (one or two sided) that contains all of the algebraic combos to be evaluated. format = x~ A + B + A:B + I(A^2/B)+ I(ln(A)), etc
#'
#' @return matrix(AlgebraicRange), where Algebraic_Range is a matrix with column names matching algebraic combos on input side of formula. Row 1 = min and row 2 = max.
#' @export
#'
#' @description Function to return min and max range for functions of base variables. Calculations are based on provided range of base variables.
#'
#' @examples
algebraic_range <- function(base_var_range, algebraic_formula){
  #Calculates: Max and min for algebraic combination of base variables when range of base variables is known
  #INPUTS
  ##base_var_range - matrix or data frame listing the range of base input variables used for algebraic model. Names must match those used in algebraic_formula
  ##			format is matrix or data frame with column name for each base input variable then minimum value in first row and maximum value in second row.
  ##algebraic_formula - a formula (one or two sided) that contains all of the algebraic combos to be evaluated. format = x~ A + B + A:B + I(A^2/B)+ I(ln(A)), etc

  #OUTPUT: matrix(AlgebraicRange), with column names matching model matrix terms on input side of formula. Row 1 = min and row 2 = max.

  #Pull names of each algebraic term and factor combination on input side of formula
  algebraic_input_terms <- colnames(model.matrix(algebraic_formula, base_var_range), "factors")[-1]

  #Create matrix to store min and max values for each algebraic term
  AlgebraicRange <- matrix(data = 0, nrow = 2, ncol = length(algebraic_input_terms))
  colnames(AlgebraicRange) <- algebraic_input_terms
  rownames(AlgebraicRange) <- c("minimum", "maximum")

  #Loop to be executed for each algebraic term found

  for (i in 1:length(algebraic_input_terms)){

    #Replace : multiplication operators with * operators to ensure expression will work correctly
    transmult <- gsub(":","*",algebraic_input_terms[i])

    #Translate algebraic term to expression
    optimexp <- parse(text = transmult)

    #Pull base variables from algebraic term
    base_inputs <- all.vars(optimexp)

    #Create function using current algebraic term for use with optimization routine
    optimfunc <- function(optimvect){

      basevarmat <- matrix(nrow = 1, ncol = length(base_inputs), data = 0)
      basevarmat <- data.frame(basevarmat)
      colnames(basevarmat) <- base_inputs
      basevarmat[1,] <- optimvect
      output <- eval(optimexp,basevarmat)
      return(output)
    }

    #Create expanded grid using all possible combinations of base input variables to use as starting guesses

    trymat <- expand.grid(base_var_range[base_inputs])

    #Add a few randomized rows between minimum and maximum to avoid potential singularities
    #trymat <- rbind.data.frame(trymat, apply(base_var_range[base_inputs], MARGIN = 2, function(x) runif(5, min = x[1], max = x[2])))

    #Define vector of minimum and maximum values for base variables used in this algebraic term
    minvect <- as.vector(base_var_range[1,base_inputs])
    maxvect <- as.vector(base_var_range[2,base_inputs])

    #Find minimum and maximum values by trying all starting combinations in expanded grid
    AlgebraicRange["minimum",algebraic_input_terms[i]] <- min(unlist(lapply(apply(trymat, MARGIN = 1, optim, fn = optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect), "[", "value")))
    AlgebraicRange["maximum",algebraic_input_terms[i]] <- max(unlist(lapply(apply(trymat, MARGIN = 1, optim, fn = optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect, control = list(fnscale = -1)), "[", "value")))

    #     #Old, less accurate code that doesn't find min and max for factor combinations
    #     #Find minimum value for algebraic term
    #     AlgebraicRange["minimum",algebraic_input_terms[i]] <- optim(minvect, optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect)$value
    #     AlgebraicRange["maximum",algebraic_input_terms[i]] <- optim(maxvect, optimfunc, method = "L-BFGS-B", lower = minvect, upper = maxvect, control = list(fnscale = -1))$value
  }

  #Return Min and Max for Algebraic Functions
  return(AlgebraicRange)}

#################################################################################################################################


#' Fast flexible filling design
#'
#' @param inputranges matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
#' @param constrainnodes matrix with nodes of constraining manifold (hull) within rectangular inputranges space. Should be a matrix of nodes as produced by the convhulln function from the geometry package. It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
#' @param k integer - number of points to create for design
#' @param nrand integer, default = NULL - number of random points to generate per design point. If not supplied, 100 random points per design point will be used. If k0means algorithm doesn't converge, try adding more random starts
#' @param scaleinputs logical, default = TRUE - whether the input ranges
#' @param clustermethod character, default = "kmeans" - method used to identify clusters, can be ward (ward clustering) or kmeans (kmeans clustering)
#' @param centermethod character, default = "centroid" - method used to identify design point from each cluster. Options are:
#'   centroid - creates a point from each cluster at the centroid
#'   MaxPro - selects the points by choosing the point from each cluster that optimizes the maximum projection criteria
#'
#' @return
#'   \item{DesignMatrix}{A matrix of the resulting design with column names corresponding to variable names in inputranges}
#' @export
#'
#' @description Function to create fast flexible filling designs that uses clustering to identify design points that should be evenly spread across design space
#'
#' @examples
fastfill <- function(inputranges, constrainnodes = NULL, k, nrand = 100, scaleinputs = TRUE, clustermethod = "kmeans", centermethod = "centroid"){
  #Function for fast flexible filling designs that uses clustering to identify design points that should be evenly spread across design space
  #INPUTS:
  ##inputranges - matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
  ##constrainnodes - matrix with nodes of constraining manifold (hull) within rectangular inputranges space. Should be a matrix of nodes as produced by the convhulln geometry package. It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
  ##k - integer - number of points to create for design
  ##nrand - integer, default = NULL - number of random points to generate per design point. If not supplied, 100 random points per design point will be used. If k0means algorithm doesn't converge, try adding more random starts
  ##scaleinputs - logical, default = TRUE - whether the input ranges
  ##clustermethod - character, default = "kmeans" - method used to identify clusters, can be ward (ward clustering) or kmeans (kmeans clustering)
  ##centermethod - character, default = "centroid" - method used to identify design point from each cluster. Options are:
  ########centroid - creates a point from each cluster at the centroid
  ########MaxPro - selects the points by choosing the point from each cluster that optimizes the maximum projection criteria
  #OUTPUT:
  ##DesignMatrix - matrix with column names corresponding to variable names in inputranges

  #Generate random data either subject to rectangular constraints or a convex hull constraint

  #Check against constraining convex manifold if provided
  if(!is.null(constrainnodes)){

    randmat <- convexuniformfill(n = nrand*k, inputranges = inputranges, constrainnodes = constrainnodes)$UniformFill

    #CODE BELOW HAS BEEN REPLACED BY FUNCTION CALL ABOVE
    # #Initialize candidate point frame
    # randmat <- matrix(ncol = ncol(inputranges), nrow = 0)
    #
    # #Set placeholder count for
    # goodnodes <- c()
    #
    # #Set starting scale factor (percent chance a candidate point will be in the hull)
    # scalefactor <- 1.1*prod(apply(inputranges, MARGIN = 2, function(x) max(x) - min(x)))/convhulln(constrainnodes, options = "FA")$vol
    #
    # #Keep generating new data until enough points have been generated in the convex manifold
    # while(nrow(randmat) < nrand*k){
    #
    #   #Randomly generate data
    #   randmat <- rbind(randmat, apply(inputranges, MARGIN = 2, function(x) runif((nrand*k-nrow(randmat))*scalefactor, min = min(x), max = max(x))))
    #
    #   ##NOTE inhull function has horrible memory handling so need to loop through in 50000 point chunks. Will likely have to cap at a certain number of points and run through loop with garbage handling
    #   #Determine which randomly generated points are in or on the surface of the convex manifold
    #   outvec <- inhull(testpts = randmat, calpts = constrainnodes)
    #
    #
    #   #Convert to logical indicating TRUE if the point is in or on the surface of the convex manifold
    #   goodnodes <- outvec != -1
    #
    #   #Only keep points in manifold
    #   randmat <- randmat[goodnodes,]
    #   gc()
    #
    #   #Update scale factor (ratio of candidate points to points actually in manifold * safety factor)
    #   scalefactor <- 1.2*length(goodnodes)/sum(goodnodes)
    #
    # }
    #
    # #Remove extra valid points
    # randmat <- randmat[1:(nrand*k),]


  }else{

    #Randomly generate data subject to rectangular constraints
    randmat <- apply(inputranges, MARGIN = 2, function(x) runif(nrand*k, min = min(x), max = max(x)))

  }


  #Scale inputs before calculating clusters if desired
  if(scaleinputs){

    randmat <- scale(randmat, center = apply(inputranges, MARGIN = 2, mean),
                     scale = apply(inputranges, MARGIN = 2, function(x) x[2] - x[1])/2)

  }

  #Calculate clusters using ward clustering (note that this will not work for datasets that are too large)
  if(clustermethod == "ward"){

    #Get clusters using ward clustering and euclidean distance
    groupid <- cutree(hclust(d = dist(randmat, method = "euclidean"), method = "ward.D"), k = k)

    #Initialize design matrix to store results
    outmat <- randmat[1:k,]
    outmat[,] <- 0

    if(centermethod == "centroid"){
      #Extract centroids
      for(i in 1:k){
        outmat[i,] <- apply(randmat[groupid == i,], MARGIN = 2, mean)
      }

    }

    if(centermethod == "MaxPro"){

      outmat <- candsetmaxpro(randmat, groupid = groupid)$DesignMat


    }

  }

  #Calculate clusters using k-means
  if(clustermethod == "kmeans"){

    #Apply kmeans clustering algorithm
    groupid <- kmeans(randmat, centers = k, iter.max = 100)

    if(centermethod == "centroid"){
      #Extract centroids
      outmat <- groupid$centers

    }

    if(centermethod == "MaxPro"){

      #Identify maximum projection design from candidate set
      outmat <- candsetmaxpro(randmat, groupid = groupid$cluster)$DesignMat

    }

  }

  #Do not cluster before calculating designs
  if(clustermethod == "None"){

    #Apply kmeans clustering algorithm
    groupid <- kmeans(randmat, centers = k, iter.max = 100)

    if(centermethod == "centroid"){
      #Extract centroids
      outmat <- "Cannot use centroid method with no clusters"

    }

    if(centermethod == "MaxPro"){

      #Identify maximum projection design from candidate set
      outmat <- candsetmaxpro(randmat, npoints = k)$DesignMat

    }

  }

  #Reverse scaling operation if it was applied
  if(scaleinputs){

    outmat <- scale(scale(data.frame(outmat), center = FALSE,
                          scale = 2/apply(inputranges, MARGIN = 2, function(x) x[2] - x[1])), center = -apply(inputranges, MARGIN = 2, mean), scale = FALSE)

  }

  return(list("DesignMatrix" = outmat))

}


#' Check if candidate points are inside, outside, or on the surface of a n-dimensional convex manifold (hull)
#'
#' @param testpts An n x p matrix of points to test for hull membership. n data points, in p dimensions.
#'   If you have many points to test, it is most efficient to call this function once with the entire set.
#' @param calpts An m x p array of vertices of the convex hull, as produced by the geometry::convhulln function.
#' @param hull (OPTIONAL) tessellation (or triangulation) generated by convhulln.
#'   If hull is left empty or not supplied, then it will be generated which will take some extra time but will provide the same results.
#' @param tol Tolerance on the tests for inclusion in the convex hull.
#'   You can think of tol as the distance a point may possibly lie outside the hull, and still be perceived as on the surface of the hull.
#'   Because of numerical slop nothing can ever be done exactly here. I might guess a semi-intelligent value of tol to be tol = 1.e-13*mean(abs(calpts(:)))
#'   In higher dimensions, the numerical issues of floating point arithmetic will probably suggest a larger value of tol.
#'
#' @return A vector of integers of length n where each entry corresponds to the same row in testpts. Value meanings are:
#'   1 = inside hull
#'   -1 = outside hull
#'   0 = on hull (to precision indicated by tol)
#'
#' @export
#'
#' @description Use to determine if candidate points lie within or on the surface of a provided convex manifold.
#'   A similar function may eventually be added to the geometry package (currently developed and being added to the github version, just waiting for publication to CRAN).
#'   Based on Matlab code by John D'Errico 04 Mar 2006 (Updated 30 Oct 2006) \url{http://www.mathworks.com/matlabcentral/fileexchange/10226-inhull} with some modifications for greatly improved speed.
#'
#' @examples
inhull <- function(testpts, calpts, hull=geometry::convhulln(calpts), tol=mean(mean(abs(as.matrix(calpts))))*sqrt(.Machine$double.eps)) {
  #++++++++++++++++++++
  # R implementation of the Matlab code by John D'Errico 04 Mar 2006 (Updated 30 Oct 2006)
  # downloaded from

  #http://www.mathworks.com/matlabcentral/fileexchange/10226-inhull
  # with some modifications and simplifications
  #
  # Efficient test for points inside a convex hull in n dimensions
  #
  #% arguments: (input)
  #% testpts - nxp array to test, n data points, in p dimensions
  #% If you have many points to test, it is most efficient to
  #% call this function once with the entire set.
  #%
  #% calpts - mxp array of vertices of the convex hull, as used by
  #% convhulln.
  #%
  #% hull - (OPTIONAL) tessellation (or triangulation) generated by convhulln
  #% If hull is left empty or not supplied, then it will be
  #% generated.
  #%
  #% tol - (OPTIONAL) tolerance on the tests for inclusion in the
  #% convex hull. You can think of tol as the distance a point
  #% may possibly lie outside the hull, and still be perceived
  #% as on the surface of the hull. Because of numerical slop
  #% nothing can ever be done exactly here. I might guess a
  #% semi-intelligent value of tol to be
  #%
  #% tol = 1.e-13*mean(abs(calpts(:)))
  #%
  #% In higher dimensions, the numerical issues of floating
  #% point arithmetic will probably suggest a larger value
  #% of tol.
  #%
  # In this R implementation default

  #tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)
  # DEFAULT: tol = 1e-6
  #
  # VALUE: Matlab returns a vector of TRUE (inside/on) or FALSE (outside)
  # This R implementation returns an integer vector of length n
  # 1 = inside hull
  # -1 = outside hull
  # 0 = on hull (to precision indicated by tol)
  #--------------------------------------------------------
  # require(geometry, quietly=TRUE) # for convhulln
  # require(MASS, quietly=TRUE) # for Null

  # ensure arguments are matrices (not data frames) and get sizes
  calpts <- as.matrix(calpts)
  testpts <- as.matrix(testpts)
  p <- dim(calpts)[2] # columns in calpts
  cx <- dim(testpts)[1] # rows in testpts
  nt <- dim(hull)[1] # number of simplexes in hull
  # find normal vectors to each simplex
  nrmls <- matrix(NA, nt, p) # predefine each nrml as NA, degenerate

  degenflag <- matrix(TRUE, nt, 1)
  for (i in 1:nt) {
    nullsp <- t(MASS::Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))

    if (dim(nullsp)[1] == 1) { nrmls[i,] <- nullsp

    degenflag[i] <- FALSE}}
  # Warn of degenerate faces, and remove corresponding normals
  if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")

  nrmls <- nrmls[!degenflag,]
  nt <- dim(nrmls)[1]
  # find center point in hull, and any (1st) point in the plane of each simplex

  center = apply(calpts, 2, mean)
  a <- calpts[hull[!degenflag,1],]
  # scale normal vectors to unit length and ensure pointing inwards
  nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
  dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
  nrmls <- nrmls*matrix(dp, nt, p)
  # if min across all faces of dot((x - a),nrml) is
  # +ve then x is inside hull
  # 0 then x is on hull
  # -ve then x is outside hull
  # Instead of dot((x - a),nrml) use dot(x,nrml) - dot(a, nrml)
  #CHECK AGAINST LINES BELOW, THEN COMMENT OUT
  #aN <- diag(a %*% t(nrmls))
  #val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min)
  #Instead of creating gigantic matrix diag(a%*%t(nrmls)) then extracting diagonals as above, be smart instead and only calculate the diagonals using rowSums(a*nrmls)
  #MASSIVE IMPROVEMENT IN MEMORY HANDLING AND SPEED OVER METHOD USED ABOVE TO CALCULATE val and aN (useless quantities not calculated anymore, calculate one row at a time instead of whole matrix, both of which would cause out of memory issues)
  aN <- rowSums(a*nrmls)
  val <- apply(testpts, MARGIN = 1, FUN = function(x) min(x%*%t(nrmls) - aN))
  # code values inside 'tol' to zero, return sign as integer
  val[abs(val) < tol] <- 0
  as.integer(sign(val))
}

########################################################################
#Function to take a square range and a constraining set of hull points and fill the space with a uniform distribution of points

#' Fill a convex hull with a uniform random sample of points
#'
#' @param n integer - number of random points to generate
#' @param inputranges matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
#' @param constrainnodes matrix with nodes of constraining manifold (hull) within rectangular inputranges space. Should be a matrix of nodes as produced by the convhulln geometry package. It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
#'
#' @return
#'   \item{UniformFill}{Matrix containing uniformly randomly distributed points in the supplied manifold. Column names are transferred from inputranges}
#' @export
#'
#' @description This function uses a hyper-rectangular proposal region and a constraining set of hull points to fill the convex hull with a uniform distribution of points
#'
#' @examples
convexuniformfill <- function(n, inputranges, constrainnodes){
  #INPUTS:
  ##n - integer - number of random points to generate
  ##inputranges - matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
  ##constrainnodes - matrix with nodes of constraining manifold (hull) within rectangular inputranges space. Should be a matrix of nodes as produced by the convhulln geometry package. It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
  #OUTPUTS:
  #UniformFill - matrix containing uniformly randomly distributed points in the supplied manifold. Column names are transferred from inputranges

  #Initialize candidate point frame
  randmat <- matrix(ncol = ncol(inputranges), nrow = 0)

  #Set placeholder count for
  goodnodes <- c()

  #Set starting scale factor (percent chance a candidate point will be in the hull)
  scalefactor <- 1.1*prod(apply(inputranges, MARGIN = 2, function(x) max(x) - min(x)))/geometry::convhulln(constrainnodes, options = "FA")$vol

  #Keep generating new data until enough points have been generated in the convex manifold
  while(nrow(randmat) < n){

    #Randomly generate data
    randmat <- rbind(randmat, apply(inputranges, MARGIN = 2, function(x) runif((n-nrow(randmat))*scalefactor, min = min(x), max = max(x))))

    #Determine which randomly generated points are in or on the surface of the convex manifold
    outvec <- inhull(testpts = randmat, calpts = constrainnodes)

    #Convert to logical indicating TRUE if the point is in or on the surface of the convex manifold
    goodnodes <- outvec != -1

    #Only keep points in manifold
    randmat <- randmat[goodnodes,]

    #Update scale factor (ratio of candidate points to points actually in manifold * safety factor)
    scalefactor <- 1.2*length(goodnodes)/sum(goodnodes)

  }

  #Remove extra valid points
  randmat <- randmat[1:(n),]

  #Return uniformly distributed random fill
  return(list("UniformFill" = randmat))

}


#' Create a maximum projection (MaxPro) design
#'
#' @param candmat matrix - A matrix of candidate points in the design space. Ideally should be a random uniform fill of the design space constructed using \code{\link{convexuniformfill}}
#' @param npoints integer - The number of points to use in the MaxPro design
#' @param groupid vector - group membership ID for all points in candmat
#'
#' @return
#'   \item{DesignMat}{A matrix with the resulting MaxPro design}
#'   \item{ObjFunc}{The value of the objective function being minimized to find the MaxPro design}
#' @export
#'
#' @description Function to find maximum projection design from set of candidate points and cluster IDs. Either npoints or groupid must be provided.
#'   TODO: ADD IN LOSS FUNCITON THAT USES nloptim (same approach and gradient function as MaxPro) AND DISTANCE FROM HULL AS ADDITIVE LOSS WITH EXTRA WEIGHT ON DISTANCE FROM HULL TO CREATE NUMERICAL OPTIMZIER THAT WILL WORK FOR MAXPRO
#'   DOUBLE CHECK LOSS FUNCTION. ALSO SEEMS TO HAVE ISSUE WITH GETTING STUCK IN LOCAL MINIMA SO WILL NOT WORK VERY WELL WITHOUT CLUSTERING FIRST.
#'
#' @examples
candsetmaxpro <- function(candmat, npoints = NULL, groupid = NULL){

  if(is.null(npoints) & is.null(groupid)){

    message("Either the number of points (npoints) or a vector of group IDs (groupid) must be supplied")
    stop()

  }

  if(!is.null(groupid)){

    npoints <- length(unique(groupid))

  }


  #Scale input matrix from 0 to 1 across all factors
  D0 <- apply(candmat,2, function(x) (x-min(x))/(max(x)-min(x)))

  #Initialize choice vector and objective function by selecting first point from each cluster

  if(is.null(groupid)){

    designvec <- sample(c(1:nrow(D0)), npoints)

  }else{
    designvec <- sapply(c(1:npoints), function(x) sample(which(groupid == x), size = 1))
  }

  objfunc <- Inf

  #Loop over candidate points equal to number of points desired times (likely to achieve convergence)
  for(i in 1:npoints){

    for(j in 1:nrow(D0)){

      tempvec <- designvec

      #Replace point from candidate in cluster if clusters are used. Otherwise try all points
      if(is.null(groupid)){
        tempvec[i] <- j

      }else{
        tempvec[groupid[j]] <- j
      }

      #Calculate MaxPro optimization minimization criterion kernel (take log to nullify multiplication by combinatorial and exponentiation of 1/p)
      tempdis=(dist(D0[tempvec,1]))^2
      for(k in 2:ncol(D0)){
        tempdis=tempdis*(dist(D0[tempvec,k]))^2
      }

      #Invert resultant distance product and get sum to get maxpro criterion kernel
      tempobj <- log(sum(1/tempdis))

      #If the design has been improved, use the new design
      if(tempobj < objfunc){

        designvec <- tempvec
        objfunc <- tempobj

      }

    }

  }

  #Calculate full MaxPro criterion
  #Calculate MaxPro optimization minimization criterion kernel
  tempdis=(dist(D0[designvec,1]))^2
  for(k in 2:ncol(D0)){
    tempdis=tempdis*(dist(D0[designvec,k]))^2
  }
  #Change invert resultant distance product per maxpro criterion
  objfunc <- (sum(1/tempdis)/choose(n = npoints, k = 2))^(1/ncol(D0))

  return(list("DesignMat" = candmat[designvec,],
              "ObjFunc" = objfunc))

}


#' Calculate moment matrix of design space
#'
#' @param inputranges matrix - A matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
#' @param modelformula formula - The formula of the model to use for calculating moment matrix. Terms must match column names of inputranges otherwise model matrix cannot be calculated appropriately
#' @param constrainnodes matrix - A matrix with nodes of constraining manifold (hull) within rectangular inputranges space.
#'   Should be a matrix of nodes as produced by the as produced by \code{\link[geometry]{convhulln}}.
#'   It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
#' @param n integer - The number of points to use to numerically approximate moment matrix integral. It is suggested that at least 1000 points per factor in model formula be used
#'
#' @return
#'   \item{MomentMatrix}{The moment matrix}
#' @export
#'
#' @description Approximates the moment matrix of the design space calculated using numeric approximation of moment matrix based on uniform distribution of random points specified by the inputs.
#'   This moment matrix estimate is typically used to calculate I-Optimality for I-Optimal design creation
#'
#' @examples
designmommat <- function(inputranges, modelformula, constrainnodes = NA, n = 1000*ncol(attributes(terms(modelformula))$factors)){
  #INPUTS:
  ##inputranges - matrix with named columns where first row is minimum and second row is maximum - range of base input variables that will be used as ranges to generate random points for design space filling (including geometric combination variables)
  ##modelformula - formula of model to use for calculating moment matrix. Terms must match column names of inputranges otherwise model matrix cannot be calculated appropriately
  ##constrainnodes - matrix with nodes of constraining manifold (hull) within rectangular inputranges space. Should be a matrix of nodes as produced by the convhulln geometry package. It is ok to provide more points than necessary to construct the bounding manifold; they will simply be removed during the manifold construction process. Only works for convex manifolds for now
  ##n - integer, default = 1000 points per factor in modelformula - number of points to use to numerically approximate moment matrix integral. Also suggest at least 1000 points per factor in model formula
  #OUTPUT:
  ##MomentMatrix - numeric approximation of moment matrix based on uniform distribution of random points specified by inputs

  #Generate random data either subject to rectangular constraints or a convex hull constraint

  #Check against constraining convex manifold if provided
  if(!is.na(constrainnodes)){

    #Generate n points inside of constraining convex manifold, if provided
    randmat <- convexuniformfill(n = n, inputranges = inputranges, constrainnodes = constrainnodes)$UniformFill

  }else{

    #Randomly generate n datapoints subject to rectangular constraints of inputranges
    randmat <- apply(inputranges, MARGIN = 2, function(x) runif(n, min = min(x), max = max(x)))

  }

  #Calculate moment matrix
  fx <- model.matrix(modelformula, data.frame(randmat))

  #Calculate moment matrix
  MomentMatrix <- (t(fx)%*%fx)/n

  return(list("MomentMatrix" = MomentMatrix))
}


#' Calculates scaled model matrix
#'
#' @param inputmat matrix - design to evaluate in terms of variables in inputranges
#' @param inputranges matrix - matrix with one column per variable in inputmat and minimum and maximum values in row 1 and 2, respectively
#' @param baseformula formula - no-intercept formula of base algebraic combination variables
#' @param baseformrange  matrix - matrix with one column per variable in baseformula and minimum and maximum values in row 1 and 2, respectively
#' @param altterms vector, optional - names of alternate terms to replace terms in baseformula with so that the model matrix function will work. If not provided, returns scaled model matrix for baseformula
#' @param fullformulareterm formula, optional - full model formula with terms replaced by altterms elements. If not provided, returns scaled model matrix for baseformula
#'
#' @return Scaled model matrix for the design
#' @export
#'
#' @description Function to return properly scaled model matrix based off of input matrix of base variables.
#'   This function is compatible with constructed variables as well (e.g. I(A/B)).
#'
#' @examples
Xmat <- function(inputmat, inputranges, baseformula, baseformrange, altterms = NULL, fullformulareterm = NULL){
  #INPUTS:
  ##inputmat - matrix - design to evaluate in terms of variables in inputranges
  ##inputranges - matrix - matrix with one column per variable in inputmat and minimum and maximum values in row 1 and 2, respectively
  ##baseformula - formula - no-intercept formula of base algebraic combination variables
  ##baseformrange - matrix - matrix with one column per variable in baseformula and minimum and maximum values in row 1 and 2, respectively
  ##altterms - vector, optional - names of alternate terms to replace terms in baseformula with so that the model matrix function will work. If not provided, returns scaled model matrix for baseformula
  ##fullformulareterm - formula, optional - full model formula with terms replaced by altterms elements. If not provided, returns scaled model matrix for baseformula

  #Create model matrix: X
  #Translate to second formula space
  X <- model.matrix(baseformula, data.frame(scale(
    scale(inputmat, center = FALSE, scale = 2/apply(inputranges, MARGIN = 2, function(x) x[2] - x[1])),
    center = -apply(inputranges, MARGIN = 2, mean), scale = FALSE)))

  #Scale to -1 to 1 range for all variables based on input range
  X <- scale(X, center = apply(baseformrange, MARGIN = 2, mean),
             scale = apply(baseformrange, MARGIN = 2, function(x) x[2] - x[1])/2)

  #If larger retermed formula was also provided, return that model matrix. Otherwise, return scaled baseformula model matrix as calculated above
  if(!(is.null(altterms) | (is.null(fullformulareterm)))){

    #Rename columns with alternate designations to make sure geometric combinations are handled
    colnames(X) <- altterms[1:ncol(X)]

    #Create full model matrix - need to update by substituting variables. Won't work right now
    X <- model.matrix(fullformulareterm, data.frame(X))

  }

  return(X)

}
