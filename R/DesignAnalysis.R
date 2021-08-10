#' Create a surface plot of prediction variance against two variables
#'
#' @param covmat The covariance matrix for the design
#' @param model_formula The formula for the model (must match formula used to create covariance matrix)
#' @param xvar The variable name to show on the x-axis of the plot. Must match base term in model_formula
#' @param yvar The variable name to show on the y-axis of the plot. Must match base term in model_formula
#' @param xlims A vector giving the limits of the plot for xvar
#' @param ylims A vector giving the limits of the plot for yvar
#' @param steps An integer specifying the number of grid-points along each axis to use to create the plot
#' @param prevplot A previous plotly plot object if this surface should be added to an existing surface created by this function.
#'
#' @return A plotly surface plot of the prediction variance against xvar and yvar.
#' @export
#'
#' @description This function creates a surface plot of prediction variance against two predictor variables.
#'   Requires the plotly library to function.
#'
#' @examples
plotlycovmat <- function(covmat, model_formula, xvar, yvar, xlims = c(-1,1), ylims = c(-1,1), steps = 100, prevplot = NULL){

  plotmat <- expand.grid(x = seq(from = min(xlims), to = max(xlims), by = (max(xlims) - min(xlims))/steps),
                         y = seq(from = min(ylims), to = max(ylims), by = (max(ylims) - min(ylims))/steps))

  colnames(plotmat)[1] <- xvar
  colnames(plotmat)[2] <- yvar

  plotmatmodmat <- stats::model.matrix(model_formula, data.frame(plotmat))
  varest <- apply(plotmatmodmat,MARGIN = 1, function(x) t(x)%*%covmat%*%(x))

  zmat <- matrix(varest, nrow = (steps+1), byrow = TRUE)

  #Create new blank plotly object if one not supplied
  if(is.null(prevplot)){
    prevplot <- plotly::plot_ly()
  }

  # Add surface to plotly object
  prevplot <- plotly::add_surface(prevplot, z = zmat, x = seq(from = min(xlims), to = max(xlims), by = (max(xlims) - min(xlims))/steps),
                                  y = seq(from = min(ylims), to = max(ylims), by = (max(ylims) - min(ylims))/steps))

  return(prevplot)

}


#' Assess discrete chioce design metrics over a range of true parameter values
#'
#' @param designlist A list of discrete choice designs created by \code{\link{optimizemodellist}}
#' @param typevec A vector of strings describing the discrete choice design types. Options are:
#' \describe{
#'   \item{"Regular"}{A standard discrete choice design without a tournament.}
#'   \item{"Tournament"}{A fixed bracket choice tournament design where the same bracket assignments are used for each tournament replicate. (Includes forward-looking designs)}
#'   \item{"RandTournament"}{A random bracket choice tournament design where bracket assignments are randomized after each tournament completion.}
#'   }
#' @param trueframe The model matrix of true parameter values to iterate across where each row is a combination of true parameter values that should have a result returned.
#' @param design_params The vector of assumed parameter values to calculate distance from for each plot. Typically the vector of assumed parameter means or modes used to create the designs being assessed.
#' @param refindex Integer giving reference design position in designlist if design efficiencies should be returned with respec to that particular design.
#' @param linetype Type of line to use for each plot (only used if returntype is "plotly")
#' @param prevplots The return object created by a previous call to ploteffs if additional results should be added to that object.
#' @param returntype A string specifying the type of results to return. Options are:
#' \describe{
#'   \item{"plotly"}{Returns a list of plotly graph objects where each design in designlist will have a corresponding line on the plot.}
#'   \item{Any other argument}{Returns a list of numeric results where each design in designlist will have a corresponding entry in the list for each metric type.}
#'   }
#' @param loess_smooth_plots If TRUE, applies loess smoothing to efficiency results
#'
#' @return Returns a list of plots or a list of results with the following entries
#'   \item{x}{The x-values used for the plot}
#'   \item{seB1}{The standard error of the first coefficient in the design.}
#'   \item{DOpt}{The D-Optimality or D-Efficiency (if refindex is provided) of each design}
#'   \item{IOpt}{The I-Optimality or I-Efficiency (if refindex is provided) of each design}
#'   \item{OpCR}{The size of a confidence region centered at the true optimum in the design space that is 1 standard error in radius or the confidence region size relative to the reference design confidence region size if refindex is provided}
#'   \item{AOpt}{The A-Optimality or A-Efficiency (if refindex is provided) of each design}
#'   \item{ProbVar}{The average probability variance of each choice set or the ratio of the average probability variance to the reference design average probability variance if refindex is provided}
#'   \item{OptVar}{The prediction variance at the true optimum in the design space or the prediction variance relative to the reference design prediction variance if refindex is provided}
#' @export
#'
#' @description Plotting function for model efficiencies across a range of parameters (hard-coded) with optional overlay.
#'
#' @examples
ploteffs <- function(designlist, typevec, trueframe, design_params, refindex = NULL, linetype = rep("solid", length(designlist)), prevplots = NULL, returntype = "plotly", loess_smooth_plots = FALSE){

  # Define plot x axis as the euclidean distance from design_params to each point in trueframe (use rounding due to calculation imprecision)
  x_plot_vals = round(proxy::dist(trueframe, matrix(design_params, nrow = 1))[,1], 14)

  # Initialize lists to hold y results
  y1 = list()
  y2 = list()
  y3 = list()
  y4 = list()
  y5 = list()
  y6 = list()
  y7 = list()
  y8 = list()
  y9 = list()

  if(is.null(prevplots)){
    #Create plots if no plots from previous call supplied
    #Plot 1:
    #Plot of standard error for one estimator vs true value of that estimator for all models
    p1 <- plotly::plot_ly()

    #Plot of D-Optimality vs true value of that estimator for all models
    p2 <- plotly::plot_ly()

    #Plot of I-Optimality vs true value of that estimator for all models
    p3 <- plotly::plot_ly()

    #Plot of size of 95% confidence region for optimal point in design space
    p4 <- plotly::plot_ly()

    #Plot of A-Optimality vs true value of that estimator for all models
    p5 <- plotly::plot_ly()

    #Plot of sum of probability variance across all questions vs true value of A for all models
    p6 <- plotly::plot_ly()

    #Plot of confidence interval size vs true value of A for all models
    p7 <- plotly::plot_ly()

    #Plot of CDF of D-Optimality across all values in trueframe
    p8 <- plotly::plot_ly()

    #Plot of CDF of I-Optimality across all values in trueframe
    p9 <- plotly::plot_ly()

    #Name plots as "...Optimality" if no reference index is used
    if(is.null(refindex)){

      p1 <- plotly::layout(p = p1, title = "SE for Beta_A Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "Standard Error for Beta_A"))
      p2 <- plotly::layout(p = p2, title = "Actual D-Optimality Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "D-Optimality"))
      p3 <- plotly::layout(p = p3, title = "Actual I-Optimality Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "I-Optimality"))
      p4 <- plotly::layout(p = p4, title = "Optimal Point Confidence Region Relative Size Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "Single Standard Error A-B Confidence Region Size"))
      p5 <- plotly::layout(p = p5, title = "Actual A-Optimality Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "A-Optimality"))
      p6 <- plotly::layout(p = p6, title = "Average Choice Set Probability Variance Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "Average Choice Set Var(p)"))
      p7 <- plotly::layout(p = p7, title = "Response Estimator Variance at Optimal Point Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = "Response Estimator Variance at Optimal Point"))
      p8 <- plotly::layout(p = p8, title = "CDF of D-Optimality Over Sampled Parameters", xaxis = list(title = "D-Optimality"), yaxis = list(title = "Fraction of Data"))
      p9 <- plotly::layout(p = p9, title = "CDF of I-Optimality Over Sampled Parameters", xaxis = list(title = "I-Optimality"), yaxis = list(title = "Fraction of Data"))

    }else{

      #Name plots relative to reference plot if a reference is used
      p1 <- plotly::layout(p = p1, title = "Beta_A Efficiency Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("Beta_A Standard Error Efficiency of", names(designlist)[refindex])))
      p2 <- plotly::layout(p = p2, title = "Actual D-Efficiency Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("D-Efficiency of", names(designlist)[refindex])))
      p3 <- plotly::layout(p = p3, title = "Actual I-Efficiency Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("I-Efficiency of", names(designlist)[refindex])))
      p4 <- plotly::layout(p = p4, title = "Optimal Point Confidence Region Efficiency Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("A-B Confidence Region Size Efficiency of", names(designlist)[refindex])))
      p5 <- plotly::layout(p = p5, title = "Actual A-Efficiency Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("A-Efficiency of", names(designlist)[refindex])))
      p6 <- plotly::layout(p = p6, title = "Average Choice Set Probability Variance Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("Ratio of Average Choice Set Var(p) to", names(designlist)[refindex])))
      p7 <- plotly::layout(p = p7, title = "Response Estimator Variance at Optimal Point Under True Model", xaxis = list(title = "Euclidean Dist from design_params"), yaxis = list(title = paste("Ratio of Estimator Variance At Optimal Point to", names(designlist)[refindex])))
      p8 <- plotly::layout(p = p8, title = "CDF of D-Efficiency Over Sampled Parameters", xaxis = list(title = "D-Efficiency"), yaxis = list(title = "Fraction of Data"))
      p9 <- plotly::layout(p = p9, title = "CDF of I-Efficiency Over Sampled Parameters", xaxis = list(title = "I-Efficiency"), yaxis = list(title = "Fraction of Data"))

    }

  }else{
    p1 <- prevplots[[1]]
    p2 <- prevplots[[2]]
    p3 <- prevplots[[3]]
    p4 <- prevplots[[4]]
    p5 <- prevplots[[5]]
    p6 <- prevplots[[6]]
    p7 <- prevplots[[7]]
    p8 <- prevplots[[8]]
    p9 <- prevplots[[9]]
  }

  #Set designs to plot. Remove reference design from list
  if(is.null(refindex)){
    plotlist <- designlist
    plottype <- typevec
  }else{
    plotlist <- designlist[-c(refindex)]
    plottype <- typevec[-c(refindex)]
    linetype <- linetype[-c(refindex)]
  }

  #Calculate design details for reference design:
  if(!is.null(refindex)){

    #get list of info matrices for fixed bracket tournament designs at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(typevec[refindex] == "Tournament"){
      ylistbasic <- apply(trueframe, MARGIN = 1, function(x){
        output <- getchoicemodinfoprobs(modmat = stats::model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]], data.frame(designlist[[refindex]]$DesignMatrix)),
                                        paramestimates = as.vector(x),
                                        matchupframe = designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe)

        #Divide information matrix by number of questions asked to normalize:
        output$info_mat <- output$info_mat/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

        #Find the optimal point based on true betas
        x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
        x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
        #Store optimal point as model matrix point for use with covariance matrix
        output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

        #Calculate gradient matrix for generation of confidence region size for A and B at optimal point
        #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
        denom <- 4*x[4]*x[5]-x[6]^2
        output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                   0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)

        #Calculate sum of probability variances weighted by likelihood of each tournament question occurring
        output$p_var <- 0
        for(j in unique(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid)){

          tempselect <- designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid == j

          output$p_var <- output$p_var + prod(output$winprobvect[tempselect]*output$occurprobvect[tempselect][1])
        }

        #Divide p_var sum by number of questions asked to normalize:
        output$p_var <- output$p_var/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

        return(output)})
    }

    #get list of info matrices for random bracket tournament designs using the average of 10 random bracket arrangements at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(typevec[refindex] == "RandTournament"){
      #get list of info matrices for random tournament designs (100 random arrangements used) at all possible values of true params and gradient matrix at optimal point that maximizes true function
      loopcount <- 10
      ylistbasic <- apply(trueframe, MARGIN = 1, function(x){

        interrimout <- list()
        for(k in 1:loopcount){
          #Get new randomization of initial tournament pairings
          neword <- sample(unique(designlist[[refindex]]$FixedObjects$extraparams[[1]]$altvect)*2)
          neword <- c(rbind(neword-1, neword))

          output <- getchoicemodinfoprobs(modmat = stats::model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]],
                                                                data.frame(designlist[[refindex]]$DesignMatrix[neword,])),
                                          paramestimates = as.vector(x),
                                          matchupframe = designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe)

          #Divide information matrix by number of questions asked to normalize:
          output$info_mat <- output$info_mat/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

          #Find the optimal point based on true betas
          x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
          x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
          #Store optimal point as model matrix point for use with covariance matrix
          output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

          #Calculate gradient matrix for generation of confidence region size for A and B at optimal point
          #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
          denom <- 4*x[4]*x[5]-x[6]^2
          output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                     0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)

          #Calculate sum of probability variances weighted by likelihood of each tournament question occurring
          output$p_var <- 0
          for(j in unique(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid)){

            tempselect <- designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid == j

            output$p_var <- output$p_var + prod(output$winprobvect[tempselect]*output$occurprobvect[tempselect][1])
          }

          #Divide p_var sum by number of questions asked to normalize:
          output$p_var <- output$p_var/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

          interrimout[[k]] <- output}

        output <- list()
        output$info_mat <- Reduce("+", lapply(interrimout, "[[", "info_mat"))/loopcount
        output$optpoint <- Reduce("+", lapply(interrimout, "[[", "optpoint"))/loopcount
        output$gradmat <- Reduce("+", lapply(interrimout, "[[", "gradmat"))/loopcount
        output$p_var <- Reduce("+", lapply(interrimout, "[[", "p_var"))/loopcount

        return(output)})
    }

    #get list of info matrices for non-tournament designs at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(typevec[refindex] == "Regular"){
      ylistbasic <- apply(trueframe, MARGIN = 1, function(x){
        output <- d_effchoice(CurrentMatrix = stats::model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]], data.frame(designlist[[refindex]]$DesignMatrix)),
                              altvect = designlist[[refindex]]$FixedObjects$extraparams[[1]]$altvect,
                              paramestimates = as.vector(x), returninfomat = TRUE)

        #Divide information matrix and variance calculated by number of questions asked to normalize:
        output$info_mat <- output$info_mat/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$altvect)
        output$p_var <- output$p_var/max(designlist[[refindex]]$FixedObjects$extraparams[[1]]$altvect)

        #Find the optimal point based on true betas
        x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
        x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
        #Store optimal point as model matrix point for use with covariance matrix
        output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

        #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
        denom <- 4*x[4]*x[5]-x[6]^2
        output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                   0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)
        return(output)})
    }
  }


  for(i in 1:length(plotlist)){

    #get list of info matrices for fixed tournament designs at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(plottype[i] == "Tournament"){
      ylistloop <- apply(trueframe, MARGIN = 1, function(x){
        output <- getchoicemodinfoprobs(modmat = stats::model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]], data.frame(plotlist[[i]]$DesignMatrix)),
                                        paramestimates = as.vector(x),
                                        matchupframe = plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe)

        #Divide information matrix by number of questions asked to normalize:
        output$info_mat <- output$info_mat/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

        #Find the optimal point based on true betas
        x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
        x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
        #Store optimal point as model matrix point for use with covariance matrix
        output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

        #Calculate gradient matrix for generation of confidence region size for A and B at optimal point
        #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
        denom <- 4*x[4]*x[5]-x[6]^2
        output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                   0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)

        #Calculate sum of probability variances weighted by likelihood of each tournament question occurring
        output$p_var <- 0
        for(j in unique(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid)){

          tempselect <- plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid == j

          output$p_var <- output$p_var + prod(output$winprobvect[tempselect]*output$occurprobvect[tempselect][1])
        }

        #Divide p_var sum by number of questions asked to normalize:
        output$p_var <- output$p_var/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

        return(output)})
    }

    #get list of info matrices for random bracket tournament designs using the average of 10 random bracket arrangements at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(plottype[i] == "RandTournament"){
      #get list of info matrices for random tournament designs (100 random arrangements used) at all possible values of true params and gradient matrix at optimal point that maximizes true function
      loopcount <- 10
      ylistloop <- apply(trueframe, MARGIN = 1, function(x){

        interrimout <- list()
        for(k in 1:loopcount){
          #Get new randomization of initial tournament pairings
          neword <- sample(unique(plotlist[[i]]$FixedObjects$extraparams[[1]]$altvect)*2)
          neword <- c(rbind(neword-1, neword))

          output <- getchoicemodinfoprobs(modmat = stats::model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]],
                                                                data.frame(plotlist[[i]]$DesignMatrix[neword,])),
                                          paramestimates = as.vector(x),
                                          matchupframe = plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe)

          #Divide information matrix by number of questions asked to normalize:
          output$info_mat <- output$info_mat/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

          #Find the optimal point based on true betas
          x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
          x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
          #Store optimal point as model matrix point for use with covariance matrix
          output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

          #Calculate gradient matrix for generation of confidence region size for A and B at optimal point
          #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
          denom <- 4*x[4]*x[5]-x[6]^2
          output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                     0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)

          #Calculate sum of probability variances weighted by likelihood of each tournament question occurring
          output$p_var <- 0
          for(j in unique(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid)){

            tempselect <- plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$uniquesetid == j

            output$p_var <- output$p_var + prod(output$winprobvect[tempselect]*output$occurprobvect[tempselect][1])
          }

          #Divide p_var sum by number of questions asked to normalize:
          output$p_var <- output$p_var/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$matchupframe$questionid)

          interrimout[[k]] <- output}

        output <- list()
        output$info_mat <- Reduce("+", lapply(interrimout, "[[", "info_mat"))/loopcount
        output$optpoint <- Reduce("+", lapply(interrimout, "[[", "optpoint"))/loopcount
        output$gradmat <- Reduce("+", lapply(interrimout, "[[", "gradmat"))/loopcount
        output$p_var <- Reduce("+", lapply(interrimout, "[[", "p_var"))/loopcount

        return(output)})
    }

    #get list of info matrices for non-tournament designs using the average of 10 random bracket arrangements at all possible values of true params and gradient matrix at optimal point that maximizes true function
    if(plottype[i] == "Regular"){
      ylistloop <- apply(trueframe, MARGIN = 1, function(x){
        output <- d_effchoice(CurrentMatrix = stats::model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]], data.frame(plotlist[[i]]$DesignMatrix)),
                              altvect = plotlist[[i]]$FixedObjects$extraparams[[1]]$altvect,
                              paramestimates = as.vector(x), returninfomat = TRUE)

        #Divide information matrix and variance calculated by number of questions asked to normalize:
        output$info_mat <- output$info_mat/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$altvect)
        output$p_var <- output$p_var/max(plotlist[[i]]$FixedObjects$extraparams[[1]]$altvect)

        #Find the optimal point based on true betas
        x2 <-as.numeric((x[3] - x[6]*x[2]/(2*x[4]))/(x[6]^2/(2*x[4]) - 2*x[5]))
        x1 <- as.numeric(-(x[2] + x[6]*x2)/(2*x[4]))
        #Store optimal point as model matrix point for use with covariance matrix
        output$optpoint <- c(1, x1, x2, x1^2, x2^2, x1*x2)

        #output$gradmat <- matrix(c(0,-1/(2*x[4]),0,x[2]/2,0,0,0,-1/(2*x[5]),0,x[3]/2), ncol = 2)
        denom <- 4*x[4]*x[5]-x[6]^2
        output$gradmat <- matrix(c(0, -2*x[5]/denom, x[6]/denom,(2*x[2]*x[5]-x[3]*x[6])*4*x[5]/denom^2, 2*x[6]*(x[2]*x[6]-2*x[3]*x[4])/denom^2,(x[3]*x[6]^2-4*x[2]*x[5]*x[6]+4*x[3]*x[4]*x[5])/denom^2,
                                   0, -2*x[4]/denom, x[6]/denom,(2*x[3]*x[4]-x[2]*x[6])*4*x[4]/denom^2, 2*x[6]*(x[3]*x[6]-2*x[2]*x[5])/denom^2,(x[2]*x[6]^2-4*x[3]*x[4]*x[6]+4*x[2]*x[5]*x[4])/denom^2), ncol = 2)
        return(output)})
    }

    ##Calculate SE of coefficient "A" normalized per number of questions asked and add to plot p1
    # p1 <- plotly::add_lines(p = p1, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) solve(x$info_mat)[2,2]))/
    #                   (sapply(ylistbasic, function(x) solve(x$info_mat)[2,2])))

    y <- sapply(ylistloop, function(x) solve(x$info_mat)[2,2])
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) solve(x$info_mat)[2,2]))}

    p1 <- plotly::add_lines(p = p1, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y1[[(length(y1) + 1)]] = y

    ##Calculate D-Optimality of design normalized per number of questions asked and add to plot p2
    y <- (sapply(ylistloop, function(x) (1/det(x$info_mat))^(1/ncol(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) (1/det(x$info_mat))^(1/ncol(x$info_mat))))}

    # Average multiple y-values at each plot x value
    plot_vals = list(x = sort(unique(x_plot_vals)))
    plot_vals$y = unlist(lapply(plot_vals$x, function(x) mean(y[x_plot_vals == x])))

    # Use loess smooth if specified
    if(loess_smooth_plots){
      plot_vals = stats::loess.smooth(plot_vals$x, plot_vals$y)
    }
    p2 <- plotly::add_lines(p = p2, name = names(plotlist)[i], x = plot_vals$x, y = plot_vals$y, line = list(dash = linetype[i]))

    y2[[(length(y2) + 1)]] = y

    # Add D-Optimality CDF to plot 8 (use 1000 evenly spaced values along the range of results)
    ecdf_fun = stats::ecdf(y)
    ecdf_x = seq(from = min(y), to = max(y), length.out = 1000)
    y = ecdf_fun(ecdf_x)
    p8 <- plotly::add_lines(p = p8, name = names(plotlist)[i], x = ecdf_x, y = y, line = list(dash = linetype[i]))
    y8[[(length(y8) + 1)]] = y


    #Calculate I-Optimality of design normalized per number of questions asked and add to plot p3
    y <- sapply(ylistloop, function(x) sum(diag(plotlist[[i]]$FixedObjects$mommatlist[[1]]%*%solve(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) sum(diag(designlist[[refindex]]$FixedObjects$mommatlist[[1]]%*%solve(x$info_mat)))))}


    # Average multiple y-values at each plot x value
    plot_vals = list(x = sort(unique(x_plot_vals)))
    plot_vals$y = unlist(lapply(plot_vals$x, function(x) mean(y[x_plot_vals == x])))

    # Use loess smooth if specified
    if(loess_smooth_plots){
      plot_vals = stats::loess.smooth(plot_vals$x, plot_vals$y)
    }
    p3 <- plotly::add_lines(p = p3, name = names(plotlist)[i], x = plot_vals$x, y = plot_vals$y, line = list(dash = linetype[i]))

    y3[[(length(y3) + 1)]] = y

    # Add I-Optimality CDF to plot 8 (use 1000 evenly spaced values along the range of results)
    ecdf_fun = stats::ecdf(y)
    ecdf_x = seq(from = min(y), to = max(y), length.out = 1000)
    y = ecdf_fun(ecdf_x)
    p9 <- plotly::add_lines(p = p9, name = names(plotlist)[i], x = ecdf_x, y = y, line = list(dash = linetype[i]))
    y9[[(length(y9) + 1)]] = y


    #Calculate determinant of covariance matrix of A and B at optimal point (proportional to confidence region volume at that point) normalized per number of questions asked
    # p4 <- plotly::add_lines(p = p4, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2)))/
    #                   (sapply(ylistbasic, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2))))
#
#     y <- sapply(ylistloop, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2))
#     if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2)))}
#
#     p4 <- plotly::add_lines(p = p4, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
#     y4[[(length(y4) + 1)]] = y


    #Calculate A-Optimality of design normalized per number of questions asked and add to plot p5
    # p5 <- plotly::add_lines(p = p5, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) sum(diag(solve(x$info_mat)))))/
    #                   (sapply(ylistbasic, function(x) sum(diag(solve(x$info_mat))))))

    y <- sapply(ylistloop, function(x) sum(diag(solve(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) sum(diag(solve(x$info_mat)))))}

    p5 <- plotly::add_lines(p = p5, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y5[[(length(y5) + 1)]] = y

    #Calculate average question variance across design (weighted by probability of each tournament matchup occurring)
    # p6 <- plotly::add_lines(p = p6, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistbasic, function(x) x$p_var))/
    #                   (sapply(ylistloop, function(x) x$p_var)))

    y <- sapply(ylistloop, function(x) x$p_var)
    if(!is.null(refindex)){y <- (sapply(ylistbasic, function(x) x$p_var))/y}

    p6 <- plotly::add_lines(p = p6, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y6[[(length(y6) + 1)]] = y


    # #Calculate fit variance at optimal point for relative size comparison
    # y <- sapply(ylistloop, function(x) as.numeric((t(x$optpoint)%*%solve(x$info_mat)%*%x$optpoint)))
    # if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) as.numeric((t(x$optpoint)%*%solve(x$info_mat)%*%x$optpoint))))}
    #
    # p7 <- plotly::add_lines(p = p7, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    # y7[[(length(y7) + 1)]] = y

  }


  #Return plots
  if(returntype == "plotly"){
    return(list(p1, p2, p3, p4, p5, p6, p7, p8, p9))
  }else{
    return(list(
      x = trueframe[,c("A")],
      seB1 = y1,
      DOpt = y2,
      IOpt = y3,
      OpCR = y4,
      AOpt = y5,
      ProbVar = y6,
      OptVar = y7,
      DOptCDF = y8,
      IOptCDF = y9
    ))
  }

}
