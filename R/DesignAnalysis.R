library(plotly)
library(mlogit)

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

  plotmatmodmat <- model.matrix(model_formula, data.frame(plotmat))
  varest <- apply(plotmatmodmat,MARGIN = 1, function(x) t(x)%*%covmat%*%(x))

  zmat <- matrix(varest, nrow = (steps+1), byrow = TRUE)

  #Create new plot if one not supplied
  if(is.null(prevplot)){
    prevplot <- plot_ly(z = zmat, x = seq(from = min(xlims), to = max(xlims), by = (max(xlims) - min(xlims))/steps),
                        y = seq(from = min(ylims), to = max(ylims), by = (max(ylims) - min(ylims))/steps)) %>% add_surface()

  }else{
    #Otherwise, add surface to existing plot
    prevplot <- add_surface(prevplot, z = zmat, x = seq(from = min(xlims), to = max(xlims), by = (max(xlims) - min(xlims))/steps),
                            y = seq(from = min(ylims), to = max(ylims), by = (max(ylims) - min(ylims))/steps))
  }

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
#' @param refindex Integer giving reference design position in designlist if design efficiencies should be returned with respec to that particular design.
#' @param linetype Type of line to use for each plot (only used if returntype is "plotly")
#' @param trueframe The matrix of true parameter values to iterate across where each line is a combination of true parameter values that should have a result returned.
#' @param prevplots The return object created by a previous call to ploteffs if additional results should be added to that object.
#' @param returntype A string specifying the type of results to return. Options are:
#' \describe{
#'   \item{"plotly"}{Returns a list of plotly graph objects where each design in designlist will have a corresponding line on the plot.}
#'   \item{Any other argument}{Returns a list of numeric results where each design in designlist will have a corresponding entry in the list for each metric type.}
#'   }
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
ploteffs <- function(designlist, typevec, refindex = NULL, linetype = rep("solid", length(designlist)), trueframe = NULL, prevplots = NULL, returntype = "plotly", design_params = NULL, loess_smooth_plots = FALSE){

  #designlist - list of designs to plot. Name will be name shown in plot
  #typevec - vector of type of design equal in length to designlist. Should be "Tournament" or "Regular"
  #refindex - integer representing which design should be used as a reference to calculate efficiency ratios. If not provided, will simply plot optimality of each model
  #linetype - vector of names for type of line plotting. Per plotly: Sets the dash style of lines. Set to a dash type string ("solid", "dot", "dash", "longdash", "dashdot", or "longdashdot") or a dash length list in px (eg "5px,10px,2px,2px").
  #design_params - if provided, will return plots vs euclidean distance from design parameters
  #returntype - character. If "plotly", returns plot objects. Otherwise, returns list of x and each y result

  #Create true prior frame if one was not provided:
  if(is.null(trueframe)){
    trueframe <- data.frame(icept = rep(1, times = 201),
                            A = rep(1, times = 201),
                            B = rep(1, times = 201),
                            A2 = rep(-2, times = 201),
                            B2 = rep(-1, times = 201),
                            AB = rep(1, times = 201))

    trueframe$A <- seq(from = -9, to = 11, by = 0.1)
    trueframe <- as.matrix(trueframe)
  }

  # Define plot x axis
  if(!is.null(design_params)){
      # If design parameters are provided, this is euclidean distance from design_params to each point in trueframe (use rounding due to calculation imprecision)
      x_plot_vals = round(proxy::dist(trueframe, matrix(design_params, nrow = 1))[,1], 14)
    }else{
      x_plot_vals = trueframe[,c("A")]
    }

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
    p1 <- plot_ly()

    #Plot of D-Optimality vs true value of that estimator for all models
    p2 <- plot_ly()

    #Plot of I-Optimality vs true value of that estimator for all models
    p3 <- plot_ly()

    #Plot of size of 95% confidence region for optimal point in design space
    p4 <- plot_ly()

    #Plot of A-Optimality vs true value of that estimator for all models
    p5 <- plot_ly()

    #Plot of sum of probability variance across all questions vs true value of A for all models
    p6 <- plot_ly()

    #Plot of confidence interval size vs true value of A for all models
    p7 <- plot_ly()

    #Plot of CDF of D-Optimality across all values in trueframe
    p8 <- plot_ly()

    #Plot of CDF of I-Optimality across all values in trueframe
    p9 <- plot_ly()

    #Name plots as "...Optimality" if no reference index is used
    if(is.null(refindex)){

      p1 <- layout(p = p1, title = "SE for Beta_A Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "Standard Error for Beta_A"))
      p2 <- layout(p = p2, title = "Actual D-Optimality Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "D-Optimality"))
      p3 <- layout(p = p3, title = "Actual I-Optimality Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "I-Optimality"))
      p4 <- layout(p = p4, title = "Optimal Point Confidence Region Relative Size Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "Single Standard Error A-B Confidence Region Size"))
      p5 <- layout(p = p5, title = "Actual A-Optimality Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "A-Optimality"))
      p6 <- layout(p = p6, title = "Average Choice Set Probability Variance Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "Average Choice Set Var(p)"))
      p7 <- layout(p = p7, title = "Response Estimator Variance at Optimal Point Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = "Response Estimator Variance at Optimal Point"))
      p8 <- layout(p = p8, title = "CDF of D-Optimality Over Sampled Parameters", xaxis = list(title = "D-Optimality"), yaxis = list(title = "Fraction of Data"))
      p9 <- layout(p = p9, title = "CDF of I-Optimality Over Sampled Parameters", xaxis = list(title = "I-Optimality"), yaxis = list(title = "Fraction of Data"))

    }else{

      #Name plots relative to reference plot if a reference is used
      p1 <- layout(p = p1, title = "Beta_A Efficiency Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("Beta_A Standard Error Efficiency of", names(designlist)[refindex])))
      p2 <- layout(p = p2, title = "Actual D-Efficiency Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("D-Efficiency of", names(designlist)[refindex])))
      p3 <- layout(p = p3, title = "Actual I-Efficiency Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("I-Efficiency of", names(designlist)[refindex])))
      p4 <- layout(p = p4, title = "Optimal Point Confidence Region Efficiency Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("A-B Confidence Region Size Efficiency of", names(designlist)[refindex])))
      p5 <- layout(p = p5, title = "Actual A-Efficiency Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("A-Efficiency of", names(designlist)[refindex])))
      p6 <- layout(p = p6, title = "Average Choice Set Probability Variance Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("Ratio of Average Choice Set Var(p) to", names(designlist)[refindex])))
      p7 <- layout(p = p7, title = "Response Estimator Variance at Optimal Point Under True Model", xaxis = list(title = "True Value of Coefficient A"), yaxis = list(title = paste("Ratio of Estimator Variance At Optimal Point to", names(designlist)[refindex])))
      p8 <- layout(p = p8, title = "CDF of D-Efficiency Over Sampled Parameters", xaxis = list(title = "D-Efficiency"), yaxis = list(title = "Fraction of Data"))
      p9 <- layout(p = p9, title = "CDF of I-Efficiency Over Sampled Parameters", xaxis = list(title = "I-Efficiency"), yaxis = list(title = "Fraction of Data"))

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
        output <- getchoicemodinfoprobs(modmat = model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]], data.frame(designlist[[refindex]]$DesignMatrix)),
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

          output <- getchoicemodinfoprobs(modmat = model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]],
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
        output <- d_effchoice(CurrentMatrix = model.matrix(designlist[[refindex]]$FixedObjects$formulalist[[1]], data.frame(designlist[[refindex]]$DesignMatrix)),
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
        output <- getchoicemodinfoprobs(modmat = model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]], data.frame(plotlist[[i]]$DesignMatrix)),
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

          output <- getchoicemodinfoprobs(modmat = model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]],
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
        output <- d_effchoice(CurrentMatrix = model.matrix(plotlist[[i]]$FixedObjects$formulalist[[1]], data.frame(plotlist[[i]]$DesignMatrix)),
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
    # p1 <- add_lines(p = p1, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) solve(x$info_mat)[2,2]))/
    #                   (sapply(ylistbasic, function(x) solve(x$info_mat)[2,2])))

    y <- sapply(ylistloop, function(x) solve(x$info_mat)[2,2])
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) solve(x$info_mat)[2,2]))}

    p1 <- add_lines(p = p1, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y1[[(length(y1) + 1)]] = y

    ##Calculate D-Optimality of design normalized per number of questions asked and add to plot p2
    y <- (sapply(ylistloop, function(x) (1/det(x$info_mat))^(1/ncol(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) (1/det(x$info_mat))^(1/ncol(x$info_mat))))}
    if(!loess_smooth_plots){
      p2 <- add_lines(p = p2, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    }else{
      smoothed_dat = loess.smooth(x_plot_vals, y)
      smoothed_dat = list(x = sort(unique(x_plot_vals)))
      smoothed_dat$y = unlist(lapply(smoothed_dat$x, function(x) mean(y[x_plot_vals == x])))

      p2 <- add_lines(p = p2, name = names(plotlist)[i], x = smoothed_dat$x, y = smoothed_dat$y, line = list(dash = linetype[i]))
    }
    y2[[(length(y2) + 1)]] = y

    # Add D-Optimality CDF to plot 8 (use 1000 evenly spaced values along the range of results)
    ecdf_fun = ecdf(y)
    ecdf_x = seq(from = min(y), to = max(y), length.out = 1000)
    y = ecdf_fun(ecdf_x)
    p8 <- add_lines(p = p8, name = names(plotlist)[i], x = ecdf_x, y = y, line = list(dash = linetype[i]))
    y8[[(length(y8) + 1)]] = y


    #Calculate I-Optimality of design normalized per number of questions asked and add to plot p3
    y <- sapply(ylistloop, function(x) sum(diag(plotlist[[i]]$FixedObjects$mommatlist[[1]]%*%solve(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) sum(diag(designlist[[refindex]]$FixedObjects$mommatlist[[1]]%*%solve(x$info_mat)))))}

    if(!loess_smooth_plots){
      p3 <- add_lines(p = p3, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    }else{
      # smoothed_dat = loess.smooth(x_plot_vals, y)
      # p3 <- add_lines(p = p3, name = names(plotlist)[i], x = smoothed_dat$x, y = smoothed_dat$y, line = list(dash = linetype[i]))
      smoothed_dat = loess.smooth(x_plot_vals, y)
      smoothed_dat = list(x = sort(unique(x_plot_vals)))
      smoothed_dat$y = unlist(lapply(smoothed_dat$x, function(x) mean(y[x_plot_vals == x])))
      print(smoothed_dat$x)
      print(dim(trueframe))
      print(length(y))

      p3 <- add_lines(p = p3, name = names(plotlist)[i], x = smoothed_dat$x, y = smoothed_dat$y, line = list(dash = linetype[i]))
    }
    y3[[(length(y3) + 1)]] = y

    # Add I-Optimality CDF to plot 8 (use 1000 evenly spaced values along the range of results)
    ecdf_fun = ecdf(y)
    ecdf_x = seq(from = min(y), to = max(y), length.out = 1000)
    y = ecdf_fun(ecdf_x)
    p9 <- add_lines(p = p9, name = names(plotlist)[i], x = ecdf_x, y = y, line = list(dash = linetype[i]))
    y9[[(length(y9) + 1)]] = y


    #Calculate determinant of covariance matrix of A and B at optimal point (proportional to confidence region volume at that point) normalized per number of questions asked
    # p4 <- add_lines(p = p4, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2)))/
    #                   (sapply(ylistbasic, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2))))

    y <- sapply(ylistloop, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) (det(t(x$gradmat)%*%solve(x$info_mat)%*%x$gradmat))^(1/2)))}

    p4 <- add_lines(p = p4, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y4[[(length(y4) + 1)]] = y


    #Calculate A-Optimality of design normalized per number of questions asked and add to plot p5
    # p5 <- add_lines(p = p5, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistloop, function(x) sum(diag(solve(x$info_mat)))))/
    #                   (sapply(ylistbasic, function(x) sum(diag(solve(x$info_mat))))))

    y <- sapply(ylistloop, function(x) sum(diag(solve(x$info_mat))))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) sum(diag(solve(x$info_mat)))))}

    p5 <- add_lines(p = p5, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y5[[(length(y5) + 1)]] = y

    #Calculate average question variance across design (weighted by probability of each tournament matchup occurring)
    # p6 <- add_lines(p = p6, name = names(plotlist)[i], x = trueframe[,c("A")],
    #                 y = (sapply(ylistbasic, function(x) x$p_var))/
    #                   (sapply(ylistloop, function(x) x$p_var)))

    y <- sapply(ylistloop, function(x) x$p_var)
    if(!is.null(refindex)){y <- (sapply(ylistbasic, function(x) x$p_var))/y}

    p6 <- add_lines(p = p6, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y6[[(length(y6) + 1)]] = y


    #Calculate fit variance at optimal point for relative size comparison
    y <- sapply(ylistloop, function(x) as.numeric((t(x$optpoint)%*%solve(x$info_mat)%*%x$optpoint)))
    if(!is.null(refindex)){y <- y/(sapply(ylistbasic, function(x) as.numeric((t(x$optpoint)%*%solve(x$info_mat)%*%x$optpoint))))}

    p7 <- add_lines(p = p7, name = names(plotlist)[i], x = x_plot_vals, y = y, line = list(dash = linetype[i]))
    y7[[(length(y7) + 1)]] = y

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


#' Simulate discrete choice results for a design using a specified true model
#'
#' @param testbase data frame - design matrix of the discrete choice design to use for simulation
#' @param trueform formula - formula used to construct model matrix for simulation
#' @param trueparams vector - vector of true coefficient values corresponding to the exact order of terms returned by a model.matrix call using trueform
#' @param altvect vector - vector of integers that assign each row in testbase to a choice set
#' @param numsims integer - number of repeats to use for simulation
#' @param fullrandomize logical - whether the order of profiles in testbase should be fully randomized before each simulation run
#'   (i.e. randomly construct choice sets by sampling without replacement from the rows of testbase for each simulation run)
#' @param tournament logical - whether a single elimination winners' tournament should be used on the results of each simulation.
#' @param randtournament logical - whether the tournament bracket used for the tournament should be randomized (even if the design is not randomized per fullrandomize)
#'
#' @return
#'   \item{simframe}{data frame containing coefficient estimates resulting from running a model against each simulation in numsims}
#'   \item{optframe}{data frame containing estimated optimal parameter value combination to maximize trueform based on estimated coefficients in simframe}
#' @export
#'
#' @description This function simulates running a discrete choice study using numsims repeats and returns the results of the analysis.
#'   It is compatible with conventional and tournament discrete choice designs.
#'
#' @examples
choicesim <- function(testbase, trueform, trueparams, altvect, numsims,
                      fullrandomize = FALSE, tournament = FALSE, randtournament = FALSE){
  #INPUTS:
  ##testbase - data frame - design matrix to use for simulation
  ##trueform - formula - formula used to construct model matrix for simulation
  ##trueparams - vector - vector of true coefficient values corresponding to the exact order of terms returned by a model.matrix call using trueform
  ##altvect - vector - vector of integers that assign each row in testbase to a question set
  ##numsims - integer - number of repeats to use for simulation
  ##fullrandomize - logical - whether the order of alternatives in testbase should be fully randomized before each simulation run
  ##tournament - logical - whether a single elimination winners' tournament should be used on the results of each simulation.
  ##randtournament - logical - whether the tournament bracket used for the tournament should be randomized (even if the design is not randomized per fullrandomize)
  #OUTPUT: named list with following attributes:
  ##simframe - data frame containing coefficient estimates resulting from running a model against each simulation in numsims
  ##optframe - data frame containing estimated optimal parameter value combination to maximize trueform based on estimated coefficients in simframe

  #Redefined below as inputs
  # testbase <- cexdesign[[length(cexdesign)]]$DesignMatrix
  # trueform <- cexdesign[[1]]$FixedObjects$formulalist[[1]]
  # #trueparams <- c(0.5, 0.5, 0.5, -1, -1)
  # trueparams <- c(-0.5, -2, 2, -3, -3)
  # altvect <- cexdesign[[1]]$FixedObjects$extraparams[[1]]$altvect
  # #Whether design should be fully randomized (random pairing of profiles)
  # fullrandomize <- FALSE
  # #Whether a winners tournament should be used
  # tournament <- FALSE
  # #Whether tournament brackets should be randomized
  # randtournament <- FALSE



  ##Automated prep work (run every time loop below is run)
  #Generate model matrix
  testdesign <- testbase
  modmat <- model.matrix(trueform, data = testdesign)
  #Recode intercept appropriately
  iceptcol <- grep("(Intercept)", colnames(modmat), fixed = TRUE)
  if(length(iceptcol) == 1){
    #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
    modmat[,iceptcol] <- rep(c(0, rep(1, length(altvect)/length(unique(altvect)) - 1)), times = length(unique(altvect)))
  }


  #Calculate utility and probability vector for each choice set based on true model
  utilvect <- exp(modmat%*%trueparams)
  modprobs <- c()
  for(i in unique(altvect)){


    #modprobs <- c(modprobs, exp(modmat[altvect == i,]%*%trueparams)/sum(exp(modmat[altvect == i,]%*%trueparams)))
    modprobs <- c(modprobs, utilvect[altvect == i]/sum(utilvect[altvect == i]))

  }

  #Restructure design frame for mlogit package estimation
  testdesign <- data.frame(testdesign)
  testdesign$Alternative <- rep(c(1:(length(altvect)/length(unique(altvect)))), times = length(unique(altvect)))
  testdesign$Chosen <- NA
  testdesign$UniqueID <- c(1:nrow(testdesign))

  #Perhaps already structured correctly?
  # logitframe <- mlogit.data(testdesign,
  #                           shape = "long",
  #                           varying = which(!(colnames(testdesign) %in% c("Chosen", "Alternative"))),
  #                           alt.var = "Alternative",
  #                           choice = "Chosen",
  #                           id.var = "UniqueID")

  ##Simulation loop for any given design (number of full runs through design)
  simframe <- matrix(data = NA, nrow = numsims, ncol = 2*length(trueparams))
  #optframe <- matrix(data = NA, nrow = 1000, ncol = ncol(testbase))
  optframe <- data.frame(matrix(data = NA, nrow = numsims, ncol = ncol(testbase)))
  colnames(optframe) <- colnames(testbase)
  optframe$opt <- NA
  optframe$var <- NA
  testdat <- data.frame()

  for(i in 1:nrow(simframe)){

    loopframe <- data.frame(testdesign)
    loopframe$Chosen <- NA

    #Completely randomize choice sets if indicated
    if(fullrandomize){
      randorder <- sample(1:nrow(loopframe), size = nrow(loopframe), replace = FALSE)

      #Randomize choice profiles
      loopframe <- loopframe[randorder,]
      loopframe$Alternative <- rep(c(1:(length(altvect)/length(unique(altvect)))), times = length(unique(altvect)))

      modmat <- model.matrix(trueform, data = loopframe)
      #Recode intercept appropriately
      iceptcol <- grep("(Intercept)", colnames(modmat), fixed = TRUE)
      if(length(iceptcol) == 1){
        #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
        modmat[,iceptcol] <- rep(c(0, rep(1, length(altvect)/length(unique(altvect)) - 1)), times = length(unique(altvect)))
      }


      #Recalculate set probabilities
      utilvect <- exp(modmat%*%trueparams)
      modprobs <- c()
      for(j in unique(altvect)){

        modprobs <- c(modprobs, utilvect[altvect == j]/sum(utilvect[altvect == j]))

      }

    }

    #Simulate results from each choice set
    for(j in unique(altvect)){

      loopframe$Chosen[altvect == j] <- as.logical(rmultinom(n = 1, size = 1, prob = modprobs[altvect == j]))

    }

    #Run tournament if applicable
    if(tournament){

      #Get the winning profiles
      winners <- loopframe$Chosen
      #winnerutils <- utilvect[winners]
      tournframe <- loopframe[winners,]
      #tournframe$Alternative <- rep(c(1,2), times = nrow(tournframe)/2)

      #Continue until a champion is crowned
      while(sum(winners) > 1){

        #Randomize tournament order if specified
        if(randtournament){
          tournframe <- tournframe[sample(c(1:nrow(tournframe)), size = nrow(tournframe)),]
        }


        #Set alternate ordering
        tournframe$Alternative <- rep(c(1,2), times = nrow(tournframe)/2)

        #Create model matrix
        tournmat <- model.matrix(trueform, data = tournframe)

        #Recode intercept appropriately
        iceptcol <- grep("(Intercept)", colnames(tournmat), fixed = TRUE)
        if(length(iceptcol) == 1){
          #If intercept is included in model, make it so it is interacted with all alternatives except for the first alternative
          tournmat[,iceptcol] <- rep(c(0, rep(1, length(altvect)/length(unique(altvect)) - 1)), times = nrow(tournmat)/2)
        }


        #Recalculate alternative utilities
        winnerutils <- exp(tournmat%*%trueparams)

        #Play each set of this tournament
        for(k in 1:(nrow(tournframe)/2)){

          tournframe$Chosen[c((k*2-1):(k*2))] <- as.logical(rmultinom(n = 1, size = 1, prob = winnerutils[c((k*2-1):(k*2))]/sum(winnerutils[c((k*2-1):(k*2))])))

        }

        #Add tournament results to loop frame
        loopframe <- rbind.data.frame(loopframe, tournframe)

        #Get the winning profiles for the next round
        winners <- tournframe$Chosen
        #winnerutils <- winnerutils[winners]
        tournframe <- tournframe[winners,]

      }

    }

    testdat <- rbind.data.frame(testdat, loopframe)
    #testdat$UniqueID <- c(1:nrow(testdat))

    #Structure data correctly for mlogit package (needs mlogit attribute even though data frame is probably correct)
    logitframe <- mlogit.data(testdat,
                              shape = "long",
                              varying = which(!(colnames(testdat) %in% c("Chosen", "Alternative"))),
                              alt.var = "Alternative",
                              choice = "Chosen")


    #Rerun model and extract coefficients
    tempmod <- tryCatch(mlogit(as.formula(paste0("Chosen", paste(trueform, collapse = ""))), logitframe), error = function(x) NULL)

    #Extract coefficient estimates and standard errors
    simframe[i,] <- tryCatch(c(tempmod$coefficients, sqrt(diag(vcov(tempmod)))), error = function(x) NA)

    #Find predicted optimal point for each
    loopoptim <- tryCatch(optim(par = c(0,0), lower = rangemat[1,], upper = rangemat[2,], method = "L-BFGS-B", fn = function(x){

      testbase[1,] <- x
      output<- -model.matrix(trueform, testbase[1,])%*%tempmod$coefficients

      return(output)
    }), error = function(x) list(par = rep(NA, times = ncol(testbase))))

    #Enter predicted maximum points into data frame
    optframe[i,1:length(loopoptim$par)] <- loopoptim$par

    #Get model prediction and prediction variance
    optframe$optpred[i] <- tryCatch(model.matrix(trueform, optframe[i,])%*%tempmod$coefficients, error = function(x) NA)
    optframe$varpred[i] <- tryCatch(model.matrix(trueform, data = optframe[i,])%*%vcov(tempmod)%*%t(model.matrix(trueform, data = optframe[i,])), error = function(x) NA)



  }

  #Get actual result for optimal prediction point using true parameters (use model.matrix.lm since it allows for the na.pass option)
  optframe$optactual <- model.matrix.lm(trueform, optframe, na.action = 'na.pass')%*%trueparams

  #Name columns
  colnames(simframe) <- c(names(tempmod$coefficients), paste0("SE_", names(tempmod$coefficients)))

  output <- list("optframe" = optframe,
                 "simframe" = simframe)
  return(output)
}
