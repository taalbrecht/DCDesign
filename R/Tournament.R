#' Create choice tournament bracket structure
#'
#' @param altvect vector - vector defining which row of the design belongs to each question set for the starting questions (not including tournament bracket that will be generated by this function)
#' @param ChoiceTournamentAlts integer, default 2 - Maximum number of alternatives that should be used for each tournament matchup
#'
#' @return
#'   \item{matchupframe}{Matrix that details all possible matchups based on tournament bracket structure.}
#'   \item{tournamentframe}{Matrix showing the tournament bracket structure.}
#' @export
#'
#' @description Creates a full tournament bracket from the starting matchups described in altvect with no matchup containing more than ChoiceTournamentAlts profiles.
#'   Tournament will always be complete (one winning profile remaining) but may include one or more choice sets with fewer than ChoiceTournamentAlts in a match if necessary.
#'
#' @examples
gentourneybracket <- function(altvect, ChoiceTournamentAlts = 2){

  #Initialize tournament matrix
  tournamentframe <- data.frame(sourcequestion =unique(altvect) , tournamentset = unique(altvect))

  #loop through to create tournament pairings
  processedrow <- 0

  while(processedrow < nrow(tournamentframe)){

    #Create one tournament set from source question winners based on ChoiceTournamentAlts set size or the number of source questions remaining
    tournamentframe$tournamentset[(1+processedrow):min(processedrow+ChoiceTournamentAlts, nrow(tournamentframe))] <- max(tournamentframe) + 1
    processedrow <- min(processedrow+ChoiceTournamentAlts, nrow(tournamentframe))

    #Add rows for winner if there is more than one champion left
    if(processedrow < nrow(tournamentframe)){

      tournamentframe[(nrow(tournamentframe) + 1), ] <- c(max(tournamentframe), 0)

    }


  }

  #Create matrix of all possible alternative matchups. Start with specified starting alternative matchups that will be asked with 100% certainty
  matchupframe <- data.frame(matrixrowid = c(1:length(altvect)), questionid = altvect)

  #Add column for a unique number to assign to each possible choice set
  matchupframe$uniquesetid <- matchupframe$questionid

  #Loop through tournament matchups to identify all possible alternative matchups based on tournament pairings
  for(i in unique(tournamentframe$tournamentset)){

    #Initialize alternativelist
    looplist <- list()

    #Loop through all alternatives in each source question for this tournament matchup
    for(j in unique(tournamentframe$sourcequestion[which(tournamentframe$tournamentset == i)])){

      looplist[[(length(looplist) + 1)]] <- unique(matchupframe$matrixrowid[which(matchupframe$questionid == j)])

    }

    #Create matrix of all possible matchups for the current tournament question
    tempframe <- expand.grid(looplist)

    #Enter all new items into matchupframe
    for(j in 1:nrow(tempframe)){

      addframe <- data.frame(matrixrowid = c(t(tempframe[j,])),
                             questionid = rep(i, times = ncol(tempframe)),
                             uniquesetid = rep(max(matchupframe$uniquesetid) + 1, times = ncol(tempframe)))

      matchupframe <- rbind.data.frame(matchupframe, addframe)

    }

  }

  #Initialize variable for source question in matchupframe (0 indicates a non-tournament question)
  matchupframe$sourcequestion <- 0

  #Loop through each tournament set
  for(i in unique(tournamentframe$tournamentset)){

    #Loop through each row tied to the current tournament set
    for(j in which(matchupframe$questionid == i)){

      #Find the corresponding source question for the current row of the matchupframe
      matchupframe$sourcequestion[j] <- unique(matchupframe$questionid[((matchupframe$matrixrowid == matchupframe$matrixrowid[j]) &
                                                                          (matchupframe$questionid %in% tournamentframe$sourcequestion[tournamentframe$tournamentset == i]))])

    }


  }


  output <- list("matchupframe" = matchupframe,
                 "tournamentframe" = tournamentframe)
  return(output)
}


#' Calculate performance metrics for choice tournament designs
#'
#' @param modmat matrix - Model matrix of design to evaluate (produced using model.matrix(formula, design matrix))
#' @param paramestimates vector - Coefficient estimates to use for calculating probabilities for the formula passed in modform.
#'   Position must match coefficient column positions in modmat.
#' @param matchupframe matrix - matrix returned by \code{\link{gentourneybracket}} that ties each row in designmat to the corresponding choice set, including all potential tournament matchups
#' @param usepriorwinprobs logical, default TRUE - whether win probabilities should be calculated using prior parameter estimates or not. If FALSE, assumes each alternative has an equal chance of winning each round (may provide better protection against prior misspecification).
#'   Otherwise, uses prior estimates to calculate win probabilities for each choice set
#'
#' @return
#'   \item{setproplist}{Information matrix for each choice set. RENAME THIS OUTPUT IN THE FUTURE TO BETTER MATCH CONTENTS.}
#'   \item{winprobvect}{Vector of win probabilities for each profile in matchupframe}
#'   \item{occurprobvect}{A vector of the probability of each choice set occurring.}
#'   \item{info_mat}{Information matrix of the tournament design}
#' @export
#'
#' @description Calculates metrics for choice tournament designs based on the prior parameter estimates and the tournament bracket.
#'   The information matrix is the sum of the information matrix for each possible choice set in the tournament multiplied by its probability of occurence.
#'
#' @examples
getchoicemodinfoprobs <- function(modmat, paramestimates, matchupframe, usepriorwinprobs = TRUE){

  #Calculate the information matrix for each question pair and place them in a list
  choiceefflist <- lapply(unique(matchupframe$uniquesetid), function(x){

    selectvect <- matchupframe$matrixrowid[matchupframe$uniquesetid == x]
    return(d_effchoice(CurrentMatrix = modmat[selectvect,], altvect = rep(1, times = length(selectvect)), paramestimates = paramestimates, returninfomat = TRUE))}
  )

  #Extract the win probabilites from each matchup and calculate resulting occurrence probabilities for each matchup
  if(usepriorwinprobs){

    #Use prior parameter estimates to calculate win probabilities
    winprobs <- c(sapply(choiceefflist, "[[", "probvect"))

  }else{

    #Assume all win probabilities are equal for every choice set (may provide better protection against prior misspecification)
    winprobs <- c(sapply(choiceefflist, function(x) rep(1/length(x$probvect), length(x$probvect))))

  }

  #Calculate full occurrence probability vector:
  occurprobs <- matchupprobs(matchupframe = matchupframe, winprobs = winprobs)$occurprobs

  #Calculate weighted information matrix by multiplying info mat for each possible matchup with chance of that matchup happening
  info_mat <- Reduce("+", Map('*', lapply(choiceefflist, "[[", "info_mat"), occurprobs[duplicated(matchupframe$uniquesetid)]))

  output <- list("setproplist" = choiceefflist,
                 "winprobvect" = winprobs,
                 "occurprobvect" = occurprobs,
                 "info_mat" = info_mat)

  return(output)

}
