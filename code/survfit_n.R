
####################### EFFECTIVE N FOR SURVIVAL AT TIME #####################

survfit_n <- function( sf, 
                       cox = NULL,
                       round = TRUE, 
                       uncensored = TRUE, 
                       bounds = TRUE,
                       mod = TRUE
){
  # Check if, when the sf object is a survfit prediction of a Cox model, the cox
  # model is also provided to the function
  if( all("survfitcox" %in% class(sf), class(cox) != "coxph") ){
    error( "Providing the coxph object is required when looking at effective
           sample size for a Cox model prediction.")
  }
  
  # First, we unwrap the survfit object into a dataframe with all relevant data
  df <- unclass(sf)[c("time", "surv", "n.risk", "n.event", "n.censor", "std.err")] |>  
    data.frame()
  
  # If it is a cox fit, we need to calculate the stabilizer
  if( "survfitcox" %in% class(sf) ){
    # Calculate the weighted risk set at each time
    w.risk <- NULL
    for( i in 1:length(df$time) ){
      index <- which( unclass(cox$y)[, 1] >= df$time[i] ) 
      w.risk[i] <- sum( exp( cox$linear.predictors[index] ) )
    }
    stabiliser <- df$n.risk / w.risk
    # Multiply by individual linear predictor
    stabiliser <- stabiliser * as.numeric(exp(cox$coefficients*cox_fit$newdata))
  } else stabiliser = 1


  # Regular Kaplan-Meier without strata
  if( length(sf$n) == 1 ){
    ess <- calculate_ess( df, sf$n, 
                          round = round, 
                          uncensored = uncensored, 
                          bounds = bounds, 
                          mod = mod)
    ess[["n.eff"]] <- ess[["n.eff"]]*stabiliser
    ess[["n.eff.mod"]] <- ess[["n.eff.mod"]]*stabiliser
    sf <- append( sf, ess )
    class(sf) <- "survfit"
    return(sf)
  }
  
  
  else{
    # Find the cut points of the different strata
    cut <- c(0, cumsum(sf$strata) )
    # Run once outside loop for list initializaion
    ess <- calculate_ess( df[(cut[1]+1):cut[2], ], sf$n[1], 
                          round = round, 
                          uncensored = uncensored, 
                          bounds = bounds, 
                          mod = mod )
    for( i in 2:length(sf$strata) ){
      ess <- Map(c, ess, calculate_ess( df[(cut[i]+1):cut[i+1], ], sf$n[i], 
                                        round = round, 
                                        uncensored = uncensored, 
                                        bounds = bounds, 
                                        mod = mod)) 
    }
    sf <- append( sf, ess )
    class(sf) <- "survfit"
    return( sf )      
  }
  class(sf) <- "survfit"
  return( sf )
}
