
##### TRIAL REPORT FUNCTION THAT SIMS REPORTING EFFECTIVE N #####

trial_report <- function( 
    
  dist = "weibull", # options are weibull, exponential or gompertz
  gammas = 1, # as in simsurv
  lambdas = 1, # as in simsurv
  groups = 1, # number of groups
  x = NULL, # dataframe containing covariate values for each individual
  n_patient = 100, # not required if x is a thing, otherwise this will be #ids
  betas = NULL, # coefficients of x, if provided
  t_enroll = NULL, # when enrollment ends
  zero_enroll = 0, # number of patients enrolled at t=0
  t_outcome = NULL, # time of primary endpoint
  iter = 10, # number of iterations to average over
  points = 20, # number of timepoints to evaluate
  title = "Modified effective sample size as a function of reporting time",
  xlab = "Reporting time after start enrollment",
  ylab = "Percentage of total sample size",
  mod = TRUE, # toggles modified and regular effective sample size
  conf = FALSE,
  conf.int = 95,
  proportion = TRUE
  
){
  
  
  ###################### WARNINGS ######################
  if( is.null(t_enroll) ) error( "The latest enrollment time is required." )
  if( is.null(t_outcome) ) error( "A primary endpoint is required." )
  if( (length(lambdas) == length(gammas)) & (length(lambdas) != groups) ) groups = length(lambdas)
  if( length(lambdas) != groups | length(gammas) != groups ) error("Number of elements in gammas and lambdas needs to equal number of groups.")
  if( "treat" %in% names(x) & nlevels(x$treat) > 2 ) error("We can only do two levels at the moment.")
  if( groups > length(dist) & length(dist) == 1 ) dist <- rep(dist, groups) 
  
  
  ############################ INITIALIZATIONS ###########################
  # If you only provide the number of patients instead of df x, create df
  if( is.null(x) ) x = data.frame(id = 1:n_patient) 
  
  # Create a timesequence of length points where it will be analysed
  time = seq(t_outcome, t_outcome+t_enroll, length.out = points+1 )
  
  # Create the matrix to save effective N in, one row for each group
  n <- array(0, dim = c( groups, length(time), iter ) ) 
  if( "treat" %in% names(x) ) n <- array(0, dim = c(2, length(time), iter) ) 
  
  
  
  ############################ LOOP ###########################
  # Loop over the different groups
  for ( g in 1:groups ){
    # Iterations to average over
    for ( i in 1:iter ){
      
      # Generate the survival time for different groups
      survival <- simsurv( dist=dist[g], 
                           gammas=gammas[g], 
                           lambdas=lambdas[g], 
                           x=x,
                           maxt=t_outcome,
                           betas = betas)
      
      # Generate enrollment time
      # First we get the number of patients that we want zero-enrolled out
      if( zero_enroll > 0 ) survival$tstart[1:zero_enroll] = rep(0, zero_enroll)
      # Then add the rest by drawing from uniform distr
      survival$tstart[(zero_enroll + 1):nrow(survival)] <- runif( nrow(survival)-zero_enroll, 0,  t_enroll)
      
      # Survival in 'real' time since start study
      survival$tstop <- survival$tstart + survival$eventtime
      
      # Add treatment if necessary and then do a different method
      if( "treat" %in% names(x) ){ 
        survival$treat = x$treat
        
        for( t in time ){
          df <- data.frame(
            time = ifelse( survival$tstop < t, survival$eventtime, survival$eventtime - (survival$tstop - t)),
            status = ifelse( survival$tstop < t, survival$status, 0),
            treat = as.numeric(as.factor(survival$treat))
          )
          obj1 <- survfit( Surv(time, status) ~ 1, data = df[ df$treat == 1, ]) |> 
            survfit_n()
          obj2 <- survfit( Surv(time, status) ~ 1, data = df[ df$treat == 2, ]) |> 
            survfit_n()
          
          # Extract effective sample size 
          if( mod ){
            n[1, which(time == t), i] <- obj1$n.eff.mod[length(obj1$n.eff.mod)]
            n[2, which(time == t), i] <- obj2$n.eff.mod[length(obj2$n.eff.mod)]
          }
          else{
            n[1, which(time == t), i] <- obj1$n.eff[length(obj1$n.eff)]
            n[2, which(time == t), i] <- obj2$n.eff[length(obj2$n.eff)]
          }
        }
      }
      
      # Create the correct survival dataframe for each analysis time based
      else{
        for( t in time ){
          df <- data.frame(
            time = ifelse( survival$tstop < t, survival$eventtime, survival$eventtime - (survival$tstop - t)),
            status = ifelse( survival$tstop < t, survival$status, 0) 
          )
          obj <- survfit( Surv(time, status) ~ 1, data = df) |> 
            survfit_n()
          
          # Extract effective sample size and add each with weight 1/iter
          if( mod ) n[g, which(time == t), i] <- obj$n.eff.mod[length(obj$n.eff.mod)]
          else n[g, which(time == t), i] <- obj$n.eff[length(obj$n.eff)]
        }
      }
    }
  }
  
  quantile <- (1-conf.int/100)/2
  
  # Initialize an empty data frame
  df_long <- data.frame()
  
  # Loop through each group dynamically and construct the long format data
  if("treat" %in% names(x)) {
    for (i in 1:nrow(n)) {
      group_label <- paste("Treat", i)
      if( proportion ){
      df_long <- rbind(df_long, data.frame(
        time = time,
        number = (rowMeans(n[i, , ])* (100/nrow(df[df$treat == i, ]))),
        group = group_label, 
        lower = (rowQuantiles(n[i, , ], probs = quantile)* (100/nrow(df[df$treat == i, ]))),
        upper = (rowQuantiles(n[i, , ], probs = 1-quantile)* (100/nrow(df[df$treat == i, ])))
      ))
      } else {
        df_long <- rbind(df_long, data.frame(
          time = time,
          number = (rowMeans(n[i, , ])),
          group = group_label,
          lower = (rowQuantiles(n[i, , ], probs = quantile) ),
          upper = (rowQuantiles(n[i, , ], probs = 1-quantile) )
        ))
      }
    }
  } else{
    for (i in 1:nrow(n)) {
      group_label <- paste("Group", i)
      df_long <- rbind(df_long, data.frame(
        time = time,
        number = (rowMeans(n[i, , ])),
        group = group_label,
        lower = (rowQuantiles(n[i, , ], probs = quantile) ),
        upper = (rowQuantiles(n[i, , ], probs = 1-quantile) )
      ))
    }
    if( proportion ){
      df_long$number <- df_long$number * (100/n_patient)
      df_long$number <- df_long$number * (100/n_patient)
      df_long$number <- df_long$number * (100/n_patient)
    }
  }

  
  
  
  ########################### PLOT ###########################
  
  # Define colors
  colors <- c("#2E7691", "#C2666B", "#c6aa2c", "#9f84af")
  linetypes <- c("Mean" = "solid",  "Confidence" = "dashed")
  
  # Initialize the ggplot object
  if( proportion ) plot <- ggplot(df_long, aes(x = time, y = percentage, color = group, linetype = "Mean"))
  else plot <- ggplot(df_long, aes(x = time, y = number, color = group, linetype = "Mean"))
  
  # Additions
  plot <- plot + 
    labs(x = xlab, 
         y = ylab,
         title = title) +
    geom_line(linewidth = 1.3) #+
    #geom_hline(yintercept = max(n), linetype = "dashed", color = "black", linewidth = 1.3)
    
  if( conf ){
    plot <- plot +
      geom_line( aes( x = time, y = upper, color = group, linetype = "Confidence")) +
      geom_line( aes( x = time, y = lower, color = group, linetype = "Confidence"))
  }
  
  # Customize the theme and scales
  plot <- plot + 
    theme_minimal() +
    scale_color_manual(values = setNames(colors[1:length(unique(df_long$group))], unique(df_long$group)), 
                       name = "Legend") +
    scale_linetype_manual(values = linetypes, 
                          name = "Legend") +
    ylim(c(0, max(n)))
  
  return( plot )
}