
##### TURN A SURVFIT OBJECT INTO A DATAFRAME WITH A CERTAIN TIME SCALE #####

sf_to_df <- function( sf, time ){
  # Turn list in to df
  df <- unclass(sf)[c("time", "surv", "n.risk", "n.uncensor", "n.eff", "n.eff.mod", "n.lower", "n.upper")] |>  
    data.frame()
  
  # merge with time grid
  df <- merge( df, data.frame(time = time ), by = "time", all = TRUE, sort = TRUE )
  
  # Fill the dataframe downwards
  df <- fill(df, 2:8, .direction = "downup")
  
  return(df[ df$time %in% time, ])
  
}