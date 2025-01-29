
####################### PLOT EFFECTIVE SAMPLE SIZE OF KM #######################
#Version based on GGSURVPLOT

plot_km_eff2 <- function( 
    sf, 
    title = NA, 
    legend = c(0.15, 0.4),
    bounds = TRUE,
    ylim = NA,
    mod = FALSE,
    both = FALSE,
    xlab = "Time",
    ylab = "Survival Probability",
    mark = F,
    col = c("#2E7691", "#C2666B", "#c6aa2c", "#37293F", "darkslategray3", "darkgoldenrod1", "chartreuse", "darkslateblue", "darkred", "bisque2"),
    risk = F,
    ess = F,
    conf.int = F,
    legend.title = "Treatment",
    grid_size = 5
){
  
  # We need the number of strata
  if( "strata"%in% names(sf) ){
    # Remove the variable name from strata
    names(sf$strata) <- gsub("^[^=]+=", "", names(sf$strata))
    
    # Select correct number of colors
    color_values <- setNames(col[1:length(sf$strata)], names(sf$strata))
    
    # Space
    if( both ){ 
      space = 5 + length(sf$strata)*2 + 2
    } else space = 5 + length(sf$strata)
    
    # Cut points
    cut <- c(0, cumsum(sf$strata) )
  } else{
    if( both ){ 
      space = 9
    } else space = 6
    cut <- c(0, sf$n )
  }
  
  # The GGSURVFIT object cannot take any value in a sf object, but only one of
  # a selection, so we need to overwrite the upper and lower conf bounds to get
  # the effective N to show up.
  sf$upper <- sf$n.eff
  sf$lower <- sf$n.eff.mod
  
  # Now we will plot the basic Kaplan-Meiers
  if( "strata"%in% names(sf) ){
    plot <- sf |>
      ggsurvfit( linewidth = 1.3 ) +
      scale_ggsurvfit() +
      theme_minimal() +
      scale_x_continuous(breaks = pretty(sf$time, n = grid_size)) +
      scale_color_manual(values = color_values) +
      labs(x = xlab, y = ylab) +
      theme(legend.position.inside = legend, panel.grid.major.x = element_blank()) +
      theme(text = element_text(size = 12), plot.title = element_text( size = 18, hjust = 0.3)) 
  }
  else{
    plot <- sf |>
      ggsurvfit( linewidth = 1.3, color = "#2E7691" ) +
      scale_ggsurvfit() +
      theme_minimal() +
      scale_x_continuous(breaks = pretty(sf$time, n = grid_size)) +
      labs(x = xlab, y = ylab) +
      theme(legend.position.inside = legend, panel.grid.major.x = element_blank()) +
      theme(text = element_text(size = 12), plot.title = element_text( size = 18, hjust = 0.3)) 
  }
    
  if( !is.na(title) ){
    plot <- plot +
      labs( title = title )
  }
  
  # Add censoring marks if required
  if( mark ) 
    plot <- plot + add_censor_mark()
  
  if( both ){
    if( mod ){
      plot <- plot +
        add_risktable(
          risktable_stats = "{conf.low} ({n.risk})",
          stats_label = c("Mod. Effective Sample Size (At Risk)")
        )
    } else {
      plot <- plot +
        add_risktable(
          risktable_stats = "{conf.high} ({n.risk})",
          stats_label = c("Effective Sample Size (At Risk)")
        )
    }
  }
  return(plot)
}
