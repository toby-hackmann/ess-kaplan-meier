
####################### PLOT EFFECTIVE SAMPLE SIZE OF KM #######################

plot_km_eff <- function( 
    sf, 
    title = "Comparing effective sample size to other data measures", 
    legend.pos = c(0.15, 0.4),
    bounds = TRUE,
    ylim = NA,
    mod = FALSE,
    both = FALSE,
    xlab = "Time",
    ylab = "Survival",
    mark = F,
    col = c("#2E7691", "#C2666B", "#c6aa2c", "#37293F", "darkslategray3", "darkgoldenrod1", "chartreuse", "darkslateblue", "darkred", "bisque2"),
    risk = F,
    ess = F,
    conf.int = F,
    legend.title = "Treatment",
    min.grid = 5
){
  
  # We need the number of strata
  if( "strata"%in% names(sf) ){
    # Remove the variable name from strata
    names(sf$strata) <- gsub("^[^=]+=", "", names(sf$strata))
    
    # Select correct number of colors
    col <- col[1:length(sf$strata)]
    
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
 
  # Starting point
  line = 5
  
  # Now we will plot the basic Kaplan-Meiers
  par(mfrow=c(1,1),mar = c(space,6,2,2)+0.1, bty = 'n') # defining plot parameters
  plot(sf, 
       xlab=xlab, 
       ylab=ylab, 
       lwd=2, 
       mark.time=mark, 
       col=col, 
       main = title, 
       conf.int = conf.int)
  # Legend
  if( "strata"%in% names(sf) )
    legend( x = legend.pos[1], 
          y = legend.pos[2], 
          title = legend.title, 
          legend = names(sf$strata), 
          col = col, 
          lwd = 2, 
          box.lty = 0, 
          bty = "n", 
          cex = 0.7 )
  
  ### FUNCTION THAT TAKES LOOKS UP THE LAST KNOWN VALUE BEFORE THE GRID TICK
  .SetCount <- function(timeindex, time, number) {
    show <- number[1]
    for (t in timeindex[-1]){
      show <- c(show, number[ max(which(time <= t))])
    }
    return(show)
  }
  
  # Number ticks on the x-axis
  grid <- axTicks(1)
  # Distance between ticks
  tick <- grid[2]-grid[1]
  
  # Now we choose how to present
  if ( both ){ # If both, we add two tables
    mtext("Effective Sample Size (At Risk):", 1, line=4, at= -tick/2, cex = 0.9, adj = 0)
    for( i in 1:length(sf$strata) ){
      mtext( paste( names(sf$strata)[i], ":"), side = 1, line = line, at = -tick/3, cex = 0.9, col = col[i], adj = 1)
      
      r = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.risk[(cut[i]+1):cut[i+1]]) 
      n = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.eff[(cut[i]+1):cut[i+1]])
      
      mtext( paste(n, "(", r,")") , side = 1, line = line, at = grid-tick/5, cex = 0.8, col = col[i], adj = 0)
      line = line + 1
    }
    line = 10
    mtext("Modified Effective Sample Size (At Risk):", 1, line=9, at= -tick/2, cex = 0.9, adj = 0)
    for( i in 1:length(sf$strata) ){
      mtext( paste( names(sf$strata)[i], ":"), side = 1, line = line, at = -tick/3, cex = 0.9, col = col[i], adj = 1)
      
      r = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.risk[(cut[i]+1):cut[i+1]]) 
      n = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.eff.mod[(cut[i]+1):cut[i+1]])
      
      mtext( paste(n, "(", r,")") , side = 1, line = line, at = grid-tick/5, cex = 0.8, col = col[i], adj = 0)
      line = line + 1
    }
  }
  
  else if( risk & ess ){ # If we want both the risk set and effective N
    if( mod ) mtext("Modified Effective Sample Size (At Risk):", 1, line=4, at= -tick/2, cex = 0.9, adj = 0)
    else mtext("Effective Sample Size (At Risk):", 1, line=4, at= -tick, cex = 0.9, adj = 0)
    for( i in 1:length(sf$strata) ){
      mtext( paste( names(sf$strata)[i], ":"), side = 1, line = line, at = -tick/3, cex = 0.9, col = col[i], adj = 1)
      
      r = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.risk[(cut[i]+1):cut[i+1]]) 
      if( mod ) n = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.eff.mod[(cut[i]+1):cut[i+1]])
      else n = .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.eff[(cut[i]+1):cut[i+1]])
        
      mtext( paste(n, "(", r,")") , side = 1, line = line, at = grid-tick/5, cex = 0.8, col = col[i], adj = 0)
      line = line + 1
    }
    
  }
  
  else if( risk ){ # If we want only risk set
    mtext("At Risk:", 1, line=4, at=-tick/2, cex = 0.9, adj = 0)
    for( i in 1:length(sf$strata) ){
      mtext( paste( names(sf$strata)[i], ":"), side = 1, line = line, at = -tick/5, cex = 0.9, col = col[i], adj = 1)
      mtext( .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.risk[(cut[i]+1):cut[i+1]]) , side = 1, line = line, at = grid, cex = 0.8, col = col[i], adj = 0)
      line = line + 1
    }
  }
  
  else if ( ess ) { # If we want only effective sample size
    mtext("Effective Sample Size", 1, line=4, at=-tick/2, cex = 0.9, adj = 0)
    for( i in 1:length(sf$strata) ){
      mtext( paste( names(sf$strata)[i], ":"), side = 1, line = line, at = -tick/5, cex = 0.9, col = col[i], adj = 1)
      mtext( .SetCount(grid, sf$time[(cut[i]+1):cut[i+1]], sf$n.eff[(cut[i]+1):cut[i+1]]) , side = 1, line = line, at = grid, adj = 0, cex = 0.8, col = col[i])
      line = line + 1
      
    }
  }
}
