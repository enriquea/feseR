
#' plotPCVariances
#'
#' This function plots the Standard Desviation, Proportion of Variances and Accumulative Variance
#' for each Principal Component of the dataset, it returns a set of plots
#' 
#' @param dat dataframe: first col: PCs index
#' @param na.rm ignore the null values
#'
plotPCVariances <- function(dat, na.rm = TRUE, ...) {
  #retrive colum names
  nm  <-  names(dat)
 
   plots <- list()  # new empty list
 
    for (i in 1:length(nm)) {
    
        if(nm[1L] != nm[i]){
         
        plot <- ggplot()
         
          if(nm[i]=='Standard_desviation'){
            plot <- ggplot(dat) + 
            geom_point(aes_string(x=nm[1L], y = nm[i]), alpha = .5, stat = 'identity', colour = 'black') +
              xlab(nm[1L]) +
              ylab(nm[i]) +
              theme_bw()
          }
          if(nm[i]=='Proportion_of_variance'){
            plot <- ggplot(dat) + 
            geom_bar(aes_string(x=nm[1L], y = nm[i]), alpha = .5,  stat = 'identity', colour = 'blue') +
              xlab(nm[1L]) +
              ylab(nm[i]) +
              theme_bw()
          }
          if(nm[i]=='Accumulative_variance'){
            plot <- ggplot(dat) + 
            geom_line(aes_string(x=nm[1L], y = nm[i]), alpha = .5, stat = 'identity', colour = 'orange') +
              xlab(nm[1L]) +
              ylab(nm[i]) +
              theme_bw()
          }
          
          plot +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = "black"))
  
          plots[[i]] <- plot # add each plot into plot list
     
     }
  }
  return (plots)
}



#' multiplot
#'
#' This function allows to plot different ggplot charts in the same grids
#'
#' @param plotlists this param allow iterate thorugth a set of plots to
#' @param file if the user wants to plot to a file, it can use the file paramterer
#' @param cols Number of columns to be use.
#' @param layout the layout to be use

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#' plotMultipleSeries
#'
#' This fucntion plots multiple error bars
#'
#' @param dat the data.frame with the follow format (four column): series, expect var, predict var, error
#'
plotMultipleSeries <- function (dat){
  # rename data
  names(dat) <- c('series', 'x', 'y', 'error')
  pd <- position_dodge(.3) # Save the dodge spec because we use it repeatedly
  plot <- ggplot(dat, aes(x=x, y=y, colour=series, group=series)) +
    geom_errorbar(aes(ymin=y-error, ymax=y+error),
                  width=.2, size=0.25, colour="black", position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=2.5)
  return (plot)
}

