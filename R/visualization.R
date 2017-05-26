library(ggbiplot)
library(ggplot2)

#' plotPCVariances
#'
#' This function plots the Standard Desviation, Proportion of Variances, Accumulative Variance
#' for each Principal Component of the dataset, it returns a set of plots
#' 
#' @param prcomp prcomp: Principal component object from prcomp R function.
#' @param groups groups: Vector with the 'class variables' as the same size of the original data.
#'

plotPCVariances <- function(prcomp = NULL, groups = NULL) {
  
  # parsing prcomp object
  summ_stat <- summary(prcomp)
  
  # getting dataframe with PCs variance information
  dat <- as.data.frame(t(summ_stat$importance))
  
  #add Principal Component index
  dat <- cbind(PCs=c(1:nrow(dat)), dat) 
  colnames(dat) <- c('PCs','Standard_Desviation','Proportion_of_Variance','Accumulative_Variance')
  
  plots <- list()  # new empty list
  
  p1 <- ggplot(dat) + 
    geom_point(aes_string(x='PCs', y = 'Standard_Desviation'), alpha = .5, stat = 'identity', colour = 'black') +
    xlab('Principal Component') +
    ylab('Standard Deviation') +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  p2 <- ggplot(dat) + 
    geom_bar(aes_string(x='PCs', y = 'Proportion_of_Variance'), alpha = .5,  stat = 'identity', colour = 'blue') +
    xlab('Principal Component') +
    ylab('Proportion of Variance') +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  p3 <- ggplot(dat) + 
    geom_line(aes_string(x='PCs', y = 'Accumulative_Variance'), alpha = .5, stat = 'identity', colour = 'orange') +
    xlab('Principal Component') +
    ylab('Accumulative Variance') +
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  p <- ggbiplot(prcomp, obs.scale = 1, var.scale = 1,
                groups = as.factor(groups), ellipse = TRUE, circle = FALSE, var.axes = FALSE) +
    scale_color_discrete(name = 'class') +
    # theme(legend.direction = 'horizontal', legend.position = 'top') +
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"))
  # add each plot into plot list
  plots[[1]] <- p
  plots[[2]] <- p1 
  plots[[3]] <- p2 
  plots[[4]] <- p3 
  
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

