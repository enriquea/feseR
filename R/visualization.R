
#' plot_pca
#'
#' This function plots the Standard Desviation, Proportion of Variances, Accumulative Variance
#' for each Principal Component, it returns a set of plots.
#'
#' @param features A numeric matrix as input.
#' @param class Response variable as numeric vector.
#' @param center Must the features be centered (default TRUE).
#' @param scale  Must the features be scaled (default TRUE).
#' @param list.plot If TRUE (default), return a list the plots summarizing the PCA process. If FALSE, return a simple plot (PC1 vs. PC2).
#'
#' @return A list of plots.
#'
#' @export

plot_pca <- function(features, class, center = TRUE, scale = TRUE, list.plot = TRUE) {

  # Matrix input validation
  valid.matrix(features)

  # compute principal componets
  features.pca <- prcomp(features, center = center, scale. = scale)

  # plot PC1 vs PC2
  p <- ggbiplot::ggbiplot(features.pca, obs.scale = 1, var.scale = 1,
                         groups = as.factor(class), ellipse = TRUE, circle = FALSE, var.axes = FALSE) +
                         scale_color_discrete(name = 'class') +
                         # theme(legend.direction = 'horizontal', legend.position = 'top') +
                         theme_bw()+
                         theme(axis.text = element_text(size=12),
                               axis.title = element_text(size=14),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))

  # if list.plot FALSE, return a simple plot
  if(!isTRUE(list.plot)){
    return(p)
  }

  # parsing prcomp object
  summ_stat <- summary(features.pca)

  # getting dataframe with PCs variance information
  dat <- as.data.frame(t(summ_stat$importance))

  #add Principal Component index
  dat <- cbind(PCs=c(1:nrow(dat)), dat)
  colnames(dat) <- c('PCs','Standard_Desviation','Proportion_of_Variance','Accumulative_Variance')

  p1 <- ggplot2::ggplot(dat) +
    geom_point(aes_string(x='PCs', y = 'Standard_Desviation'), alpha = .5, stat = 'identity', colour = 'black') +
    xlab('Principal Component') +
    ylab('Standard Deviation') +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  p2 <- ggplot2::ggplot(dat) +
    geom_bar(aes_string(x='PCs', y = 'Proportion_of_Variance'), alpha = .5,  stat = 'identity', colour = 'blue') +
    xlab('Principal Component') +
    ylab('Proportion of Variance') +
    theme_bw() +
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  p3 <- ggplot2::ggplot(dat) +
    geom_line(aes_string(x='PCs', y = 'Accumulative_Variance'), alpha = .5, stat = 'identity', colour = 'orange') +
    xlab('Principal Component') +
    ylab('Accumulative Variance') +
    theme_bw()+
    theme(axis.text = element_text(size=12),
          axis.title = element_text(size=14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))

  # add each plot into plot list
  plots <- list(pc.plot = p, st.dev.plot = p1, prop.var = p2, acc.var = p3)

  return (plots)
}

#' plot_corr
#'
#' Simple function to plot a correlation matrix between input features.
#'
#' @param features A numeric matrix as input.
#' @param corr.method Specifies the type of correlation to be computed (pearson, spearman or kendall)
#'
#' @return A correlation plot.
#'
#' @export

plot_corr <- function(features, corr.method = 'pearson') {

  # Matrix input validation
  valid.matrix(features)

  # compute correlation
  featureCorr <- cor(features, method = corr.method)

  # plotting corr matrix
  plot <- corrplot::corrplot(featureCorr,
                             order = "hclust",
                             tl.pos = 'n',
                             tl.cex = 0.8,
                             cl.cex = 1.0) # Make plot
  return(invisible(plot))

}


#' multiplot
#'
#' This function allows to plot different ggplot charts in the same grid.
#'
#' @param plotlist list of ggplot objects
#' @param file if the user wants to plot to a file, it can use the file paramterer
#' @param cols Number of columns to be use.
#' @param layout the layout to be use
#' 
multiplot <- function(plotlist = NULL, file, cols = 1, layout = NULL) {

  # plots <- c(list(...), plotlist)
  plots <- plotlist
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

