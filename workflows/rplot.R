library(dplyr)
library(plyr)
library(reshape2)
library(pIR)

#' plotMultipleSeries
#'
#' This fucntion plots multiple error bars
#'
#' @param dat the data.frame with the follow format (four column): expect var, predict var, error, series
#'
plotMultipleSeries <- function (dat){
    # rename data
    names(dat) <- c('experimental', 'predicted', 'error', 'method')
    pd <- position_dodge(.3) # Save the dodge spec because we use it repeatedly
    plot <- ggplot(dat, aes(x=experimental, y=predicted, colour=method, group=method)) +
        geom_errorbar(aes(ymin=predicted-error, ymax=predicted+error),
                      width=.8, size=0.25, colour="black", position=pd) +
        geom_line(position=pd, size=0.2) +
        geom_point(position=pd, size=1.5) +
        theme_bw() +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour="black"))
    return (plot)
}

df <- data.frame(peptideClass, predSVMFS3, predSVMFS7, predSVMFS10, prednoCSVM3, prednoCSVM7, prednoCSVM10, predSVMnoFS, predSVMnoFSC)
colnames(df) <- c("experimental", "X2-CM-RBE-SVM-CV3", "X2-CM-RBE-SVM-CV7", "X2-CM-RBE-SVM-CV10", "RBE-SVM-CV3", "RBE-SVM-CV7", "RBE-SVM-CV10", "SVM", "X2-CM-SVM")



dat <- df
dat$`X2-CM-RBE-SVM-CV7` <- NULL
dat$`X2-CM-RBE-SVM-CV10`<- NULL
dat$`RBE-SVM-CV7` <- NULL
dat$`RBE-SVM-CV10` <- NULL

mergedata <- data.frame()
nm <- colnames(dat)
for (i in 1:length(nm)) {
    if(nm[1] != nm[i]){
        method <- nm[i]
        newData <- data.frame(x = dat[1L], y = dat[i])
        colnames(newData) <- c('experimental', 'predicted')
        fractions_stats <- ddply(newData, ~experimental, summarise,
                                 predicted = mean(predicted),
                                 error = rmse(experimental, predicted),
                                 method = method)
        mergedata <- rbind(mergedata, fractions_stats)
    }
}

png(filename = 'svm_methods_plot.png', width = 800, height = 800)
plotMultipleSeries(mergedata)
dev.off()

df <- data.frame(peptideClass, predSVMFS3, predSVMFS7, predSVMFS10, prednoCSVM3, prednoCSVM7, prednoCSVM10, predSVMnoFS, predSVMnoFSC)

colnames(df) <- c("experimental", "SVM-FS-CV3", "SVM-FS-CV7", "SVM-FS-CV10", "SVM-NONCF-CV3", "SVM-NONCF-CV7", "SVM-NONCF-CV10", "SVM-NONFS-NOCV", "SVM-NONBFS-NOCV")

save(df,file="SVMPredcitions.Rda")


