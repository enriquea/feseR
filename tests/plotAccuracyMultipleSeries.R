library(plyr)
library(ggplot2)
library(pIR)

metrics.list <- list(RFE_CV10, PCA_RFE_C10, X2_RFE_C10, X2_PCA_RFE_C10)
metrics.merged <- do.call(rbind, metrics.list)
metrics.merged <- subset(metrics.merged, Variables <= 100)

write.table (metrics.merged, file = 'analysis/benchmark_accuracy_random_forest_classifiers.txt', row.names = FALSE, sep = '\t')

dat <- metrics.merged

# rename data
# names(dat) <- c('experimental', 'predicted', 'error', 'method')

pd <- position_dodge(.3) # Save the dodge spec because we use it repeatedly
plot <- ggplot(dat, aes(x=Variables, y=Accuracy, colour=method, group=method)) +
  #geom_errorbar(aes(ymin=Accuracy-AccuracySD, ymax=Accuracy-AccuracySD), width=.2, size=0.25, colour="black", position=pd) +
  stat_smooth(size=0.5, span = 0.5, alpha = 0.2) +
  #geom_point(size=0.5) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"))

png(filename = 'figures/benchmark_accuracy_genearray_dataset.png', width = 800, height = 800)
plot
dev.off()


