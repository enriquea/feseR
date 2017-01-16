
library(caret)


# # Generalized partial least squares
# plsProfile <- rfe(x = trainDescr, 
#                y = as.factor(trainClass), 
#                sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
#                rfeControl = rfeControl(functions = caretFuncs, method="cv", number=5),
#                method = 'gpls')
# 
# # Generalized linear model
# lmProfile <- rfe(x = trainDescr, 
#                 y = as.factor(trainClass), 
#                 sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
#                 rfeControl = rfeControl(functions = lmFuncs, method="cv", number=5))
# 
# # K-nearest neighbors analysis
# knnProfile <- rfe(x = trainDescr, 
#                  y = as.factor(trainClass), 
#                  sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
#                  rfeControl = rfeControl(functions = caretFuncs, method="cv", number=5),
#                  method = 'knn')
# 
# # Linear discriminant analysis
# ldaProfile <- rfe(x = trainDescr, 
#                   y = as.factor(trainClass), 
#                   sizes = c(1:5,10,20,30,50,100,200,300,400,500), 
#                   rfeControl = rfeControl(functions = ldaFuncs, method="cv", number=5))