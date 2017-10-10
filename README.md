[![Travis build status](https://travis-ci.org/enriquea/feseR.svg?branch=master)](https://travis-ci.org/enriquea/feseR)

# feseR: Feature Selection in R

## Introduction

We provide here a R package which combine multiple Feature Selection (FS) methods in a workflow for analizing high-dimentional omics. The different feature selection steps can be classificated in: Univariate (Correlation filter and Gain Information), Multivariate (Principal Component Analysis and Matrix Correlation based) and Recursive Feature Elimination (wrapped up with a Machine Learning algorithm, i.e. Random Forest). The goal is to essemble the different steps in a efficient workflow to perform feature selection task in the contex of classification and regression problems.

## How to install

The first step is to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then, we can install the package using: 

    install_github("enriquea/feseR")
    library(feseR)


## How to use

We provide [here](https://github.com/enriquea/feseR/blob/master/vignettes/feseR.html) some examples for illustrating how to use the package.

## This library has been used in:

Enrique Audain, Yassel Ramos, Henning Hermjakob, Darren R. Flower, Yasset Perez-Riverol; Accurate estimation of isoelectric point of protein and peptide based on amino acid sequences, Bioinformatics, Volume 32, Issue 6, 15 March 2016, Pages 821–827, [article](https://academic.oup.com/bioinformatics/article/32/6/821/1744386/Accurate-estimation-of-isoelectric-point-of)

## How to cite

If you find useful this tool in your work, you could want citing us:
Yasset Perez-Riverol, Max Kuhn, Juan Antonio Vizcaíno, Marc-Phillip Hitz, Enrique Audain. 2017. *Accurate and Fast feature selection workflow for high-dimensional omics data*. doi: https://doi.org/10.1101/144162.