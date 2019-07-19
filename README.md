
[![Travis build status](https://travis-ci.org/enriquea/feseR.svg?branch=master)](https://travis-ci.org/enriquea/feseR)


# feseR: Feature Selection in R


## Introduction

We provide here a R package which combine multiple Feature Selection (FS) methods in a workflow for analizing high-dimentional omics data. The different feature selection steps can be classificated in: i) Univariate (Correlation filter and Gain Information), ii) Multivariate (Principal Component Analysis and Matrix Correlation based) and iii) Recursive Feature Elimination (wrapped up with a Machine Learning algorithm, e.g. Random Forest). The goal is to essemble the different steps in an efficient workflow to perform feature selection in the contex of classification and regression problems.

## How to install

The first step is to install `devtools`:  

    install.packages("devtools")
    library(devtools)
   
Then, we can install the package using: 

    install_github("enriquea/feseR")
    library(feseR)


## How to use

We provide [here](https://github.com/enriquea/feseR/blob/master/vignettes/feser.pdf) some examples for illustrating how to use the package.

## This library has been used in:

Enrique Audain, Yassel Ramos, Henning Hermjakob, Darren R. Flower, Yasset Perez-Riverol; Accurate estimation of isoelectric point of protein and peptide based on amino acid sequences, Bioinformatics, Volume 32, Issue 6, 15 March 2016, Pages 821–827, [article](https://academic.oup.com/bioinformatics/article/32/6/821/1744386/Accurate-estimation-of-isoelectric-point-of)

## How to cite

If you find useful this tool in your work, you could want citing us:
Perez-Riverol Y, Kuhn M, Vizcaíno JA, Hitz M-P, Audain E (2017) Accurate and fast feature selection workflow for high-dimensional omics data. PLoS ONE 12(12): e0189875. https://doi.org/10.1371/journal.pone.0189875

## Mainteiner

Enrique Audain ([enriquea](https://github.com/enriquea))

Dmitry Rychkov ([drychkov](https://github.com/drychkov))

Yasset Perez-Riverol ([ypriverol](https://github.com/ypriverol))
