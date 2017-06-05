# Feature Selection in R (feseR)

## Introduction

We provide here some Feature Selection (FS) workflows to analyze high-dimentional omics data in R environment. The different feature selection steps can be classificated in: Univariate (Correlation filter and Gain Information), Multivariate (Principal Component Analysis and Matrix Correlation based) and Recursive Feature Elimination (wrapped up with a Machine Learning algorithm). The goal is to essemble the different steps in a efficient workflow to perform feature selection task in the contex of classification and regression problems.

## How to use

The users could use the tool as simple R script at the moment. Thus, for example, for combining Univariate Correlation (X2) -> Multivariate correlation (MC) -> Recursive Feature Selection wrapped with Random Forest (RFE-RF), see implementation [here](https://github.com/enriquea/feseR/blob/master/workflows/X2_CM_RFE_RF.R). Optionally, other combination to test could be Univariate Correlation (X2) -> Principal Component Analysis (PCA) -> Recursive Feature Selection wrapped with Random Forest (RFE-RF) [here](https://github.com/enriquea/feseR/blob/master/workflows/X2_PCA_RFE_RF.R). Each workflow generates a report summarizing stats from training and testing phases.

## Dataset

We provide some example datasets (Transcriptomics and Proteomics, [here](https://github.com/enriquea/feseR/tree/master/data)). Some general description are listed bellow:

#### Classification

- GSE5325 (Analysis of breast cancer tumor samples using 2-color cDNA microarrays)
- TNBC (Label-free deep proteome analysis of 44 (samples and technical replicates) human breast specimens)
- GSE48760 (Transcriptomics analysis of left ventricles of mouse subjected to an isoproterenol challenge)
- GSE6919/GPL8300 (Expression data from normal and prostate tumor tissues)
- GSE6919/GPL92 (Expression data from normal and prostate tumor tissues)
- GSE6919/GPL93 (Expression data from normal and prostate tumor tissues)

#### Regression
- Isopoint (Peptides fractionated by IPG-HPLC and analyzed by Mass spectrometry)

## How to cite
If you find useful this tool in your work, you could want citing us:
Yasset Perez-Riverol, Max Kuhn, Juan Antonio Vizcaíno, Marc-Phillip Hitz, Enrique Audain. 2017. *Accurate and Fast feature selection workflow for high-dimensional omics data*. doi: https://doi.org/10.1101/144162.