This repository contains the R scripts used in Dapas et al. (2019). Phenotypic clustering reveals distinct subtypes of polycystic ovary syndrome with novel genetic associations. (bioRxiv 814210; doi: https://doi.org/10.1101/814210)

The aforementioend manuscript detailing these methods and applications is currently under peer-review.


# PCOS_phenotype_clustering.R
This script is designed to read in a tab-delimited text file with trait values and assay method codes for the following traits by default: Age, BMI, Testosterone, SHBG, Insulin, Glucose, DHEAS, LH, and FSH.

The default traits, column names, and other variables are defined in the "Variables" section of this script. These variable can be modified as desired. 

This clustering method log-normalizes trait distributions, then adjusts by age and assay, then applies an inverse normal transformation (INT), in which the trait distributions are fit onto a normal distribution. Clusters are therefore primarily driven by multicollinearity of input variables.


# PCOS_subtype_classifer.R
This script builds and compares subtype classifiers modeled from adjusted quantitative trait data and clusters produced by PCOS_phenotype_clustering.R. The script compares the following methods: support vector machine, random forest, Ggussian mixed-model, and quadratic discriminant analysis. Each classifier is compared using 10-fold cross-validation. The classifier with the lowest error rate was then applied to a family-based PCOS cohort to identify and compare subtype membership.
