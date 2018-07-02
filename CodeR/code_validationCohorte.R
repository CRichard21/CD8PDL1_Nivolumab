######################################
## -------------------------------- ##
## ------ Analyse des RNASeq ------ ##
## -------------------------------- ##
######################################


## Set path to working directory

path.wd <- "my.working.directory"
setwd(path.wd)

## Load required packages

library(xlsx)
library(survminer)
library(survivalROC)
source("CodeR/code_cutoff.R")


## Import data

myData <- read.xlsx(file = "dataValidation.xlsx",
                    sheetIndex = 1, stringsAsFactors = F)

#myData <- myData[!(myData[,1]%in%c(1,4,16,15,11,56,35,34,23)),]
# myData <- myData[!(myData[,1]%in%c(1,4,16,15,11,23)),]

attach(myData)

## ------------------------------------------
## Analysis of  RNA data --------------------
## ------------------------------------------

## Association with PFS ----------------------

## CD8A RNA

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("./cutoffFinder/CD8A_valid"),
                         biomarker=CD8A_RNA,
                         outcome=NULL,
                         time=PFS,
                         event=evtPFS,
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier", "HR"),
                         nmin=10)

## CD274 RNA

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("./cutoffFinder/CD274_valid"),
                         biomarker=CD274_RNA,
                         outcome=NULL,
                         time=PFS,
                         event=evtPFS,
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier", "HR"),
                         nmin=10)


## Multivariate model  ---------------------
## -----------------------------------------

detach(myData)
myData <- subset(myData, Histology!="other")
attach(myData)


## Clinical model --------------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology)
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC


## CD8 PDL1 model --------------------------

myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeDGM+Histology+as.factor(CD8A_RNA>2.189)+as.factor(CD274_RNA>4.799))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS, evtPFS, marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

Roc <- val_pred

## IFNG model --------------------------

myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeDGM+Histology+as.factor(IFNG_RNA>8.985))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS, evtPFS, marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

Roc <- cbind(Roc, IFNG = val_pred)

## EIGS model --------------------------

myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeDGM+Histology+as.factor(EIGS_RNA>10.37))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS, evtPFS, marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

Roc <- cbind(Roc, EIGS = val_pred)

write.xlsx(cbind(PFS, Roc), file = "Roccurve.xlsx")
