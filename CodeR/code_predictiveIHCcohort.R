#############################################
## --------------------------------------- ##
## -- Analysis of IHC predictive cohort -- ##
## --------------------------------------- ##
#############################################


## Set path to working directory

path.wd <- "my.working.directory"
setwd(path.wd)

## Load required packages

library(xlsx)
library(survminer)
library(survivalROC)
source("CodeR/code_cutoff.R")

## Import data

myData <- read.xlsx(file = "dataIHCRNApredictive.xlsx",
                    sheetIndex = 1, stringsAsFactors = F)
attach(myData)
bestRECIST[bestRECIST=="SD"] <- "PR"


## ------------------------------------------
## Analysis of clinical data ----------------
## ------------------------------------------


## Sex

table(Sex)
prop.table(table(Sex))

table(Sex, bestRECIST)
prop.table(table(Sex, bestRECIST),2)
chisq.test(table(Sex, bestRECIST))

## Age

summary(AgeNivo)
sd(AgeNivo)

summary(AgeNivo[bestRECIST=="PD"])
sd(AgeNivo[bestRECIST=="PD"])

summary(AgeNivo[bestRECIST!="PD"])
sd(AgeNivo[bestRECIST!="PD"])

wilcox.test(AgeNivo~bestRECIST)

## WHO performance status

table(WHO.performance.status)
prop.table(table(WHO.performance.status))

table(WHO.performance.status, bestRECIST)
prop.table(table(WHO.performance.status, bestRECIST),2)
fisher.test(table(WHO.performance.status, bestRECIST))

## Histology

table(Histology)
prop.table(table(Histology))

table(Histology, bestRECIST)
prop.table(table(Histology, bestRECIST),2)
chisq.test(table(Histology[Histology!="other"], bestRECIST[Histology!="other"]))

## Tumor stage

table(Tumor.stage)
prop.table(table(Tumor.stage))

table(Tumor.stage, bestRECIST)
prop.table(table(Tumor.stage, bestRECIST),2)
chisq.test(table(Tumor.stage, bestRECIST))

## Smoker

table(Smoker)
prop.table(table(Smoker))

table(Smoker, bestRECIST)
prop.table(table(Smoker, bestRECIST),2)
fisher.test(table(Smoker, bestRECIST))


## ------------------------------------------
## Analysis of IHC and RNA data -------------
## ------------------------------------------


PDL1_Sp142 <- factor(ifelse(PDL1_IC_Sp142<2&PDL1_TC_Sp142<2,"low","high"), levels = c("low", "high"))

## Association with PFS ----------------------

## CD8 IHC

fit <- survfit(Surv(PFS, evtPFS)~as.factor(CD8_IHC>7.676), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

chisq.test(table(CD8_IHC>7.676, bestRECIST))

## CD8A RNA

fit <- survfit(Surv(PFS, evtPFS)~as.factor(CD8A_RNA>5.0), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fisher.test(table(CD8A_RNA>5.0, bestRECIST))

## PDL1 SP142

fit <- survfit(Surv(PFS, evtPFS)~PDL1_Sp142, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fisher.test(table(PDL1_Sp142, bestRECIST))

## CD274 RNA

fit <- survfit(Surv(PFS, evtPFS)~as.factor(CD274_RNA>5.453), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

fisher.test(table(PDL1_Sp142, bestRECIST))

## CD8 PDL1 IHC

fit <- survfit(Surv(PFS, evtPFS)~as.factor(CD8_IHC>7.676&PDL1_Sp142=="high"), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## CD8A CD274 RNA

fit <- survfit(Surv(PFS, evtPFS)~as.factor(CD8A_RNA>5.0&CD274_RNA>5.453), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

detach(myData)

## Multivariate model  ---------------------
## -----------------------------------------

PDL1_Sp142 <- PDL1_Sp142[myData$Histology!="other"]
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

## CD8 IHC model --------------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(CD8_IHC))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(CD8_IHC>7.676))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

## CD8A RNA model --------------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(CD8A_RNA))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(CD8A_RNA>5.0))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

## PDL1 Sp142 model --------------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(PDL1_Sp142))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(PDL1_Sp142=="high"))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

## CD274 RNA model --------------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(CD274_RNA))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(CD274_RNA>5.453))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

## CD8A CD274 RNA model ---------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(CD8A_RNA)|is.na(CD274_RNA))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(CD8A_RNA>5.0)+as.factor(CD274_RNA>5.453))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/pred_CD8PDL1_RNA_PFS"),
                         biomarker=val_pred,
                         outcome=NULL,
                         time=PFS[-indNA],
                         event=evtPFS[-indNA],
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)





detach(myData)
myData <- read.xlsx(file = "dataIHCRNApredictive.xlsx",
                    sheetIndex = 1, stringsAsFactors = F)
attach(myData)

## Analysis of IFNG and EIGS ---------------
## -----------------------------------------

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/pred_IFNG_RNA_PFS"),
                         biomarker=IFNG_RNA,
                         outcome=NULL,
                         time=PFS,
                         event=evtPFS,
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)

fit <- survfit(Surv(PFS, evtPFS)~as.factor(IFNG_RNA>10.57), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)


my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/pred_EIG_RNA_PFS"),
                         biomarker=EIG_RNA,
                         outcome=NULL,
                         time=PFS,
                         event=evtPFS,
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)

fit <- survfit(Surv(PFS, evtPFS)~as.factor(EIG_RNA>11.67), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

detach(myData)
myData <- subset(myData, Histology!="other")
attach(myData)

## Analysis of IFNG and EIGS ---------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(IFNG_RNA))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(IFNG_RNA>10.57))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/pred_multi_IFNG_RNA_PFS"),
                         biomarker=val_pred,
                         outcome=NULL,
                         time=PFS[-indNA],
                         event=evtPFS[-indNA],
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)


## EIG signature model --------------------

indNA <- which(is.na(PFS)|is.na(Sex)|is.na(AgeNivo)|is.na(WHO.performance.status)|is.na(Histology)|is.na(EIG_RNA))
myCox <- coxph(Surv(PFS, evtPFS)~Sex+AgeNivo+as.factor(WHO.performance.status)+Histology+as.factor(EIG_RNA>11.67))
summary(myCox)
extractAIC(myCox)
val_pred <- myCox$linear.predictors
time_dep <- survivalROC(PFS[-indNA], evtPFS[-indNA], marker=val_pred,
                        predict.time=6, method="KM")
time_dep$AUC

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/pred_multi_EIGS_RNA_PFS"),
                         biomarker=val_pred,
                         outcome=NULL,
                         time=PFS[-indNA],
                         event=evtPFS[-indNA],
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)
