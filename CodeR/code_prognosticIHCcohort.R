#############################################
## --------------------------------------- ##
## -- Analysis of IHC prognostic cohort -- ##
## --------------------------------------- ##
#############################################


## Set path to working directory

path.wd <- "my.working.directory"
setwd(path.wd)

## Load required packages

library(xlsx)
library(survminer)
source("CodeR/code_cutoff.R")

## Import data

myData <- read.xlsx(file = "Data/dataIHCprognostic.xlsx",
                    sheetIndex = 1)
attach(myData)


## ------------------------------------------
## Analysis of clinical data ----------------
## ------------------------------------------


## Sex

table(Sex)
prop.table(table(Sex))

## Age

summary(AgeDGM)
sd(AgeDGM)

## WHO performance status

table(WHO.performance.status)
prop.table(table(WHO.performance.status))

## Histology

table(Histology)
prop.table(table(Histology))

## Tumor stage

table(Tumor.stage)
prop.table(table(Tumor.stage))

## Smoker

table(Smoker)
prop.table(table(Smoker))


## ------------------------------------------
## Analysis of IHC data ---------------------
## ------------------------------------------


## 22C3 and Sp142 labeling concordance

PDL1_Sp142 <- factor(ifelse(PDL1_IC_Sp142<2&PDL1_TC_Sp142<2,"low","high"), levels = c("low", "high"))
PDL1_22C3 <- factor(PDL1_22C3, levels = c("<1%", "1-49%", ">50%"))
table(PDL1_Sp142, PDL1_22C3)
fisher.test(table(PDL1_Sp142, PDL1_22C3))

## Correlation between CD8 and PDL1 Sp142

wilcox.test(CD8~PDL1_Sp142)


## Association with OS ----------------------

## CD8

my_cutoff_tfh=get.cutoff(type="survival_significance",
                         paste0("cutoffFinder/CD8_IHC_OS"),
                         biomarker=CD8,
                         outcome=NULL,
                         time=OS,
                         event=evtOS,
                         cutoff=NULL,
                         threshold=NULL,
                         plots=c("kaplanmeier"),
                         nmin=10)
fit <- survfit(Surv(OS, evtOS)~as.factor(CD8>7.807), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## PDL1 SP142

fit <- survfit(Surv(OS, evtOS)~PDL1_Sp142, data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## CD8+ PDL1 SP142+

fit <- survfit(Surv(OS, evtOS)~as.factor(CD8<7.807&PDL1_Sp142=="high"), data = myData)
ggsurvplot(fit, surv.median.line = "hv", pval = TRUE)

## Multivariate model  ---------------------

summary(coxph(Surv(OS, evtOS)~Sex+AgeDGM+as.factor(WHO.performance.status)+Histology+as.factor(!(CD8<7.807&PDL1_Sp142=="high"))))

detach(myData)
