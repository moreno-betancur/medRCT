---
title: "`medRCT`: Estimating mediation effects that emulate a target randomized controlled trial (RCT)"
author: "Margarita Moreno-Betancur"
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{"medRCT"}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
  
This vignette provides a worked example on simulated data showing how to use the R function `medRCT_4med` provided in [this](https://github.com/moreno-betancur/medRCT) repository for the paper (click [here](https://raw.githack.com/moreno-betancur/medRCT/master/medRCT_4med.R) to download it): 

>Moreno-Betancur M, Moran P, Becker D, Patton GC, Carlin JB. "Defining mediation effects for multiple mediators using the concept of the target randomized trial". [https://arxiv.org/abs/1907.06734](https://arxiv.org/abs/1907.06734)



## Loading and looking at the example dataset

To illustrate, we use a dataset `dat4med.csv` (available for download [here](https://raw.githack.com/moreno-betancur/medRCT/master/dat4med.csv)) that was simulated roughly based on one imputed version of the real data from the Victorian Adolescent Health Cohort Study (VAHCS) that was analysed in the paper.
  
The dataset consists of the following variables, none of which has missing data:

* ID: identifier variable 
* C1,...,C6: set of binary confounders
* A: binary exposure
* M1,...,M4: set of binary mediators
* Y: binary outcome


```r
# Load and look
dat <- read.csv("dat4med.csv")
head(dat)
```

```
##   ID C1 C2 C3 C4 C5 C6 A M1 M2 M3 M4 Y
## 1  1  1  1  1  0  1  1 1  0  0  1  0 1
## 2  2  0  0  0  0  0  0 0  0  0  1  0 0
## 3  3  0  0  0  0  1  0 0  0  0  0  0 0
## 4  4  0  0  0  1  0  0 0  0  0  1  0 0
## 5  5  1  1  0  1  0  1 1  0  1  1  0 0
## 6  6  0  0  1  0  0  0 0  0  0  1  1 0
```

```r
# Summarise
apply(dat[, 2:ncol(dat)], 2, mean, na.rm = F)
```

```
##     C1     C2     C3     C4     C5     C6      A     M1     M2     M3 
## 0.3560 0.1370 0.1800 0.2470 0.5115 0.3425 0.0985 0.2235 0.1270 0.6615 
##     M4      Y 
## 0.1095 0.2495
```

## Loading function and required libraries

The function `medRCT_4med` is provided in [this](https://github.com/moreno-betancur/medRCT) repository to conduct analyses as in the paper (click [here](https://raw.githack.com/moreno-betancur/medRCT/master/medRCT_4med.R) to download it). 

While there is a plan to make the function more general in the future, this version of the function assumes a setting with 4 binary interdependent mediators, with exposure and outcome being binary as well. During the estimation process, the function includes all 2-way interactions amongst exposure and mediators in the parametric models that are the building pieces for the simulation-based g-computation estimation procedure. The package `zoo` is required by the function and it is set up so that it can be called using the `boot` function from the `boot` package, which also needs to be loaded to run the analysis using the bootstrap. 


```r
library(boot)
library(zoo)
source("medRCT_4med.R")
```

The function `medRCT_4med` takes the following arguments:

* `dat`: Dataset (data.frame)
* `ind`: Indices of records in `dat` to conduct the analysis on (numerical vector), which is required by `boot`; default is to use all of `dat`, as is
* `exposure`: Exposure name (character vector, length 1)
* `outcome`: Outcome name (character vector, length 1)
* `mediators`: Names of mediators (character vector, length 4), in the desired order for reporting; this order also defines the order of the sequence for effects under sequential policies
* `confounders`: Names of the confounders in the dataset (Character vector, length>0)
* `sim`: Number of desired Monte Carlo simulation runs; default is to use 200 runs
* `RCT`: Type of RCT to emulate (charactervector, length 1); must be either "one-policy" to estimate mediation effects emulating an RCT under the one-policy premise,  "sequential" to estimate mediation effects emulating an RCT under sequential policies, or "both" to obtain both sets of estimates


## Using the function

### Estimate mediation effects under "one-policy premise"

The following code shows how to obtain estimates of mediation effects that emulate an RCT under a "one-policy premise"

```r
#set seed for reproducibility and also for estimates of TCE, IDE and IIE_1/IIE_seq1 obtained with the two RCT options to coincide
set.seed(4750) 

# Estimate the effects with the bootrstrap
bstrap<-boot(data=dat, statistic=medRCT_4med, 
             exposure="A", outcome="Y", mediators=c("M1","M2","M3","M4"),
             confounders=c("C1","C2","C3","C4","C5","C6"), mcsim=200, RCT="one-policy",
             stype="i", R=2)

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE
RES$pvalue<-2*pnorm(abs(RES$Estimate/RES$SE),lower.tail = F)
RES$PropTCE<-round(100*(RES$Estimate/RES$Estimate[1]),0)
RES<-round(RES[,-2],2)
RES<-cbind(Estimand=c("TCE","IDE","IIE_1","IIE_2","IIE_3","IIE_4","IIE_int"),RES)

print(RES,row.names = F)
```

```
##  Estimand Estimate CIlow CIupp pvalue PropTCE
##       TCE     0.12  0.05  0.20   0.00     100
##       IDE     0.07  0.02  0.12   0.00      56
##     IIE_1     0.01  0.00  0.01   0.04       4
##     IIE_2     0.01  0.00  0.03   0.03      11
##     IIE_3     0.04  0.03  0.05   0.00      33
##     IIE_4     0.01  0.00  0.02   0.08       7
##   IIE_int    -0.01 -0.02 -0.01   0.00     -11
```

```r
# Can check bootstrap distribution is approx normal using (not run here)
# par(mfrow=c(3,3))
# for(i in 1:7) hist(bstrap$t[,i])
```

To learn more about the interpretation of these and the above results refer to the paper (reference below).

### Estimate mediation effects under sequential policies

The following code shows how to obtain estimates of mediation effects that emulate an RCT under sequential policies

```r
#set seed for reproducibility and also for estimates of TCE, IDE and IIE_1/IIE_seq1 obtained with the two RCT options to coincide
set.seed(4750) 

# Estimate the effects with the bootrstrap
bstrap<-boot(data=dat, statistic=medRCT_4med, 
             exposure="A", outcome="Y", mediators=c("M1","M2","M3","M4"),
             confounders=c("C1","C2","C3","C4","C5","C6"), mcsim=200, RCT="sequential",
             stype="i", R=2)

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE
RES$pvalue<-2*pnorm(abs(RES$Estimate/RES$SE),lower.tail = F)
RES$PropTCE<-round(100*(RES$Estimate/RES$Estimate[1]),0)
RES<-round(RES[,-2],2)
RES<-cbind(Estimand=c("TCE","IDE","IIE_seqfull","IIE_seq1","IIE_seq2","IIE_seq3","IIE_seq4","IIE_seqint"),RES)

print(RES,row.names = F)
```

```
##     Estimand Estimate CIlow CIupp pvalue PropTCE
##          TCE     0.12  0.05  0.20   0.00     100
##          IDE     0.07  0.02  0.12   0.00      56
##  IIE_seqfull     0.06  0.03  0.08   0.00      46
##     IIE_seq1     0.01  0.00  0.01   0.04       4
##     IIE_seq2     0.01  0.00  0.03   0.10      10
##     IIE_seq3     0.04  0.03  0.05   0.00      30
##     IIE_seq4     0.00  0.00  0.00   0.08       2
##   IIE_seqint     0.00  0.00  0.00   0.00      -2
```

```r
# Can check bootstrap distribution is approx normal using (not run here)
# par(mfrow=c(3,3))
# for(i in 1:7) hist(bstrap$t[,i])
```



## Cite this as:

Moreno-Betancur M, Moran P, Becker D, Patton GC, Carlin JB. "Defining mediation effects for multiple mediators using the concept of the target randomized trial". [https://arxiv.org/abs/1907.06734](https://arxiv.org/abs/1907.06734)


