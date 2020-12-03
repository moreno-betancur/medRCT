---
title: "`medRCT`: Estimating mediation effects that emulate a target randomized controlled trial (RCT)"
author: "Margarita Moreno-Betancur"
date: December 3, 2020
output: 
  rmarkdown::html_vignette:
    keep_md: true
vignette: >
  %\VignetteIndexEntry{"medRCT"}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
  
This vignette provides a worked example on simulated data showing how to use the R function `medRCT_4med` downloadable [here](https://raw.githack.com/moreno-betancur/medRCT/master/medRCT_4med.R) and stored [here](https://github.com/moreno-betancur/medRCT), which accompanies the paper: 

>Moreno-Betancur M, Moran P, Becker D, Patton GC, Carlin JB. "Mediation effects that emulate a target randomised trial: Simulation-based evaluation of ill-defined interventions on multiple mediators". [https://arxiv.org/abs/1907.06734](https://arxiv.org/abs/1907.06734)



## Loading and looking at the example dataset

To illustrate how to use the function, we use a dataset `dat4med.csv` (downloadable [here](https://raw.githack.com/moreno-betancur/medRCT/master/dat4med.csv)) that was simulated roughly based on one imputed version of the real data that was analysed in the paper, from the Victorian Adolescent Health Cohort Study (VAHCS).
  
The dataset consists of the following variables, none of which has missing data:

* ID: identifier variable 
* C1,...,C6: set of binary confounders
* A: binary exposure
* M1,...,M4: set of binary mediators
* Y: binary outcome


```r
#Load and look
dat<-read.csv("dat4med.csv")
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
#Summarise
apply(dat[,2:ncol(dat)],2,mean,na.rm=F)
```

```
##     C1     C2     C3     C4     C5     C6      A     M1     M2     M3     M4 
## 0.3560 0.1370 0.1800 0.2470 0.5115 0.3425 0.0985 0.2235 0.1270 0.6615 0.1095 
##      Y 
## 0.2495
```

## Loading the function and required libraries

The function `medRCT_4med` is dowloadable [here](https://raw.githack.com/moreno-betancur/medRCT/master/medRCT_4med.R) and can be used to conduct analyses as in [the paper](https://arxiv.org/abs/1907.06734). 

While the plan is to make the function more general in the future (watch [this](https://github.com/moreno-betancur/medRCT) space!), the version of the function available for now assumes a setting with 4 binary interdependent mediators, with exposure and outcome being binary as well. The function depends on the R package `zoo` and it is set up so that it can be called using the `boot` function from the `boot` package, which also needs to be loaded to run the analysis using the bootstrap: 


```r
library(boot)
library(zoo)
source("medRCT_4med.R")
```

The function `medRCT_4med` takes the following arguments:

* `dat`: Dataset (data.frame)
* `ind`: Indices of records in `dat` on which to conduct the analysis (numerical vector), which is required by `boot`; default is to use `dat` as is
* `exposure`: Exposure name (character vector, length 1)
* `outcome`: Outcome name (character vector, length 1)
* `mediators`: Names of mediators (character vector, length 4), in the desired order for reporting; this order also defines the order of the sequence for effects under sequential policies
* `confounders`: Names of the confounders in the dataset (character vector, length>0)
* `sim`: Number of desired Monte Carlo simulation runs; default is to use 200 runs
* `RCT`: Type of RCT to emulate (character vector, length 1); must be either "one-policy_A" or "one-policy_B" to estimate mediation effects emulating an RCT under the one-policy premise using approches (a) and (b) in the paper, respectively;  "sequential" to estimate mediation effects emulating an RCT under sequential policies; or "all" to obtain all sets of estimates

Of note, during the estimation process, the function includes all 2-way interactions amongst exposure and mediators in the parametric models that are the building pieces for the simulation-based g-computation estimation procedure. 

## Using the function `medRCT_4med` 

### Estimate mediation effects under "one-policy premise"

The following code shows how to obtain estimates of mediation effects that emulate an RCT under a "one-policy premise", using approach (a) and 100 bootstrap runs (for illustrative purposes - might want to use more in practice)

```r
# Set seed for reproducibility but also so that estimates of common effects obtained using either of the RCT options coincide
set.seed(4750) 

# Estimate the effects with the bootrstrap
bstrap<-boot(data=dat, statistic=medRCT_4med, 
             exposure="A", outcome="Y", mediators=c("M1","M2","M3","M4"),
             confounders=c("C1","C2","C3","C4","C5","C6"), mcsim=200, RCT="one-policy_A",
             stype="i", R=100)

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE
RES$pvalue<-2*pnorm(abs(RES$Estimate/RES$SE),lower.tail = F)
RES$PropTCE<-round(100*(RES$Estimate/RES$Estimate[1]),0) #Express effects as a proportion of the TCE
RES<-round(RES[,-2],2)
RES<-cbind(Estimand=c("TCE","IDE","IIE_1","IIE_2","IIE_3","IIE_4","IIE_int"),RES)

print(RES,row.names = F)
```

```
##  Estimand Estimate CIlow CIupp pvalue PropTCE
##       TCE     0.12  0.06  0.18   0.00     100
##       IDE     0.07  0.01  0.13   0.02      56
##     IIE_1     0.01 -0.01  0.02   0.31       5
##     IIE_2     0.01  0.00  0.03   0.07      11
##     IIE_3     0.04  0.01  0.07   0.00      32
##     IIE_4     0.01 -0.01  0.03   0.38       7
##   IIE_int    -0.01 -0.03  0.00   0.05     -11
```

```r
# Can check bootstrap distribution is approx normal using (not run here)
# par(mfrow=c(3,3))
# for(i in 1:7) hist(bstrap$t[,i])
```

To obtain estimates of mediation effects under a "one-policy premise" using approach (b), one can use a similar code:


```r
# Set seed for reproducibility but also so that estimates of common effects obtained using either of the RCT options coincide
set.seed(4750) 

# Estimate the effects with the bootrstrap
bstrap<-boot(data=dat, statistic=medRCT_4med, 
             exposure="A", outcome="Y", mediators=c("M1","M2","M3","M4"),
             confounders=c("C1","C2","C3","C4","C5","C6"), mcsim=200, RCT="one-policy_B",
             stype="i", R=100)

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE
RES$pvalue<-2*pnorm(abs(RES$Estimate/RES$SE),lower.tail = F)
RES$PropTCE<-round(100*(RES$Estimate/RES$Estimate[1]),0) #Express effects as a proportion of the TCE
RES<-round(RES[,-2],2)
RES<-cbind(Estimand=c("TCE","IDE","IIE_1_prime","IIE_2_prime","IIE_3_prime","IIE_4_prime","IIE_int_prime"),RES)

print(RES,row.names = F)
```

```
##       Estimand Estimate CIlow CIupp pvalue PropTCE
##            TCE     0.12  0.06  0.18   0.00     100
##            IDE     0.07  0.01  0.13   0.02      56
##    IIE_1_prime     0.00 -0.01  0.01   0.71       1
##    IIE_2_prime     0.01  0.00  0.03   0.17       9
##    IIE_3_prime     0.04  0.01  0.06   0.01      30
##    IIE_4_prime     0.01 -0.01  0.03   0.38       7
##  IIE_int_prime     0.00 -0.01  0.01   0.46      -3
```

```r
# Can check bootstrap distribution is approx normal using (not run here)
# par(mfrow=c(3,3))
# for(i in 1:7) hist(bstrap$t[,i])
```

To learn more about the interpretation of these results and how the two approaches above compare refer to the paper (reference below).

### Estimate mediation effects under sequential policies

The following code shows how to obtain estimates of mediation effects that emulate an RCT under sequential policies

```r
# Set seed for reproducibility but also so that estimates of common effects obtained using either of the RCT options coincide
set.seed(4750) 

# Estimate the effects with the bootrstrap
bstrap<-boot(data=dat, statistic=medRCT_4med, 
             exposure="A", outcome="Y", mediators=c("M1","M2","M3","M4"),
             confounders=c("C1","C2","C3","C4","C5","C6"), mcsim=200, RCT="sequential",
             stype="i", R=100)

# Set-up results table
RES<-data.frame(Estimate=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
RES$CIlow<-RES$Estimate-1.96*RES$SE
RES$CIupp<-RES$Estimate+1.96*RES$SE
RES$pvalue<-2*pnorm(abs(RES$Estimate/RES$SE),lower.tail = F)
RES$PropTCE<-round(100*(RES$Estimate/RES$Estimate[1]),0) #Express effects as a proportion of the TCE
RES<-round(RES[,-2],2)
RES<-cbind(Estimand=c("TCE","IDE","IIE_seqfull","IIE_seq1","IIE_seq2","IIE_seq3","IIE_seq4","IIE_seqint"),RES)

print(RES,row.names = F)
```

```
##     Estimand Estimate CIlow CIupp pvalue PropTCE
##          TCE     0.12  0.06  0.18   0.00     100
##          IDE     0.07  0.01  0.13   0.02      56
##  IIE_seqfull     0.06  0.02  0.09   0.00      45
##     IIE_seq1     0.01 -0.01  0.02   0.31       5
##     IIE_seq2     0.01  0.00  0.02   0.09       9
##     IIE_seq3     0.04  0.01  0.06   0.00      31
##     IIE_seq4     0.00 -0.02  0.02   0.98       0
##   IIE_seqint     0.00 -0.01  0.00   0.30      -2
```

```r
# Can check bootstrap distribution is approx normal using (not run here)
# par(mfrow=c(3,3))
# for(i in 1:7) hist(bstrap$t[,i])
```



## Cite this as:

Moreno-Betancur M, Moran P, Becker D, Patton GC, Carlin JB. "Mediation effects that emulate a target randomised trial: Simulation-based evaluation of ill-defined interventions on multiple mediators". [https://arxiv.org/abs/1907.06734](https://arxiv.org/abs/1907.06734)


