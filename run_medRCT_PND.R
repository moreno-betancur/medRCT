#################################################################################
#  run_medRCT_PND.R                                                             # 
#                                                                               #
#  Author: Margarita Moreno-Betancur                                            #
#  Date: 6 April 2020                                                           #
#                                                                               #
#  Code running function medRCT_PND.R on imputed data using bootstrap for       #
#  the paper:                                                                   #
#                                                                               #
#    Spry et al. "Preventing perinatal depression : a causal mediation analysis #
#                 of a 20-year preconception cohort"                            # 
#                                                                               #
#                                                                               #
#################################################################################



tt<-proc.time()
args<-c("m","Y_8wk","T4") # arguments: imputation number, outcome ("Y_8wk" or "Y_1yr"), table ("T4" "T5" or "S3")
m<-as.numeric(args[1])
outcome<-args[2]     
scenario<-args[3]   


seedn<-ifelse(outcome=="Y_8wk", 93847, 45263)

library(boot)
library(zoo)
options(warn=0)

source("medRCT_PND.R")


set.seed(seedn*m)
datMI<-read.csv("datPND.csv",row.names=NULL)

if(scenario=="T4") meds<-c("L1","L2","M1_10","M2_10","M3_10","M4_un") else
  if(scenario=="T5") meds<-c("L1","L2","M1a","M1b","M2_10","M3_10","M4_un") else
    if(scenario=="S3") meds<-c("L1","L2","M1_25","M2_25","M3_25","M4_un")

confsL<-
  c("parents_divorce", "parents_noyr12", "cauc_3gen",             ### Confounders - Origin
    "smkd_any", "binge_any", "thcw_any", "anti_any", "edu_noyr12",  ### /// Confounders - Adolescent
    "div_any", "livebirth_any", "miscar_any", "termin_any")         ### /// Confounders - Adult

datMI$Y<-datMI[,outcome]

for(s in 1:(length(meds)))
  datMI[,paste("M",s,sep="")]<-datMI[,meds[s]]

S<-length(meds)

dat<-datMI[datMI$"imp"==m,][,c("A",paste("M",1:S,sep=""),"Y",confsL)]


bstrap<-boot(data=dat, statistic=medRCT_PND, 
             exposure="A",outcome="Y",mediators=paste("M",1:S,sep=""),confounders=confsL,mcsim=200, RCT="one-policy",
             K=S,first=3,redconf=3,whichred=S,
             stype="i", R=1000,parallel="multicore",ncpus=5)

RES<-data.frame(EST=bstrap$t0,SE=apply(bstrap$t,2,sd,na.rm=T))
write.csv(RES,paste("RES_",m,"_",outcome,"_",scenario,".csv",sep=""))
write.csv(bstrap$t,paste("DIST_",m,"_",outcome,"_",scenario,".csv",sep=""))

png(paste("bs_",m,"_",outcome,"_",scenario,".png",sep=""))
par(mfrow=c(3,3))
for(i in 1:7) hist(bstrap$t[,i])
dev.off()


print((proc.time()-tt)/60)
