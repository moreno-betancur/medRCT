#################################################################################
#  medRCT_PND.R                                                                 # 
#                                                                               #
#  Author: Margarita Moreno-Betancur                                            #
#  Date: 6 April 2020                                                           #
#                                                                               #
#  Function to estimate mediation effects that emulate a target RCT             #
#  under the "one-policy premise" as per the paper:                             #
#                                                                               #
#    Moreno-Betancur et al. "Mediation effects that emulate a target randomised #
#                          trial: Simulation-based evaluation of ill-defined    #
#                          interventions on multiple mediators"                 # 
#                         (https://arxiv.org/abs/1907.06734)                    #
#                                                                               #
#  Version: Adaptation of original medRCT_4med function for perinatal           #
#           depression (PND) paper:                                            #
#                                                                               #
#     Spry et al. "Preventing perinatal depression : a causal mediation analysis#
#                  of a 20-year preconception cohort"                           # 
#                                                                               #
#  Differences with medRCT_4med:                                                #
#  - Number of mediators, including some only used as post-exposure confounders #
#  - Models without interactions and, for one mediator, excluding a confounder  #
#  - No sequential policies                                                     #  
#  - Three-level exposure                                                       #
#                                                                               #
#################################################################################


medRCT_PND<-function(dat, ind=1:nrow(dat),exposure, outcome, mediators, confounders, mcsim=200, RCT="one-policy",K=6,first=3,redconf=3,whichred=K)
{

  #Rename all variables & prepare dataset
  dat$A<-dat[,exposure]
  dat$Y<-dat[,outcome]
  
  for(k in 1:K) dat[,paste("M",k,sep="")]<-dat[,mediators[k]]

  dat<-dat[,c("A",paste("M",1:K,sep=""),"Y",confounders)]
  
  #Take boostrap sample
  data<-dat[ind,]
  
  #Set flag to capure bootstrap samples to reject
  flag<-FALSE 
  
  #Prepare confounders for formulae
  confs<-paste(confounders,collapse="+")
  confs2<-paste(confounders[-redconf],collapse="+")
  
  #Replicate dataset for simulations
  dat2<-data
  dat2[,1:(2+K)]<-NA_integer_
  dat2<-coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
  n<-nrow(dat2)
  
  ## ESTIMATE DISTRIBUTIONS ##
  
  #### Joint of M1 to MK
  
  for(k in 1:K)
  {
    if(k==1)
      fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",confs,sep="")),data=data,family=binomial) else
        if(k!=whichred)
        fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",1:(k-1),sep=""),collapse="+"),"+",confs,sep="")),
                 data=data,family=binomial)
        else    fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",1:(k-1),sep=""),collapse="+"),"+",confs2,sep="")),
                         data=data,family=binomial)
        
      #print(summary(fit)  )
      if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
      
      for(a in c(0,1))
      {
        
        dat2$A<-a
        
        if(k!=1)
        {for(l in 1:(k-1))
          dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),(l-1)),rep("m",K-(l-1))),collapse=""),sep=""))}
        
        assign(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""),
               rbinom(n,1,predict(fit,newdata=dat2,type="response")))
      }
  }
  
  
  #### For p3, p4, p5, p6
  
  #Marginals under A=0
  
  for(k in first:K)
  {
    if(k!=whichred)
          fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",confs,sep="")),data=data,family=binomial)
    else       fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",confs2,sep="")),data=data,family=binomial)  #confs2 for last med
    
    if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
    
    for(a in c(0))
    {dat2$A<-a
    assign(paste("m",k,"_",a,"_",paste(rep("m",K),collapse=""),sep=""),rbinom(n,1,predict(fit,newdata=dat2,type="response")))
    }
  }
  
  ### Joint of others under A=1
  
  for(MM in first:K)
  {
  for(k in setdiff(first:K,MM))
  {
    if(k!=whichred)
    fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",setdiff(1:(k-1),MM),sep=""),collapse="+"),"+",confs,sep="")),
             data=data,family=binomial)
    else
      fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",setdiff(1:(k-1),MM),sep=""),collapse="+"),"+",confs2,sep="")),
               data=data,family=binomial)
    
    if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
    
    for(a in 1)
    {
      
      dat2$A<-a
      
      if(k!=1)
      {for(l in setdiff(1:(k-1),MM))
        dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),(l-1)),rep("m",K-(l-1))),collapse=""),sep=""))}
      
        assign(paste("m",k,"_",a,"_",paste(c(rep(paste(a),min(k-1,MM-1)),"m",rep(paste(a),max(k-1-MM,0)),rep("m",K-1-min(k-1,MM-1)-max(k-1-MM,0))),collapse=""),sep=""),
               rbinom(n,1,predict(fit,newdata=dat2,type="response")))

    }
  }
  }
  
  #### For pall
  
  # Joint of main ones under A=0
  
  for(k in (first+1):K)
  {
    if(k!=whichred)
      fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",first:(k-1),sep=""),collapse="+"),"+",confs,sep="")),
               data=data,family=binomial)
    else
      fit<-glm(as.formula(paste("M",k,"~as.factor(A)+",paste(paste("M",first:(k-1),sep=""),collapse="+"),"+",confs2,sep="")),
               data=data,family=binomial)
    
      
               if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
               
               for(a in 0)
               {
                 
                 dat2$A<-a
                 
                 for(l in first:(k-1))
                   dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),(l-first)),rep("m",K-l+1)),collapse=""),sep=""))
                 
                   assign(paste("m",k,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),k-first), rep("m",K+1-k)),collapse=""),sep=""),
                          rbinom(n,1,predict(fit,newdata=dat2,type="response")))
                 }
               
               
  }
  
  
  
  ### outcome
  
  # Y
  
  fit<-glm(as.formula(paste("Y~as.factor(A)+",paste(paste("M",1:K,sep=""),collapse="+"),"+",confs,sep="")),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  #print(summary(fit)  )
  
  ### ESTIMATE OUTCOME EXPECTATION IN EACH ARM ###
  
  #p_ctr
  
  a<-0
  dat2$A<-a
  for(k in 1:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
    #print(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_ctr<-mean(y0)
  
  
  #p_trt
  
  a<-1
  dat2$A<-a
  for(k in 1:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
    #print(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_trt<-mean(y1)  
  
  
  #p_all
  
  a<-1
  dat2$A<-a
  for(k in 1:(first-1))
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
    #print(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  a<-0
  for(k in first:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),k-first),rep("m",K+1-k)),collapse=""),sep=""))
   # print(paste("m",k,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),k-first),rep("m",K+1-k)),collapse=""),sep=""))
    }
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_all<-mean(y1)
  
  
  #p_first....p_K
  
  a<-1
  dat2$A<-a
  
    for(k in 1:(first-1))
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
   # print(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  for(MM in first:K)
  {
    
   dat2[,paste("M",MM,sep="")]<-get(paste("m",MM,"_",0,"_",paste(paste(rep("m",K),sep=""),collapse=""),sep=""))
   #print(paste("m",MM,"_",0,"_",paste(paste(rep("m",K),sep=""),collapse=""),sep="")) 

  for(k in setdiff(first:K,MM))
  {

    dat2[,paste("M",k,sep="")]<- get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),min(k-1,MM-1)),"m",rep(paste(a),max(k-1-MM,0)),rep("m",K-1-min(k-1,MM-1)-max(k-1-MM,0))),collapse=""),sep=""))
    #print(paste("m",k,"_",a,"_",paste(c(rep(paste(a),min(k-1,MM-1)),"m",rep(paste(a),max(k-1-MM,0)),rep("m",K-1-min(k-1,MM-1)-max(k-1-MM,0))),collapse=""),sep=""))
  }
  
   y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
   
   assign(paste("p_",MM,sep=""),mean(y0))
  
  }
  ### ESTIMATE EFFECTS ###
  
  ### Effects under "one-policy premise"
  
   for(k in first:K)
  assign(paste("IIE_",k,sep=""),p_trt-get(paste("p_",k,sep="")))
 
  IIE_int<-p_trt-p_all
  for(k in first:K) IIE_int<-IIE_int-get(paste("IIE_",k,sep=""))
  
  ### Interventional direct effect and total causal effect
  
  IDE<-p_all-p_ctr  
  TCE<-p_trt-p_ctr
  
  ### Collect and return results
  
  if(RCT=="one-policy")
  {res<-c(TCE,IDE)
   for(k in first:K) res<-c(res,get(paste("IIE_",k,sep="")))
   res<-c(res,IIE_int)} else   stop("RCT argument must be 'one-policy'")
             
  if(!flag)return(res) else
    return(rep(NA,length(res)))
}