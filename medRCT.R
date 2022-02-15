#################################################################################
#  medRCT.R                                                                     #
#  Version 0.1                                                                  #
#  Author: Margarita Moreno-Betancur                                            #
#  Date: 14 February 2022                                                       # 
#                                                                               #
#  Function to estimate mediation interventional effects that emulate a target  #
#  trial with any number of multiple mediators, including some not of primary   #
#  interest but that are intermediate confounders                               #
#                                                                               #
#  Reference:                                                                   #
#  Moreno-Betancur et al. "Mediation effects that emulate a target randomised   #
#  trial: Simulation-based evaluation of ill-defined interventions on multiple  #
#  mediators" Statistical Methods in Medical Research 2021. 30(6):1395-1412     #
#                                                                               #
#                                                                               #
#################################################################################

### Arguments for function:
# dat: A data.frame with the data for analysis
# ind: A vector of indices that define the sample from dat on which to conduct the analysis. Defaults to all rows of dat.
#      This facilitates the use of this function within the boot() function from the boot package
# exposure: Character string with name of the exposure in dat, which must be a binary 0/1 variable, with 1=exposed and 0=unexposed
# outcome: Character string with name of the outcome in dat, which must be a binary 0/1 variable, with 1=case and 0=non-case
# mediators: Character vector with the names of the mediators in dat, which must be binary 0/1 variables in the order required for  
#            estimation of the interventional effects accounting for flow-on effects. The index of the first mediator of interest is 
#            defined using the argument "first", the previous ones are considered intermediate confounders
# first: Index of first mediator of interest (the previous ones are treated as intermediate confounders, not of interest)  
# confounders: Character vector with names of the confounders in dat, which must be of the required class (e.g. factor if appropriate)
# Xconfsint: Character string specifying the exposure-confounder or confounder-confounder interaction terms to include in all 
#            regression models in the procedure. Defaults to include all two-way exposure-confounder interactions and no confounder-confounder 
#            interactions. All models also include all relevant two-way exposure-mediator and mediator-mediator interactions 
#            (cannot be changed through an argument to the function)
# mcsim: Number of Monte Carlo simulations to conduct     


medRCT<-function(dat, ind=1:nrow(dat),exposure, outcome, mediators, first, confounders, Xconfsint=NULL, mcsim)
{
  
  K<-length(mediators)
  
  #Rename all variables & prepare dataset
  dat$X<-dat[,exposure]
  dat$Y<-dat[,outcome]
  for(k in 1:K) dat[,paste("M",k,sep="")]<-dat[,mediators[k]]
  
  dat<-dat[,c("X",paste("M",1:K,sep=""),"Y",confounders)]
  
  #Take boostrap sample
  data<-dat[ind,]
  
  #Set flag to capure bootstrap samples to reject
  flag<-FALSE 
  
  #Prepare confounder terms for formulae (defaults to all exposure-confounder interactions if not provided)
  if(is.null(Xconfsint)) Xconfsint<-paste(paste(rep("X",length(confounders)),confounders,sep="*"),collapse="+")

  #Replicate dataset for simulations
  dat2<-data
  dat2[,1:(2+K)]<-NA_integer_
  dat2<-coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
  n<-nrow(dat2)
  
  ## ESTIMATE DISTRIBUTIONS ##
  
  #### Joint of M1 to MK under X=0 and X=1
  
  for(k in 1:K)
  {
    if(k==1)
      fit<-glm(as.formula(paste("M",k,"~X+",Xconfsint,sep="")),data=data,family=binomial) else
         fit<-glm(as.formula(paste("M",k,"~(X+",paste(paste("M",1:(k-1),sep=""),collapse="+"),")^2+",Xconfsint,sep="")),data=data,family=binomial)
  
       if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
        
        for(a in c(0,1))  
        {
          
          dat2$X<-a
          
          if(k!=1)
          {for(l in 1:(k-1))
            dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),(l-1)),rep("m",K-(l-1))),collapse=""),sep=""))}
          
          assign(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""),
                 rbinom(n,1,predict(fit,newdata=dat2,type="response")))
        }
  }
  
  
  #### For p_first,..., p_K
  
  #Marginals under X=0
  
  for(k in first:K)
  {
    
    fit<-glm(as.formula(paste("M",k,"~X+",Xconfsint,sep="")),data=data,family=binomial)
    
    if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
    
    a<-0
    dat2$X<-a
    assign(paste("m",k,"_",a,"_",paste(rep("m",K),collapse=""),sep=""),rbinom(n,1,predict(fit,newdata=dat2,type="response")))
    
  }
  
  ### Joint of others under X=1
  
  for(MM in first:K)
  {
    for(k in setdiff(first:K,MM))
    {
     
        fit<-glm(as.formula(paste("M",k,"~(X+",paste(paste("M",setdiff(1:(k-1),MM),sep=""),collapse="+"),")^2+",Xconfsint,sep="")),
                 data=data,family=binomial)
      
      if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
      
        a<-1
        dat2$X<-a
        
        if(k!=1)
        {for(l in setdiff(1:(k-1),MM))
          dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),(l-1)),rep("m",K-(l-1))),collapse=""),sep=""))}
        
        assign(paste("m",k,"_",a,"_",paste(c(rep(paste(a),min(k-1,MM-1)),"m",rep(paste(a),max(k-1-MM,0)),rep("m",K-1-min(k-1,MM-1)-max(k-1-MM,0))),collapse=""),sep=""),
               rbinom(n,1,predict(fit,newdata=dat2,type="response")))
      
    }
  }
  
  
  #### For p_first_prime,...., p_K_prime
  
  ### Conditionals under X=1
  
  for(MM in first:(K-1))
  {
    
    for(k in (MM+1):K)
    {
      
      fit<-glm(as.formula(paste("M",k,"~(X+",paste(paste("M",1:(k-1),sep=""),collapse="+"),")^2+",Xconfsint,sep="")),
               data=data,family=binomial)
      
      if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
      
      a<-1
      dat2$X<-a
      
      if(MM!=1)
      {
      for(l in 1:(MM-1))
        dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),(l-1)),rep("m",K-(l-1))),collapse=""),sep=""))
      }
      dat2[,paste("M",MM,sep="")]<-get(paste("m",MM,"_",0,"_",paste(rep("m",K),collapse=""),sep=""))
      
      if(k>(MM+1))
      {
        for(l in (MM+1):(k-1))
        dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep(paste(a),MM-1),0,rep(paste(a),max(l-1-MM,0)),rep("m",K-MM-max(l-1-MM,0))),collapse=""),sep=""))
      }
      
      assign(paste("m",k,"_",a,"_",paste(c(rep(paste(a),MM-1),0,rep(paste(a),max(k-1-MM,0)),rep("m",K-MM-max(k-1-MM,0))),collapse=""),sep=""),
             rbinom(n,1,predict(fit,newdata=dat2,type="response")))
      
    }
  }
  
  
  #### For p_all
  
  # Joint of main ones under X=0
  
  for(k in (first+1):K)
  {
     fit<-glm(as.formula(paste("M",k,"~(X+",paste(paste("M",first:(k-1),sep=""),collapse="+"),")^2+",Xconfsint,sep="")),
               data=data,family=binomial)
   
    
    if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
    
      a<-0
      dat2$X<-a
      
      for(l in first:(k-1))
      dat2[,paste("M",l,sep="")]<-get(paste("m",l,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),(l-first)),rep("m",K-l+1)),collapse=""),sep=""))
      
      assign(paste("m",k,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),k-first), rep("m",K+1-k)),collapse=""),sep=""),
             rbinom(n,1,predict(fit,newdata=dat2,type="response")))
    
    
    
  }
  
  
  
  ### outcome
  
  # Y
  
  fit<-glm(as.formula(paste("Y~(X+",paste(paste("M",1:K,sep=""),collapse="+"),")^2+",Xconfsint,sep="")),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
 
  
  ### ESTIMATE OUTCOME EXPECTATION IN EACH ARM ###
  
  #p_ctr
  
  a<-0
  dat2$X<-a
  for(k in 1:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  y0<-predict(fit,newdata=dat2,type="response")
  
  p_ctr<-mean(y0)
  
  
  #p_trt
  
  a<-1
  dat2$X<-a
  for(k in 1:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  y1<-predict(fit,newdata=dat2,type="response")
  
  p_trt<-mean(y1)  
  
  
  #p_all
  
  a<-1
  dat2$X<-a
  for(k in 1:(first-1))
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  a<-0
  for(k in first:K)
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep("m",first-1),rep(paste(a),k-first),rep("m",K+1-k)),collapse=""),sep=""))
  }
  
  y1<-predict(fit,newdata=dat2,type="response")
  
  p_all<-mean(y1)
  
  
  #p_first....p_K
  
  a<-1
  dat2$X<-a
  
  for(k in 1:(first-1))
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  
  for(MM in first:K)
  {
    
    dat2[,paste("M",MM,sep="")]<-get(paste("m",MM,"_",0,"_",paste(paste(rep("m",K),sep=""),collapse=""),sep=""))

    for(k in setdiff(first:K,MM))
    {
      dat2[,paste("M",k,sep="")]<- get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),min(k-1,MM-1)),"m",rep(paste(a),max(k-1-MM,0)),rep("m",K-1-min(k-1,MM-1)-max(k-1-MM,0))),collapse=""),sep=""))
    }
    
    y0<-predict(fit,newdata=dat2,type="response")
    
    assign(paste("p_",MM,sep=""),mean(y0))
    
  }
  
  #p_first_prime....p_Kminus1_prime
  
  a<-1
  dat2$X<-a
  
  for(MM in first:(K-1))
  {
  
  if(MM!=1)
  {  
  for(k in 1:(MM-1))
  {
    dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),(k-1)),rep("m",K-(k-1))),collapse=""),sep=""))
  }
  }
  
  dat2[,paste("M",MM,sep="")]<-get(paste("m",MM,"_",0,"_",paste(paste(rep("m",K),sep=""),collapse=""),sep=""))
    
  if((MM+1)<=K)
  {
    for(k in (MM+1):K)
  {  
  dat2[,paste("M",k,sep="")]<-get(paste("m",k,"_",a,"_",paste(c(rep(paste(a),MM-1),0,rep(paste(a),max(k-1-MM,0)),rep("m",K-MM-max(k-1-MM,0))),collapse=""),sep=""))
  }
  } 
  
  y0<-predict(fit,newdata=dat2,type="response")
    
  assign(paste("p_",MM,"_prime",sep=""),mean(y0))
    
  }

  
  ### ESTIMATE EFFECTS ###
  
  ### Interventional effects 
  
  for(k in first:K)
    assign(paste("IIE_",k,sep=""),p_trt-get(paste("p_",k,sep="")))
  
  for(k in first:(K-1))
    assign(paste("IIE_",k,"_prime",sep=""),p_trt-get(paste("p_",k,"_prime",sep="")))
  
  
  IIE_all<-p_trt-p_all

  ### Collect and return results
  
  res<-vector()
  for(k in first:K) res<-c(res,get(paste("IIE_",k,sep="")))
  for(k in first:(K-1)) res<-c(res,get(paste("IIE_",k,"_prime",sep="")))
  res<-c(res,IIE_all, p_trt, p_ctr,p_all)
  for(k in first:K) res<-c(res, get(paste("p_",k,sep="")))
  for(k in first:(K-1)) res<-c(res, get(paste("p_",k,"_prime",sep="")))
  
  names(res)<-c(paste("IIE_",first:K,sep=""),paste("IIE_",first:(K-1),"_prime",sep=""), "IIE_all", 
                "p_trt", "p_ctr","p_all", paste("p_",first:K,sep=""), paste("p_",first:(K-1),"_prime",sep=""))
  
    if(!flag)return(res) else
    return(rep(NA,length(res)))
}