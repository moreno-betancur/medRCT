#################################################################################
#  medRCT_4med.R                                                                # 
#                                                                               #
#  Function to estimate mediation effects that emulate a target RCT             #
#  with 4 mediators as per the paper:                                           #
#                                                                               #
#    Moreno-Betancur et al."Mediation effects that emulate a target randomised  #
#                          trial: Simulation-based evaluation of ill-defined    #
#                          interventions on multiple mediators"                 # 
#                          (https://arxiv.org/abs/1907.06734)                   #
#                                                                               #
#                                                                               #
#  Margarita Moreno-Betancur, 3 Dec 2020                                        #
#                                                                               #
#################################################################################


medRCT_4med<-function(dat, ind=1:nrow(dat),exposure, outcome, mediators, confounders, mcsim=200, RCT)
{
  K<-4 
  
  #Rename all variables & prepare dataset
  dat$A<-dat[,exposure]
  dat$Y<-dat[,outcome]
  dat$M1<-dat[,mediators[1]]
  dat$M2<-dat[,mediators[2]]
  dat$M3<-dat[,mediators[3]]
  dat$M4<-dat[,mediators[4]]
  
  dat<-dat[,c("A","M1","M2","M3","M4","Y",confounders)]
  
  #Take boostrap sample
  data<-dat[ind,]
  
  #Set flag to capure bootstrap samples to reject
  flag<-FALSE 
  
  #Prepare confounders for formulae
  confs<-paste(confounders,collapse="+")
  
  #Replicate dataset for simulations
  dat2<-data
  dat2[,1:(2+K)]<-NA_integer_
  dat2<-coredata(dat2)[rep(seq(nrow(dat2)),mcsim),]
  n<-nrow(dat2)
  
  ## ESTIMATE DISTRIBUTIONS ##
  
  ### marginals
  
  #M1 - MARGINAL
  
  fit<-glm(as.formula(paste("M1~A+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  m1_0_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  m1_1_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  #M2 - MARGINAL
  
  fit<-glm(as.formula(paste("M2~A+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  m2_0_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  m2_1_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  #M3 - MARGINAL
  
  fit<-glm(as.formula(paste("M3~A+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  m3_0_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  m3_1_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  #M4 - MARGINAL
  
  fit<-glm(as.formula(paste("M4~A+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  m4_0_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  m4_1_mmmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  ### order 2 conditionals
  
  #M2 - GIVEN M1
  
  fit<-glm(as.formula(paste("M2~A*M1+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  m2_0_0mmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  m2_1_1mmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  m2_1_0mmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  #M3 - GIVEN M1
  
  fit<-glm(as.formula(paste("M3~A*M1+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  m3_0_0mmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  m3_1_1mmm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  #M3 - GIVEN M2
  
  fit<-glm(as.formula(paste("M3~A*M2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M2<-m2_0_mmmm
  m3_0_m0mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M2<-m2_1_mmmm
  m3_1_m1mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  #M4 - GIVEN M3
  
  fit<-glm(as.formula(paste("M4~A*M3+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M3<-m3_0_mmmm
  m4_0_mm0m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M3<-m3_1_mmmm
  m4_1_mm1m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  ### order 3 conditionals
  
  #M3 - GIVEN M1 and M2
  
  fit<-glm(as.formula(paste("M3~(A+M1+M2)^2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_0mmm
  m3_0_00mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  m3_1_11mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_1_0mmm
  m3_1_01mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_0_mmmm
  m3_1_10mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  #M4 - GIVEN M2 and M3
  
  fit<-glm(as.formula(paste("M4~(A+M2+M3)^2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_0_m0mm
  m4_0_m00m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M2<-m2_1_mmmm
  dat2$M3<-m3_1_m1mm
  m4_1_m11m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  
  #M4 - GIVEN M1 and M3
  
  fit<-glm(as.formula(paste("M4~(A+M1+M3)^2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  dat2$M3<-m3_0_0mmm
  m4_0_0m0m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M3<-m3_1_1mmm
  m4_1_1m1m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  #M4 - GIVEN M1 and M2
  
  fit<-glm(as.formula(paste("M4~(A+M1+M2)^2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_0mmm
  m4_0_00mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  m4_1_11mm<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  ### order 4 conditionals
  
  #M4 - GIVEN M1, M2 and M3
  
  fit<-glm(as.formula(paste("M4~(A+M1+M2+M3)^2+",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_0mmm
  dat2$M3<-m3_0_00mm
  m4_0_000m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  dat2$M3<-m3_1_11mm
  m4_1_111m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_1_0mmm
  dat2$M3<-m3_1_01mm
  m4_1_011m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_1_10mm
  m4_1_101m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  dat2$M3<-m3_0_mmmm
  m4_1_110m<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  ### outcome
  
  # Y
  
  fit<-glm(as.formula(paste("Y~(A+M1+M2+M3+M4)^2 +",confs)),data=data,family=binomial)
  if((!fit$converged)|any(is.na(fit$coefficients))) flag<-TRUE
  
  
  ### ESTIMATE OUTCOME EXPECTATION IN EACH ARM ###
  
  #p_ctr
  
  dat2$A<-0
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_0mmm
  dat2$M3<-m3_0_00mm
  dat2$M4<-m4_0_000m
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_ctr<-mean(y0)
  
  #p_all
  
  dat2$A<-1
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_all<-mean(y1)
  
  #p_trt
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm   
  dat2$M3<-m3_1_11mm   
  dat2$M4<-m4_1_111m   
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_trt<-mean(y1)
  
  #p_1
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_1_mmmm   
  dat2$M3<-m3_1_m1mm   
  dat2$M4<-m4_1_m11m   
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_1<-mean(y0)
  
  #p_2
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_1_1mmm   
  dat2$M4<-m4_1_1m1m   
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_2<-mean(y0)
  
  #p_3
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  dat2$M3<-m3_0_mmmm
  dat2$M4<-m4_1_11mm   
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_3<-mean(y0)
  
  #p_4
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  dat2$M3<-m3_1_11mm
  dat2$M4<-m4_0_mmmm
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_4<-mean(y0)
  
  #p_1_prime
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_1_0mmm #  
  dat2$M3<-m3_1_01mm #
  dat2$M4<-m4_1_011m #
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_1_prime<-mean(y0)
  
  #p_2_prime
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_1_10mm #   
  dat2$M4<-m4_1_101m #    
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_2_prime<-mean(y0)
  
  #p_3_prime
  
  dat2$A<-1
  dat2$M1<-m1_1_mmmm
  dat2$M2<-m2_1_1mmm
  dat2$M3<-m3_0_mmmm
  dat2$M4<-m4_1_110m #  
  
  y0<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_3_prime<-mean(y0)
  
  
  #p_{1,2}
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_1_mmmm   
  dat2$M4<-m4_1_mm1m   
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_12<-mean(y1)
  
  #p_{1,2,3}
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_0_mmmm   
  dat2$M4<-m4_1_mmmm   
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_123<-mean(y1)
  
  #p_{1,2,3,4}
  
  dat2$A<-1
  dat2$M1<-m1_0_mmmm
  dat2$M2<-m2_0_mmmm
  dat2$M3<-m3_0_mmmm   
  dat2$M4<-m4_0_mmmm   
  
  y1<-rbinom(n,1,predict(fit,newdata=dat2,type="response"))
  
  p_1234<-mean(y1)
  
  ### ESTIMATE EFFECTS ###
  
  ### Effects under "one-policy premise"
  # approach (a) 
  
  IIE_1<-p_trt-p_1
  IIE_2<-p_trt-p_2
  IIE_3<-p_trt-p_3
  IIE_4<-p_trt-p_4
  IIE_int<-p_trt-p_all-IIE_1-IIE_2-IIE_3-IIE_4
  
  # approach (b)
  IIE_1_prime<-p_trt-p_1_prime
  IIE_2_prime<-p_trt-p_2_prime
  IIE_3_prime<-p_trt-p_3_prime
  IIE_4_prime<-p_trt-p_4
  IIE_int_prime<-p_trt-p_all-IIE_1_prime-IIE_2_prime-IIE_3_prime-IIE_4_prime
  
  ### Effects under sequential policies
  
  IIE_seq1<-p_trt-p_1
  IIE_seq2<-p_1-p_12
  IIE_seq3<-p_12-p_123
  IIE_seq4<-p_123-p_1234
  IIE_seqfull<-p_trt-p_1234
  IIE_seqint<-p_1234-p_all
  
  ### Interventional direct effect and total causal effect
  
  IDE<-p_all-p_ctr  
  TCE<-p_trt-p_ctr
  
  ### Collect and return results
  
  if(RCT=="one-policy_A")
    res<-c(TCE,IDE,IIE_1,IIE_2,IIE_3,IIE_4, IIE_int) else
      if(RCT=="one-policy_B") 
        res<-c(TCE,IDE,IIE_1_prime,IIE_2_prime,IIE_3_prime,IIE_4_prime, IIE_int_prime) else
          if(RCT=="sequential")
            res<-c(TCE,IDE, IIE_seqfull,IIE_seq1,IIE_seq2,IIE_seq3,IIE_seq4,IIE_seqint) else
              if(RCT=="all")
                res<-c(TCE,IDE, IIE_1,IIE_2,IIE_3,IIE_4, IIE_int, 
                       IIE_1_prime,IIE_2_prime,IIE_3_prime,IIE_4_prime, IIE_int_prime,
                       IIE_seqfull,IIE_seq1,IIE_seq2,IIE_seq3,IIE_seq4,IIE_seqint) else
                         stop("RCT argument must be either 'one-policy_A' or 'one-policy_B' or 'sequential' or 'all'")
  
  if(!flag)return(res) else
    return(rep(NA,length(res)))
}