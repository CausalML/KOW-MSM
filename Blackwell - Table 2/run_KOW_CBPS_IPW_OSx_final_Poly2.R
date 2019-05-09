rm(list=ls())
library(CBPS)
data(Blackwell)
#Blackwell <- Blackwell[Blackwell$time<=4,]









calc_ipw_5 <- function(data){
  
  n <- length(which(data$time==1))
  
  data$Treatment <- data$d.gone.neg
  
  a1 <-as.integer(data$Treatment[data$time==1]==1)
  a0 <-as.integer(data$Treatment[data$time==1]==0)
  
  
  ppr0 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
                camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office
              ,family="binomial",data=data[data$time==1,])$fitted
  ipw0<-rep(0,n)
  ipw0[which(a1==1)] <- 1/ppr0[which(a1==1)]
  ipw0[which(a0==1)] <- 1/(1-ppr0[which(a0==1)])
  ipw1 <-1
  
  ppr1 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
                camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office
              ,family="binomial",data=data[data$time==2,])$fitted
  
  ipw1<-rep(0,n)
  ipw1[which(data$Treatment[data$time==2]==1)] <- 1/ppr1[which(data$Treatment[data$time==2]==1)]
  ipw1[which(data$Treatment[data$time==2]==0)] <- 1/(1-ppr1[which(data$Treatment[data$time==2]==0)])
  
  
  ppr2 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
                camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office
              ,family="binomial",data=data[data$time==3,])$fitted
  
  ipw2<-rep(0,n)
  ipw2[which(data$Treatment[data$time==3]==1)] <- 1/ppr1[which(data$Treatment[data$time==3]==1)]
  ipw2[which(data$Treatment[data$time==3]==0)] <- 1/(1-ppr1[which(data$Treatment[data$time==3]==0)])
  
  
  ppr3 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
                camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office
              ,family="binomial",data=data[data$time==4,])$fitted
  
  ipw3<-rep(0,n)
  ipw3[which(data$Treatment[data$time==4]==1)] <- 1/ppr1[which(data$Treatment[data$time==4]==1)]
  ipw3[which(data$Treatment[data$time==4]==0)] <- 1/(1-ppr1[which(data$Treatment[data$time==4]==0)])
  
  
  ppr4 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
                camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office
              ,family="binomial",data=data[data$time==5,])$fitted
  
  ipw4<-rep(0,n)
  ipw4[which(data$Treatment[data$time==5]==1)] <- 1/ppr1[which(data$Treatment[data$time==5]==1)]
  ipw4[which(data$Treatment[data$time==5]==0)] <- 1/(1-ppr1[which(data$Treatment[data$time==5]==0)])
  
  
  
  sppr0 <- glm(d.gone.neg ~ 1
               ,family="binomial",data=data[data$time==1,])$fitted
  
  sipw0<-rep(0,n)
  sipw0[which(data$Treatment[data$time==1]==1)] <- sppr0[which(data$Treatment[data$time==1]==1)]/ppr0[which(data$Treatment[data$time==1]==1)]
  sipw0[which(data$Treatment[data$time==1]==0)] <- (1-sppr0[which(data$Treatment[data$time==1]==0)])/(1-ppr0[which(data$Treatment[data$time==1]==0)])
  
  
  sppr1 <- glm(d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3
               ,family="binomial",data=data[data$time==2,])$fitted
  
  sipw1<-rep(0,n)
  sipw1[which(data$Treatment[data$time==2]==1)] <- sppr1[which(data$Treatment[data$time==2]==1)]/ppr1[which(data$Treatment[data$time==2]==1)]
  sipw1[which(data$Treatment[data$time==2]==0)] <- (1-sppr1[which(data$Treatment[data$time==2]==0)])/(1-ppr1[which(data$Treatment[data$time==2]==0)])
  
  sppr2 <- glm(d.gone.neg ~ d.gone.neg.l1+ d.gone.neg.l2 + d.neg.frac.l3
               ,family="binomial",data=data[data$time==3,])$fitted
  
  sipw2<-rep(0,n)
  sipw2[which(data$Treatment[data$time==3]==1)] <- sppr2[which(data$Treatment[data$time==3]==1)]/ppr2[which(data$Treatment[data$time==3]==1)]
  sipw2[which(data$Treatment[data$time==3]==0)] <- (1-sppr2[which(data$Treatment[data$time==3]==0)])/(1-ppr2[which(data$Treatment[data$time==3]==0)])
  
  sppr3 <- glm(d.gone.neg ~ d.gone.neg.l1
               ,family="binomial",data=data[data$time==4,])$fitted
  
  sipw3<-rep(0,n)
  sipw3[which(data$Treatment[data$time==4]==1)] <- sppr3[which(data$Treatment[data$time==4]==1)]/ppr3[which(data$Treatment[data$time==4]==1)]
  sipw3[which(data$Treatment[data$time==4]==0)] <- (1-sppr3[which(data$Treatment[data$time==4]==0)])/(1-ppr3[which(data$Treatment[data$time==4]==0)])
  
  sppr4 <- glm(d.gone.neg ~ d.gone.neg.l1+ d.gone.neg.l2 + d.neg.frac.l3
               ,family="binomial",data=data[data$time==5,])$fitted
  
  sipw4<-rep(0,n)
  sipw4[which(data$Treatment[data$time==5]==1)] <- sppr4[which(data$Treatment[data$time==5]==1)]/ppr4[which(data$Treatment[data$time==5]==1)]
  sipw4[which(data$Treatment[data$time==5]==0)] <- (1-sppr4[which(data$Treatment[data$time==5]==0)])/(1-ppr4[which(data$Treatment[data$time==5]==0)])
  
  ipw <- ipw0*ipw1*ipw2*ipw3*ipw4
  sipw <- sipw0*sipw1*sipw2*sipw3*sipw4
  
  return(list(ipw=ipw,sipw=sipw))
  
  
}  

ipw <- calc_ipw_5(Blackwell)$ipw
sipw <- calc_ipw_5(Blackwell)$sipw

####################################################################################################

# Matlab

####################################################################################################
#options(matlab="/usr/local/MATLAB/R2017a/bin/matlab")
options(matlab="/Volumes/Untitled/Applications/MATLAB_R2017a.app/bin/matlab")

PORT <- 9991

Matlab$startServer(port=PORT)

matlab <- Matlab(port=PORT)
isOpen <- open(matlab)
isOpen


#setwd("/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19")
setwd("/Users/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19")

evaluate(matlab, "
         disp('executing gpml startup script...')
         mydir = fileparts(mfilename('fullpath'));                   % where am I located
         addpath(mydir), dirs = {'cov','doc','inf','lik','mean','prior','util'};
         for d = dirs, addpath([mydir,'/',d{:}]), end
         
         "
)

evaluate(matlab, "path")

KGram <- function(XX,AA,t,scale1,scale){
  
  if(t==1){
    spl <- polydot(degree = 2, scale = scale, offset = 1)
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KX
  }
  
  if(t>1){
    #rbf <- rbfdot(sigma = 1)
    lin <- polydot(degree = 1, scale = 1, offset = 1)
    spl <- polydot(degree = 2, scale = scale, offset = 1)
    KA <- kernelMatrix(lin,x = scale(AA))
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KA*KX
  }
  
  return(K)
  
}


gpml <- function(y,XX,AA,t){
  
  #evaluate(matlab, "clear all;")
  
  if(t==1){
    
    XX <- scale(XX)
    
    setVariable(matlab, x=XX)
    setVariable(matlab, y=y)
    
    evaluate(matlab, "
             meanfunc = [];                                               % empty: don't use a mean function
             covfunc = {@covPoly,2}            %  Polynomial kernel
             likfunc = @likGauss;                                         % Gaussian likelihood
             
             hyp = struct('mean', [], 'cov', [1 1], 'lik', -1);
             
             hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
             
             "
    )
    
    hyp2 <- getVariable(matlab, "hyp2")
    
    scale <- (exp(hyp2$hyp2[[2]][1,2]))
    offset <- exp(hyp2$hyp2[[2]][1,1])
    
    gamma <- exp(hyp2$hyp2[[2]][1,2])
    gamma2 <- (exp(gamma))
    
    #sigma2 <- (exp(hyp2$hyp2[[3]][1,1]))^2
    sigma2 <- exp(hyp2$hyp2[[3]])^2
    
    covparam <- hyp2$hyp2[[2]]
    
    evaluate(matlab, "
             clear hyp2;
             clear x;
             clear y;
             ")  
    
    scale1 <- scale
    
  }
  
  
  
  if(t>1){
    
    AA <- cbind(AA)
    AA <- scale(AA)
    XX <- scale(XX)
    XA <- cbind(AA,XX)
    
    setVariable(matlab, x=XX)
    setVariable(matlab, xa=XA)
    setVariable(matlab, y=y)
    setVariable(matlab, a=AA)
    
    maskA <- c(rep(1,dim(AA)[2]),rep(0,dim(XX)[2]))
    maskX <- c(rep(0,dim(AA)[2]),rep(1,dim(XX)[2]))
    
    setVariable(matlab, maskA=maskA)
    setVariable(matlab, maskX=maskX)
    
    evaluate(matlab, "
             meanfunc = [];                                               % empty: don't use a mean function
             covfunc = {@covProd, {{ @covMask,{maskA,{@covPoly,1}}}, {@covMask,{maskX,{@covPoly,2}}}}}             % Product of cov functions: 1. RBF 2. Polynomial
             likfunc = @likGauss;                                         % Gaussian likelihood
             
             hyp = struct('mean', [], 'cov', [1 1 1 1], 'lik', 1);
             
             hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, xa, y);
             
             "
    )
    
    
    
    hyp2 <- getVariable(matlab, "hyp2")
    
    scale <- (exp(hyp2$hyp2[[2]][1,4]))
    scale1 <- (exp(hyp2$hyp2[[2]][1,2]))
    offset <- exp(hyp2$hyp2[[2]][1,3])
    
    # scale <- (exp(hyp2$hyp2[[2]][1,2]))^2
    # offset <- exp(hyp2$hyp2[[2]][1,1])
    
    
    gamma <- exp((hyp2$hyp2[[2]][1,2]))*scale
    gamma2 <- (exp(gamma))^2
    
    #sigma2 <- (exp(hyp2$hyp2[[3]][1,1]))^2
    sigma2 <- exp(hyp2$hyp2[[3]])^2
    
    covparam <- hyp2$hyp2[[2]]
    
    evaluate(matlab, "
             clear hyp2;
             clear x;
             clear y;
             ")  
    
  }
  
  
  
  return(list(gamma=gamma,sigma2=sigma2,lambda=sigma2/gamma2,scale=scale,scale1=scale1,offset=offset))
  
  
}



####################################################################################################

# Optimization Problem

####################################################################################################

#rm(list=ls())

ML <- 1

scale_ll <- LAMB <- NULL
n <- length(Blackwell$d.gone.neg[Blackwell$time==1])

#loop over time to make the Q, A matrix and the c vector for the QP
for(t in 1:max(unique(Blackwell$time))){
  
  if(t==1){
    a1 <-as.integer(Blackwell$d.gone.neg[Blackwell$time==t]==1)
    a0 <-as.integer(Blackwell$d.gone.neg[Blackwell$time==t]==0)
    I1 <- diag(a1)
    I0 <- diag(a0)
    
    
    X <- cbind(  Blackwell$camp.length[which(Blackwell$time==t)],
                 Blackwell$deminc[which(Blackwell$time==t)],
                 Blackwell$base.poll[which(Blackwell$time==t)],
                 Blackwell$year.2002[which(Blackwell$time==t)],
                 Blackwell$year.2004[which(Blackwell$time==t)],
                 Blackwell$year.2006[which(Blackwell$time==t)],
                 Blackwell$base.und[which(Blackwell$time==t)],
                 Blackwell$office[which(Blackwell$time==t)]
    )
    #K <- KGram(X,1,t)
    
    X <- scale(X)
    
    gpmlres <-  gpml(Blackwell$demprcnt[which(Blackwell$time==t)],X,1,t)
    lambda0 <- ML*as.numeric(gpmlres$sigma2)
    K <- KGram(X,1,t,scale1=gpmlres$scale1,scale= gpmlres$scale)
    
    M1 <- eigenMapMatMult(I1,K)
    M0 <- eigenMapMatMult(I0,K)
    #M1 <- I1%*%K
    #M0 <- I0%*%K
    
    tempMa1 <- eigenMapMatMult(M1,I1)
    tempMa0 <- eigenMapMatMult(M0,I0)
    
    #Update Q
    #Q1 <- ((1/n)*(M1%*%I1 + lambda0*I1) + (1/n)*(M0%*%I0 + lambda0*I0))
    Q1 <- ((1/n)*(tempMa1 + lambda0*I1) + (1/n)*(tempMa0 + lambda0*I0))
    
    #Update c
    ones <- rep(1,n)
    tempMa11 <- eigenMapMatMult(M1,ones)
    tempMa10 <- eigenMapMatMult(M0,ones)
    #c1 <- (2/n)*(M1%*%rep(1,n)) + (2/n)*(M0%*%rep(1,n))
    c1 <- (2/n)*(tempMa11) + (2/n)*(tempMa10)
    
    Qt <- ct <- 0    
    Xt <- X
    
    scale_ll[t] <- gpmlres$scale
    LAMB[t] <- lambda0
    
  }#end if t
  
  
  if(t>1){
    
    at1 <-as.integer((Blackwell$d.gone.neg[Blackwell$time==t]==1))
    at0 <-as.integer((Blackwell$d.gone.neg[Blackwell$time==t]==0))
    It1 <- diag(at1)
    It0 <- diag(at0)
    
    Xt <- cbind(Xt,
                Blackwell$camp.length[which(Blackwell$time==t)],
                Blackwell$deminc[which(Blackwell$time==t)],
                Blackwell$base.poll[which(Blackwell$time==t)],
                Blackwell$year.2002[which(Blackwell$time==t)],
                Blackwell$year.2004[which(Blackwell$time==t)],
                Blackwell$year.2006[which(Blackwell$time==t)],
                Blackwell$base.und[which(Blackwell$time==t)],
                Blackwell$office[which(Blackwell$time==t)]
    )
    
    Xt <- scale(Xt)
    
    if(t==2){At <-cbind( Blackwell$d.gone.neg.l1[which(Blackwell$time==t)]
                         
                         
    )}
    if(t>2){At <- cbind(At,Blackwell$d.gone.neg.l1[which(Blackwell$time==t)],
                        Blackwell$d.gone.neg.l2[which(Blackwell$time==t)],
                        Blackwell$d.neg.frac.l3[which(Blackwell$time==t)]
    )}
    
    At <- scale(At)
    
    #K <- KGram(Xt,At,t)
    
    gpmlres <- gpml(Blackwell$demprcnt[which(Blackwell$time==t)],Xt,At,t)
    K <- KGram(Xt,At,t,gpmlres$scale1,gpmlres$scale)
    lambda0 <- ML*as.numeric(gpmlres$sigma2)
    
    ee <- diag(rep(1,n))
    # It1K <- It1%*%K
    # It0K <- It0%*%K
    It1K <- eigenMapMatMult(It1,K)
    It0K <- eigenMapMatMult(It0,K)
    
    tempMa1 <- eigenMapMatMult(It1K,It1)
    tempMa1e <- eigenMapMatMult(It1K,ee)
    tempMa0 <- eigenMapMatMult(It0K,It0)
    tempMa0e <- eigenMapMatMult(It0K,ee)
    tempMaee <- eigenMapMatMult(ee,K)
    tempMaeeK <- eigenMapMatMult(tempMaee,ee)
    
    #Update Q
    # Qt <-  Qt + (   (1/n)*(It1K%*%It1 -2*It1K%*%ee + ee%*%K%*%ee + lambda0*It1) +
    #                   (1/n)*(It0K%*%It0 -2*It0K%*%ee + ee%*%K%*%ee + lambda0*It0)
    # )
    Qt <-  Qt + (   (1/n)*(tempMa1 -2*tempMa1e + tempMaeeK + lambda0*It1) + 
                      (1/n)*(tempMa0 -2*tempMa0e + tempMaeeK + lambda0*It0)
    )
    
    ct <- 2*((lambda0*a1/n) + (lambda0*a0/n))
    
    scale_ll[t] <- gpmlres$scale
    LAMB[t] <- lambda0
    
  }
}



Q <- Q1 + Qt

c <- -(c1+ct)

Ic <- diag(1,n)
IIc<- 1+diag(0,n)
ec<-rep(1,n)
eenc<-(ec/n)%*%t(ec)
Q1c<-(1/(n-1))*(Ic + IIc/n -2*eenc)

model <- list()
model$Q <- Q
model$A <- matrix(1/n, nrow = 1, ncol=n, byrow=T)
#model$A <- matrix(0, nrow = 1, ncol=n, byrow=T)
model$sense <- "="
model$rhs <- c(1)
#model$rhs <- c(0)
model$vtypes <- "C"
model$obj <- c
params <- list(Presolve=2,OutputFlag=0)
#model$lb <- rep(-Inf,n)
model$lb <- rep(0,n)

result <- try(gurobi(model,params))

W <- result$x

close(matlab)


form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 +
camp.length + deminc + base.poll + year.2002 + year.2004 + year.2006 + base.und + office"



fit1<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
            data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE,
            msm.variance = "full", time.vary = FALSE)

fit2<-CBMSM(formula = form1, time=Blackwell$time,id=Blackwell$demName,
            data=Blackwell, type="MSM",  iterations = NULL, twostep = TRUE,
            msm.variance = "approx", time.vary = FALSE)


#################KOW
fit_kow <- lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,weights=W)
fit_kow$coefficients
sqrt(diag(sandwich(fit_kow)))
fit_kow$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_kow)))
fit_kow$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_kow)))


fit_kowc <- lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,weights=W)
fit_kowc$coefficients
sqrt(diag(sandwich(fit_kowc)))
fit_kowc$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_kowc)))
fit_kowc$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_kowc)))
summary(fit_kowc)



#################IPTW
fit_ipw <- lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,weights=ipw)
fit_ipw$coefficients
sqrt(diag(sandwich(fit_ipw)))
fit_ipw$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_ipw)))
fit_ipw$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_ipw)))


fit_ipwc <- lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,weights=ipw)
fit_ipwc$coefficients
sqrt(diag(sandwich(fit_ipwc)))
fit_ipwc$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_ipwc)))
fit_ipwc$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_ipwc)))
summary(fit_ipwc)

#################SIPW
fit_ipws <- lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,weights=sipw)
fit_ipws$coefficients
sqrt(diag(sandwich(fit_ipws)))
fit_ipws$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_ipws)))
fit_ipws$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_ipws)))

fit_ipwsc <- lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,weights=sipw)
fit_ipwsc$coefficients
sqrt(diag(sandwich(fit_ipwsc)))
fit_ipwsc$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_ipwsc)))
fit_ipwsc$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_ipwsc)))
summary(fit_ipwsc)


#################CBPS FULL
fit_cbps <- lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,weights=fit1$weights)
fit_cbps$coefficients
sqrt(diag(sandwich(fit_cbps)))
fit_cbps$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_cbps)))
fit_cbps$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_cbps)))


fit_cbpsc <- lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,weights=fit1$weights)
fit_cbpsc$coefficients
sqrt(diag(sandwich(fit_cbpsc)))
fit_cbpsc$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_cbpsc)))
fit_cbpsc$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_cbpsc)))
summary(fit_cbpsc)



#################CBPS APPROX
fit_cbps2 <- lm(demprcnt[time==1]~fit1$treat.hist,data=Blackwell,weights=fit2$weights)
fit_cbps2$coefficients
sqrt(diag(sandwich(fit_cbps2)))
fit_cbps2$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_cbps2)))
fit_cbps2$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_cbps2)))


fit_cbpsc2 <- lm(demprcnt[time==1]~fit1$treat.cum,data=Blackwell,weights=fit2$weights)
fit_cbpsc2$coefficients
sqrt(diag(sandwich(fit_cbpsc2)))
fit_cbpsc2$coefficients+qnorm(.975)*sqrt(diag(sandwich(fit_cbpsc2)))
fit_cbpsc2$coefficients-qnorm(.975)*sqrt(diag(sandwich(fit_cbpsc2)))
summary(fit_cbpsc2)