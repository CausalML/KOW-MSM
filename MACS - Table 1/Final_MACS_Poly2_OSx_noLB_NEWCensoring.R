


library(foreign)
library(survival)
library(ipw)
library(sandwich)
library(kernlab)
library(gurobi)
library(qpcR)
library(R.matlab)
library(VN)
library(ipw)
rm(list=ls())
set.seed(12345)


####################################################################################################

# Matlab

####################################################################################################
#options(matlab="/usr/local/MATLAB/R2017a/bin/matlab")
options(matlab="/Volumes/Untitled/Applications/MATLAB_R2017a.app/bin/matlab")

PORT <- 9987

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

#setwd("~/research/Research/projects/Cornell/Longitudinal/Rcode/case study/MACS/data/stata/analysis/")


#dta <- read.dta("~/research/Research/projects/Cornell/Longitudinal/Rcode/case study/MACS/data/stata/analysis/complete_2001_011217.dta")
dta <- read.dta("/Users/micsan/Desktop/Prova Infcare/complete_2001_180418.dta")


#dta <- dta[dta$firstvisityy>=2001&dta$lastvisityy<=2010&dta$tempv<17,]
dta <- dta[dta$firstvisityy>=2001&dta$lastvisityy<=2005,]
table(dta$death)


scale_ll <- LAMB <- NULL

ite1 <- dta$wbc*dta$treat_ipw
ite2 <- dta$plate*dta$treat_ipw
ite3 <- dta$leu3n*dta$treat_ipw
ite4 <- dta$age*dta$treat_ipw
ite5 <- dta$rbc*dta$treat_ipw


ite12 <- dta$wbc^2*dta$treat_ipw
ite22 <- dta$plate^2*dta$treat_ipw
ite32 <- dta$leu3n^2*dta$treat_ipw
ite42 <- dta$age^2*dta$treat_ipw
ite52 <- dta$rbc^2*dta$treat_ipw






temp <- ipwtm(
  exposure = treat,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ wbc+rbc+plate+leu3n+age+treat_ipw,
  id = caseid,
  timevar = visit,
  type = "all",
  data = dta)
summary(temp$ipw.weights)


temps <- ipwtm(
  exposure = treat,
  family = "binomial",
  link = "logit",
  numerator = ~ treat_ipw,
  denominator = ~ wbc+rbc+plate+leu3n+age+treat_ipw,
  id = caseid,
  timevar = visit,
  type = "all",
  data = dta)
summary(temps$ipw.weights)





tempc <- ipwtm(
  exposure = censor,
  family = "binomial",
  link = "logit",
  numerator = ~ 1,
  denominator = ~ wbc+rbc+plate+leu3n+age+treat_ipw,
  id = caseid,
  timevar = visit,
  type = "all",
  data = dta)
summary(tempc$ipw.weights)



tempcs <- ipwtm(
  exposure = censor,
  family = "binomial",
  link = "logit",
  numerator = ~ treat_ipw,
  denominator = ~ wbc+rbc+plate+leu3n+age+treat_ipw,
  id = caseid,
  timevar = visit,
  type = "all",
  data = dta)
summary(tempc$ipw.weights)

# 
# +I(wbc^2)+I(rbc^2)+I(plate^2)+I(leu3n^2)+I(age^2)
# +ite1+ite2+ite3+ite4+ite5+I(rbc^2*wbc^2)+I(rbc^2*plate^2)+I(rbc^2*leu3n^2)+I(rbc^2*age^2)
# +I(plate^2*wbc^2)+I(plate^2*leu3n^2)+I(plate^2*age^2)+I(leu3n^2*wbc^2)+I(leu3n^2*age^2)+ite12+ite22+ite32+ite42+ite52




ipw     <- temp$ipw.weights
ct_ipw  <- temp$ipw.weights*tempc$ipw.weights


sipw     <- temps$ipw.weights
cts_ipw  <- temps$ipw.weights*tempcs$ipw.weights

trunc_cts_ipw       <- ifelse(ipw< quantile(cts_ipw,.995), cts_ipw, quantile(cts_ipw,.995))





KGram <- function(XX,AA,t,scale1,scale){
  
  if(t==1){
    spl <- polydot(degree = 2, scale = scale, offset = 1)
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KX
  }
  
  if(t>1){
    rbf <- rbfdot(sigma = 1)
    lin <- polydot(degree = 1, scale = 1, offset = 1)
    spl <- polydot(degree = 2, scale = scale, offset = 1)
    #KC <- kernelMatrix(lin,x = scale(CC))
    KA <- kernelMatrix(lin,x = scale(AA))
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KA*1*KX
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




















#ll <- 10
tol <- 1e-08
#tolub <- max(cts_ipw)
ML <- 1
SCALE <- 1

for(t in 1:max(unique(dta$tempv))){
  print(paste("#################################################################################################################",t, " out of ", max(unique(dta$tempv))))
  
  censoring <- 1-dta$death
  
  if(t==1){
    a1 <-as.integer(dta$treat[dta$tempv==t]==1)
    a0 <-as.integer(dta$treat[dta$tempv==t]==0)
    
    c1 <-as.integer(dta$censor[dta$tempv==t]==1)
    c0 <-as.integer(dta$censor[dta$tempv==t]==0)
    
    I1 <- diag(a1)
    I0 <- diag(a0)
    
    C1 <- diag(c1)
    C0 <- diag(c0)
    
    n1 <- dim(diag(as.integer(dta$treat[dta$tempv==t]==1)))[1]
    #wbc+rbc+plate+leu3n+age
    X <- cbind(dta$leu3n[which(dta$tempv==t)],dta$wbc[which(dta$tempv==t)],dta$rbc[which(dta$tempv==t)],
               dta$plate[which(dta$tempv==t)],dta$age[which(dta$tempv==t)])
    X <- scale(X)
    
    #K <- KGram(X,1,t)
    
    gpmlres <-  gpml(dta$time[which(dta$tempv==t)],X,1,t)
    lambda0 <- ML*as.numeric(gpmlres$sigma2)
    K <- KGram(X,1,t,scale1=gpmlres$scale1,scale= gpmlres$scale)
    # lambda0 <- ML
    # K <- KGram(X,1,1,t,scale1=SCALE,scale=SCALE)
    
    
    
    #gpmlres <-  gpml(data$t.censor[which(data$time==t)],X,1,t)
    #print("OK GPML 1")
    #K <- gpmlres$K
    #Sigma2[ii,t] <- gpmlres$sigma2
    IC11 <- eigenMapMatMult(I1,C1) 
    IC10 <- eigenMapMatMult(I1,C0) 
    IC01 <- eigenMapMatMult(I0,C1) 
    IC00 <- eigenMapMatMult(I0,C0) 
    
    M11 <- eigenMapMatMult(IC11,K)
    M10 <- eigenMapMatMult(IC10,K)
    M01 <- eigenMapMatMult(IC01,K)
    M00 <- eigenMapMatMult(IC00,K)
    
    #M1 <- I1%*%K
    #M0 <- I0%*%K
    
    tempMa11 <- eigenMapMatMult(M11,IC11)
    tempMa10 <- eigenMapMatMult(M10,IC10)
    tempMa01 <- eigenMapMatMult(M01,IC01)
    tempMa00 <- eigenMapMatMult(M00,IC00)
    
    #Update Q
    #Q1 <- ( (1/n1)*(tempMa11 + lambda0*IC11) + (1/n1)*(tempMa10 + lambda0*IC10) + (1/n1)*(tempMa01 + lambda0*IC01) + (1/n1)*(tempMa00 + lambda0*IC00))
    Q1 <- ( (1/n1)*(tempMa10 + lambda0*IC10)  + (1/n1)*(tempMa00 + lambda0*IC00))
    
    #Update c
    ones <- rep(1,n1)
    tempMa11 <- eigenMapMatMult(M11,ones)
    tempMa10 <- eigenMapMatMult(M10,ones)
    tempMa01 <- eigenMapMatMult(M01,ones)
    tempMa00 <- eigenMapMatMult(M00,ones)
    #c1 <- (2/n)*(M1%*%rep(1,n)) + (2/n)*(M0%*%rep(1,n))
    #c1 <- -( (2/n1)*(tempMa11) + (2/n1)*(tempMa01) + (2/n1)*(tempMa10)  + (2/n1)*(tempMa00))
    c1 <- -( (2/n1)*(tempMa10)  + (2/n1)*(tempMa00))
    
    Qt <- ct <- 0    
    Xt <- X
    
    scale_ll[t] <- gpmlres$scale
    LAMB[t] <- lambda0
    
    model <- list()
    # model$A          <- matrix(0, nrow=1, ncol=n1, byrow=T)
    # model$rhs        <- c(0)
    model$A          <- matrix(1/n1, nrow=1, ncol=n1, byrow=T)
    model$rhs        <- c(1)
    model$modelsense <- "min"
    model$Q          <- Q1
    model$obj        <- c1
    model$sense      <- c("=")
    model$lb <- rep(tol,n1)
    #model$ub <- rep(tolub,n1)
    model$vtypes <- "C"
    params <- list(Presolve=2,OutputFlag=0)
    
    result <- try(gurobi(model,params))
    if (class(result) != "try-error") {
      WW <-  result$x
      
      Qt <- Q1    
      ct <- 0
      Xt <- X
      length(Xt) <- n1
    }
  }#end if t
  
  if (class(result) != "try-error") {
    if(t>1){
      
      a1 <-as.integer((dta$treat[dta$tempv==t]==1))
      a0 <-as.integer((dta$treat[dta$tempv==t]==0))
      I1 <- diag(a1)
      I0 <- diag(a0)
      
      c1 <-as.integer((dta$censor[dta$tempv==t]==1))
      c0 <-as.integer((dta$censor[dta$tempv==t]==0))
      C1 <- diag(c1)
      C0 <- diag(c0)
      
      nt <- length(a1)
      
      
      Xt <- qpcR:::cbind.na(dta$leu3n[which(dta$tempv==t)],dta$wbc[which(dta$tempv==t)],dta$rbc[which(dta$tempv==t)],
                            dta$plate[which(dta$tempv==t)],dta$age[which(dta$tempv==t)])[1:nt,]
      
      # Xt <- qpcR:::cbind.na(Xt,
      #                       dta$leu3n[which(dta$tempv==t)],dta$wbc[which(dta$tempv==t)],dta$rbc[which(dta$tempv==t)],
      #                       dta$plate[which(dta$tempv==t)],dta$age[which(dta$tempv==t)])[1:nt,]
      
      Xt <- scale(Xt)
      
      if(t==2){
        At <- cbind(dta$treat_ipw[which(dta$tempv==t)])
        At <- scale(At)
        #Ct <- cbind(dta$censor_ipw[which(dta$tempv==t)])
      }
      
      if(t>2){
        At <- qpcR:::cbind.na(dta$treat_ipw[which(dta$tempv==t)])[1:nt,]
        At <- scale(At)
        # At <- qpcR:::cbind.na(At,
        #                       dta$treat_ipw[which(dta$tempv==t)])[1:nt,]
        #Ct <- qpcR:::cbind.na(Ct,dta$censor_ipw[which(dta$tempv==t)])[1:nt,]
      }
      
      #K <- KGram(Xt,At,t)
      
      gpmlres <- gpml(dta$time[which(dta$tempv==t)],Xt,At,t)
      K <- KGram(Xt,At,t,gpmlres$scale1,gpmlres$scale)
      lambda0 <- ML*as.numeric(gpmlres$sigma2)
      # K <- KGram(Xt,At,1,t,SCALE,SCALE)
      # lambda0 <- ML
      
      ee <- diag(rep(1,nt))
      
      IC11 <- eigenMapMatMult(I1,C1) 
      IC10 <- eigenMapMatMult(I1,C0) 
      IC01 <- eigenMapMatMult(I0,C1) 
      IC00 <- eigenMapMatMult(I0,C0) 
      
      M11 <- eigenMapMatMult(IC11,K)
      M10 <- eigenMapMatMult(IC10,K)
      M01 <- eigenMapMatMult(IC01,K)
      M00 <- eigenMapMatMult(IC00,K)
      
      #M1 <- I1%*%K
      #M0 <- I0%*%K
      
      tempMa11 <- eigenMapMatMult(M11,IC11)
      tempMa10 <- eigenMapMatMult(M10,IC10)
      tempMa01 <- eigenMapMatMult(M01,IC01)
      tempMa00 <- eigenMapMatMult(M00,IC00)
      
      
      tempMa11e <- eigenMapMatMult(M11,ee)
      tempMa10e <- eigenMapMatMult(M10,ee)
      tempMa01e <- eigenMapMatMult(M01,ee)
      tempMa00e <- eigenMapMatMult(M00,ee)
      
      tempMaee <- eigenMapMatMult(ee,K)
      tempMaeeK <- eigenMapMatMult(tempMaee,ee)
      
      scale_ll[t] <- gpmlres$scale
      LAMB[t] <- lambda0
      
      Qt <-   (     
                      (1/nt)*(tempMa10 -2*tempMa10e + tempMaeeK + lambda0*IC10) +
                      
                      (1/nt)*(tempMa00 -2*tempMa00e + tempMaeeK + lambda0*IC00)
      )
      
      # Qt <-   Qt[1:nt,1:nt] + (     (1/nt)*(tempMa11 -2*tempMa11e + tempMaeeK + lambda0*IC11) + 
      #                                 (1/nt)*(tempMa10 -2*tempMa10e + tempMaeeK + lambda0*IC10) +
      #                                 (1/nt)*(tempMa01 -2*tempMa01e + tempMaeeK + lambda0*IC01) +
      #                                 (1/nt)*(tempMa00 -2*tempMa00e + tempMaeeK + lambda0*IC00) 
      # )
      
      
      
      ct <- -(2*( (lambda0*a1*c0/nt) + (lambda0*a0*c0/nt) ))
      
      c <- (ct)
      # c <- (c1[1:nt]+ct)
      
      #Update Q
      Q <- Qt
      
      
      model <- list()
      model$Q <- Q
      #model$A <- matrix(1/n, nrow = 1, ncol=n, byrow=T)
      model$A <- matrix(1/nt, nrow = 1, ncol=nt, byrow=T)
      model$sense <- "="
      #model$rhs <- c(2)
      model$rhs <- c(1)
      model$vtypes <- "C"
      model$obj <- c
      params <- list(Presolve=2,OutputFlag=0)
      #model$lb <- rep(-Inf,n)
      model$lb <- rep(tol,nt)
      #model$ub <- rep(tolub,nt)
      params <- list(Presolve=2,OutputFlag=0,QCPDual=0)
      
      
      result <- try(gurobi(model,params))
      
      if (class(result) != "try-error") {
        WW <- append(WW,result$x)
      }
    }
  }#end if result from t=1
}

summary(WW)
close(matlab)






###############################IPTCW
fitiptcw <- coxph(Surv(timestart,timeend,death)~treat,data=dta,weights = ct_ipw,robust=T)
exp(fitiptcw$coeff)
sqrt(fitiptcw$var)
sqrt(diag(sandwich(fitiptcw)))
summary(fitiptcw)$conf.int[3]
summary(fitiptcw)$conf.int[4]
summary(ct_ipw)
var(ct_ipw)


###############################Stable IPTCW
fitsiptcw <- coxph(Surv(timestart,timeend,death)~treat,data=dta,weights = cts_ipw,robust=T)
exp(fitsiptcw$coeff)
sqrt(fitsiptcw$var)
sqrt(diag(sandwich(fitsiptcw)))
summary(fitsiptcw)$conf.int[3]
summary(fitsiptcw)$conf.int[4]
summary(cts_ipw)
var(cts_ipw)


###############################KOW
fitkow <- coxph(Surv(timestart,timeend,death)~treat,data=dta,weights = WW,robust=T)
exp(fitkow$coeff)
sqrt(fitkow$var)
sqrt(diag(sandwich(fitkow)))
summary(fitkow)$conf.int[3]
summary(fitkow)$conf.int[4]
summary(WW)
var(WW)



#save.image("Wcomplete3.RData")

