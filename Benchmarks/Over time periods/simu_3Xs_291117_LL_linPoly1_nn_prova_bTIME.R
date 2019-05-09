

set.seed(1)

#library(StatMatch)
library(kernlab)
library(MASS)
library(gurobi)
library(KRLS)
library(beepr)
library(sandwich)
library(R.matlab)
library(VN)
library(rbenchmark)

rm(list=ls())


source("~/research/Research/projects/Cornell/Longitudinal/Rcode/Linear/calc_ipw2.R")

#######################################################################################################################
#
# Load Matlab and open session
#
# Notes: before these following steps: open matlab ---> from terminal matlab -nodesktop -nosplash
# and run MatlabServer inside Matlab. 
#
#######################################################################################################################
options(matlab="/usr/local/MATLAB/R2017a/bin/matlab")

PORT <- 9989

Matlab$startServer(port=PORT)

matlab <- Matlab(port=PORT)
isOpen <- open(matlab)
isOpen

setwd("/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/cov',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/doc',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/inf',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/lik',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/mean',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/prior',oldpath)")

evaluate(matlab, "oldpath = path;
path('/home/micsan/Documents/MATLAB/gpml-matlab-v4.1-2017-10-19/util',oldpath)")

evaluate(matlab, "
         disp('executing gpml startup script...')
         mydir = fileparts(mfilename('fullpath'));                   % where am I located
         addpath(mydir), dirs = {'cov','doc','inf','lik','mean','prior','util'};
         for d = dirs, addpath([mydir,'/',d{:}]), end
         
         "
)

evaluate(matlab, "path")



##########################################
#Negative log likelihood
##########################################



KGram <- function(XX,AA,t,scale1,scale){
  
  if(t==1){
    spl <- polydot(degree = 1, scale = 1, offset = 1)
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KX
  }
  
  if(t>1){
    #rbf <- rbfdot(sigma = 1)
    lin <- polydot(degree = 1, scale = 1, offset = 1)
    spl <- polydot(degree = 1, scale = 1, offset = 1)
    KA <- kernelMatrix(lin,x = scale(AA))
    KX <- kernelMatrix(spl,x = scale(XX))
    K <- KA*KX
  }
  
  return(K)
  
}


gpml <- function(y,XX,AA,t){
  
  #evaluate(matlab, "clear all;")
  
  if(t==1){
   
    setVariable(matlab, x=XX)
    setVariable(matlab, y=y)

    evaluate(matlab, "
        meanfunc = [];                                               % empty: don't use a mean function
        covfunc = {@covPoly,1}            %  Polynomial kernel
        likfunc = @likGauss;                                         % Gaussian likelihood
        
        hyp = struct('mean', [], 'cov', [2.71 2.71], 'lik', -1);
        
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
        covfunc = {@covProd, {{ @covMask,{maskA,{@covPoly,1}}}, {@covMask,{maskX,{@covPoly,1}}}}}             % Product of cov functions: 1. RBF 2. Polynomial
        likfunc = @likGauss;                                         % Gaussian likelihood

        hyp = struct('mean', [], 'cov', [0 0 2.71 2.71], 'lik', 1);

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




time_gpml_t <- NULL

optimize <- function(data){
  print("CIAO")
  #loop over time to make the Q, A matrix and the c vector for the QP
  for(t in 1:max(unique(data$time))){
    
    
    if(t==1){
      a1 <-as.integer(data$Treatment[data$time==t]==1)
      a0 <-as.integer(data$Treatment[data$time==t]==0)
      I1 <- diag(a1)
      I0 <- diag(a0)
      
      
      X <- cbind(data$X1[which(data$time==t)],data$X2[which(data$time==t)],data$X3[which(data$time==t)])
      time_gpml_t[t] <- benchmark(gpmlres <-  gpml(data$y[which(data$time==t)],X,1,t),replications=1)$elapsed
      
      lambda0 <- as.numeric(gpmlres$sigma2)
      
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
      
      scale_ll[ii] <- gpmlres$scale
      LAMB1[ii] <- lambda0
      
    }#end if t
    
    
    if(t>1){
      
      at1 <-as.integer((data$Treatment[data$time==t]==1))
      at0 <-as.integer((data$Treatment[data$time==t]==0))
      It1 <- diag(at1)
      It0 <- diag(at0)
      
      Xt <- cbind(Xt,
                  data$X1[which(data$time==t)],
                  data$X2[which(data$time==t)],
                  data$X3[which(data$time==t)])
      
      At <- data$TreatmentIPW[which(data$time==t)]
      
      
      time_gpml_t[t] <- benchmark(gpmlres <- gpml(data$y[which(data$time==t)],Xt,At,t),replications=1)$elapsed
      K <- KGram(Xt,At,t,gpmlres$scale1,gpmlres$scale)
      
      lambda0 <- as.numeric(gpmlres$sigma2)
      
      scale_ll2[ii] <- gpmlres$scale
      LAMB2[ii] <- lambda0
      
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
      
    }
  }
  
  
  #Update Q
  Q <- Q1 + Qt
  
  c <- -(c1+ct)
  
  
  model <- list()
  # Ic <- diag(1,n)
  # IIc<- 1+diag(0,n)
  # ec<-rep(1,n)
  # eenc<-(ec/n)%*%t(ec)
  # Q1c<-(1/(n-1))*(Ic + IIc/n -2*eenc)
  Q1c <- diag(1/n,n)
  
  #vvww <- t(ipw_1)%*%diag(1/n,n)%*%(ipw_1)
  
  
  model$A          <- matrix(rep(0,n), nrow=1, byrow=T)
  model$rhs        <- c(0)
  model$modelsense <- "min"
  model$Q          <- Q
  model$obj        <- c
  model$sense      <- c("=")
  model$lb <- rep(0,n)
  model$vtypes <- "C"
  # qc1 <- list()
  # qc1$Qc <- Q1c
  # qc1$rhs <- ll*vvww
  # qc1$sense  <- c("<=")
  
  params <- list(Presolve=2,OutputFlag=0,QCPDual=0)
  
  #model$quadcon <- list(qc1)
  
  
  # model <- list()
  # model$Q <- Q
  # #model$A <- matrix(1/n, nrow = 1, ncol=n, byrow=T)
  # model$A <- matrix(0, nrow = 1, ncol=n, byrow=T)
  # model$sense <- "="
  # #model$rhs <- c(2)
  # model$rhs <- c(0)
  # model$vtypes <- "C"
  # model$obj <- c
  # params <- list(Presolve=2,OutputFlag=0)
  # #model$lb <- rep(-Inf,n)
  # model$lb <- rep(0,n)
  
  time_gurobi <- benchmark(result <- try(gurobi(model,params)),replications = 1)$elapsed
  # 
  # ipw <- calc_ipw(data)$ipw
  # sipw <- calc_ipw(data)$sipw
  
  return(list(result=result,time_gpml=time_gpml_t,time_gurobi=time_gurobi))
  
}





a <- .5
b <- .1
e <- .5
f <- .2



MW <- MW2 <- X1 <- X0 <- B <- WW <- NULL
#X1 <- X0 <- matrix(0,nrow = ite,ncol=repp-1)
tol <- 1e-8
jj <- 0

MERRSIPW <- MERRIPW <- MERR <- MER <- NULL
MSESIPW <- MSEIPW <- MSEMW <- MSEM <- NULL
MBSIPW <- MBIPW <- MB <- MBM <- NULL
MMSIPW <- MMIPW <- MMW2 <- MMM <- NULL


MMW <- BB <- MX1 <- MX0 <- MM <- SDW <- MW3<- ER<- BM <- NULL 
#WX0 <- WX1 <- WA0 <- WA1 <- XXX <- AAA <- matrix(0,nrow=ite,ncol=(repp-1))
d11<-d01<-d00<-d1<-d0<-NULL

MP1<-MP0<-MP00<-MP11<-MP01<-MP10<-ERR<-ERRIPW<-MIPW<-BIPW<-Blm <- BC <- ERRC <- ERRlm <-MC<-SDIPW<-NULL
ERRSIPW <- SDSIPW <- BSIPW <- MSIPW <- MIPWP1 <- MIPWP0 <- MIPWP11 <- MIPWP10 <- MIPWP01 <- MIPWP00 <- NULL
ERRTIPW <- SDTIPW <- BTIPW <- MTIPW <- MTIPWP1 <- MTIPWP0 <- MTIPWP11 <- MTIPWP10 <- MTIPWP01 <- MTIPWP00 <- NULL
MSIPWP1 <- MSIPWP0 <- MSIPWP11 <- MSIPWP10 <- MSIPWP01 <- MSIPWP00 <- NULL
MTIPWP1 <- MTIPWP0 <- MTIPWP11 <- MTIPWP10 <- MTIPWP01 <- MTIPWP00 <- NULL

MMTIPWP1 <- MMTIPWP0 <- MMTIPWP11 <- MMTIPWP00 <- MMTIPWP10 <- MMTIPWP01 <- NULL
MMIPWP1 <- MMIPWP0 <- MMIPWP11 <- MMIPWP00 <- MMIPWP10 <- MMIPWP01 <- NULL
MMSIPWP1 <- MMSIPWP0 <- MMSIPWP11 <- MMSIPWP00 <- MMSIPWP10 <- MMSIPWP01 <- NULL
MMWP1 <- MMWP0 <- MMWP11 <- MMWP00 <- MMWP10 <- MMWP01 <- NULL
MMlmc <- BBC <- SEC <- EERRC <- MMlm <- BBlm <- SElm <- EERRlm <- MMTIPW <- BBTIPW <- SETIPW <- EERRTIPW <- NULL
MMIPW <- BBIPW <- SEIPW <- EERRIPW <- MMSIPW <- BBSIPW <- SESIPW <- EERRSIPW <- MMW <- BB <- SEW <- EERR <- NULL

lambda0 <- theta0 <- sigma0 <- lambda1 <- theta1 <- sigma1 <- NULL

MIPW_2 <- MSIPW_2 <- BIPW_2 <- BSIPW_2 <- ERRIPW_2 <- ERRSIPW_2 <- NULL
MERRIPW_2 <- MERRSIPW_2 <- MSEIPW_2 <- MSESIPW_2 <- MBIPW_2 <- MBSIPW_2 <- MMIPW_2 <- MMSIPW_2 <- NULL

var_opw <- var_ipw <- var_sipw <- NULL
SDM <- SDIPW_2 <- SDSIPW_2 <- COVM <- COVW <- COVIPW <- COVSIPW <- COVIPW_2 <- COVSIPW_2 <- NULL
MCOVM <- MCOVW <- MCOVIPW <- MCOVSIPW <- MCOVIPW_2 <- MCOVSIPW_2 <- NULL


delta <- .8

n <- 500
#n <- 1000
REP <- 3


ite <- 100
seqq <- c(seq(0,100,length.out=25),seq(101,1000,length.out = 25))

#seqq <- c(seq(0,1000,length.out=10))
# seqqn <- seq(900,1000,length.out = 8)

#seqqn <- c(500,600,700)

scale_ll <- scale_ll2 <- LAMB1 <- LAMB2 <- MSCALE <- MSCALE2 <- MLA1 <- MLA2<- time_temp <- NULL
time_gpml <- time_gpml2 <- time_gurobi <- NULL

Sigma2 <- matrix(0,nrow=ite,ncol=REP)

seqq <- 1

#seqqn <- seq(100,1000,length.out = 10)
# seqqn <- c(500,600,700,800,900,1000)
# seqqn <- c(900,1000)

seqqtt <- seq(2,6)


MA_TIME_TEMP  <- MA_TIME_GPML <- MA_TIME_GUROBI  <- matrix(NA,nrow=ite,ncol=length(seqqtt))
TIME_TEMP <- TIME_GPML <- TIME_GUROBI <- NULL

jj <- 0

nn <- 100

for (ttt in seqqtt) {

  REP <- ttt
  jj <- jj+1
  print(nn)
  print(paste("############################################################################################################################################################### Iteration ", jj, "out of ", length(seqqtt)))
  for(ii in 1:ite){
    print(paste("###############################################################################################################################################################", ii))
    #data <-  read.table(paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_cubic_nointe_031017",nn,"_repp",REP,"_",ii,".csv",sep=""), header=T)	
    #data <- read.table(paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_linear_nointe_081117",n,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    data <- read.table(file=paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_211117_TIME",ttt,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    
    n <- dim(data)[1]/max(unique(data$time))
    
    
    #Variables initialization for gurobi
    Q <- matrix(tol,nrow=n,ncol=n)
    c <- 0
    en <- rep(1/n,n)
    
    ee <- diag(rep(1,n))
    
    #lama <- lam0 <- lambda0 <- 0
    IIa <- II0 <- matrix(0,nrow=n,ncol=n)
    A <- NULL
    XL <- XA <- XAIPW <- NULL
    

    
    
    ####################################################################################################
    
    # Optimization Problem
    
    ####################################################################################################
    
    time_temp[ii] <- benchmark(result_temp <- optimize(data),replications=1)$elapsed #end system.time
    
    result <- result_temp$result
    time_gpml[ii] <- mean(result_temp$time_gpml)
    time_gurobi[ii] <- result_temp$time_gurobi
    
    save.image("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6.RData") 
      
    
  }  
  
  TIME_TEMP[jj] <- mean(time_temp,na.rm=T) 
  TIME_GPML[jj] <- mean(time_gpml,na.rm=T) 
  TIME_GUROBI[jj] <- mean(time_gurobi,na.rm=T) 
  
  MA_TIME_TEMP[,jj] <- time_temp  
  MA_TIME_GPML[,jj] <- time_gpml
  MA_TIME_GUROBI[,jj] <- time_gurobi
 
  save.image("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6.RData")
  
}     

close(matlab)

library(beepr)
beep(5)

