

set.seed(1)

#library(StatMatch)
library(kernlab)
library(MASS)
library(CBPS)
library(rbenchmark)



rm(list=ls())


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
MCBPSW <- BCBPS <- ERRCBPS <- MMCBPSW <- MBCBPS <- MSECBPS <- MERCBPS <-  NULL


delta <- .8

n <- 500
#n <- 1000
REP <- 3
nn <- n

ite <- 100
seqq <- c(seq(0,100,length.out=25),seq(101,1000,length.out = 25))

#seqq <- c(seq(0,1000,length.out=10))
# seqqn <- seq(300,1000,length.out = 8)

#seqqn <- c(500,600,700)

scale_ll <- scale_ll2 <- LAMB1 <- LAMB2 <- MSCALE <- MSCALE2 <- MLA1 <- MLA2<- time_temp <- time_tempA <-   NULL


Sigma2 <- matrix(0,nrow=ite,ncol=REP)

seqq <- 1

seqqn <- seq(100,1000,length.out = 10)

seqqn <- 100

jj <- 0


seqqtt <- seq(2,6)

MA_TIME_TEMP  <- MA_TIME_TEMPA  <- matrix(NA,nrow=ite,ncol=length(seqqtt))
TIME_TEMP  <- TIME_TEMPA  <- NULL

for (ttt in seqqtt) {
  
  REP <- ttt
  jj <- jj+1
  #print(jj)
  print(paste("Iteration ", jj, "out of ", length(seqqtt)))
  for(ii in 1:ite){
    print(paste("###############################################################################################################################################################", ii))
    #data <-  read.table(paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_cubic_nointe_031017",nn,"_repp",REP,"_",ii,".csv",sep=""), header=T)	
    #data <- read.table(paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_linear_nointe_081117",nn,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    #data <- read.table(file=paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_211117_nn",nn,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    #data <- read.table(file=paste("~/KOW/data/Longdata3X_L_L_211117_nn",nn,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    data <- read.table(file=paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_211117_TIME",ttt,"_repp",REP,"_",ii,".csv",sep=""), header=T)
    
    
    n <- dim(data)[1]/max(unique(data$time))
    

    
    dat1 <- data.frame(as.factor(data$ID),data$time,data$X1,data$X2,data$X3,data$TreatmentIPW,data$Treatment)
    colnames(dat1) <- c("ID","time","X1","X2","X3","TreatmentIPW","Treatment")
    form0 <- "Treatment ~ TreatmentIPW + X1 + X2 + X3 + TreatmentIPW:X1 + TreatmentIPW:X2 + TreatmentIPW:X3"
    
    time_temp[ii] <- benchmark(fit0 <- CBMSM(formula=form0, time=dat1$time,id=dat1$ID,time.vary=FALSE,data=dat1,type="MSM",iterations = NULL,twostep = TRUE,msm.variance = "full"), replications=1)$elapsed
    time_tempA[ii] <- benchmark(fit0 <- CBMSM(formula=form0, time=dat1$time,id=dat1$ID,time.vary=FALSE,data=dat1,type="MSM",iterations = NULL,twostep = TRUE,msm.variance = "approx"),replications = 1)$elapsed
    
    
    save.image("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6_CBPS.RData")

    
  }  
  
  MA_TIME_TEMP[,jj] <- time_temp  
  MA_TIME_TEMPA[,jj] <- time_tempA  
  
  TIME_TEMP[jj] <- mean(time_temp,na.rm=T)
  TIME_TEMPA[jj] <- mean(time_tempA,na.rm=T)
  
  save.image("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6_CBPS.RData")
  
}     

time_temp

