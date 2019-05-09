
set.seed(12345)


rm(list=ls())


generate <- function(n,repp,tef0,xef1,xef2,xef3,xef4,xef5,a,e,b1,b2,b3,b4,b5,f1,f2,f3,f4,f5){
  
      
    #I generate an ID for each patient repeated 
    df        <- data.frame(ID = c(1:n), count = rep(repp,n))
    df_long   <- data.frame(ID = rep(df$ID, df$count))
    df_long   <- data.frame(df_long,time = rep(0:(repp-1),n),mounth = ((rep(0:(repp-1),n)-1)*6))
    
    #For subject i, i = 1, . . . , N, we first generated the baseline covariate Li(1) ??? Lognormal(6, 1).
    L1 <- L2 <- L3 <- L4 <- L5 <- rep(NA,repp)
    L1IPW <- L2IPW <- L3IPW <- L4IPW <- L5IPW <- AIPW <- rep(NA,repp)
    A <- rep(NA,repp)
    L1v <- L1vipw <- NULL
    L2v <- L2vipw <- NULL
    L3v <- L3vipw <- NULL
    L4v <- L4vipw <- NULL
    L5v <- L5vipw <- NULL
    Av <- Avipw <- NULL
    Dv <- Rv <-  NULL
    

    ####################################################################################################
    
    # Time-varying variable and treatment
    
    ####################################################################################################
    for (i in 1:n){
      
      L1[1] <- NA
      L1[2] <- rnorm(1,0,1)
      L2[1] <- NA
      L2[2] <- rnorm(1,0,1)
      L3[1] <- NA
      L3[2] <- rnorm(1,0,1)
      L4[1] <- NA
      L4[2] <- rnorm(1,0,1)
      L5[1] <- NA
      L5[2] <- rnorm(1,0,1)
      
      A[1] <- 0
      AIPW[1] <- 0
      
      for (j in 2:(repp)){
        
        #pr    <- (1+exp(-(  a + e*A[j-1] + .3*b1*(L1[j]^3) + .3*b2*(L2[j]^3) + .3*b3*(L3[j]^3)) ))^(-1)
        pr    <- (1+exp(-(  a + e*A[j-1] + b1*(L1[j]) + b2*((L2[j])) + b3*(L3[j]) + b4*(L4[j]) + b5*(L5[j])
                            + .2*A[j-1]*L1[j] + .2*A[j-1]*L2[j] + .2*A[j-1]*L3[j]
                            )))^(-1)
        A[j]  <- rbinom(1,1,pr)
        AIPW[j] <- A[j-1]
        if(j<(repp)){ 
          #AIPW[j-1] <- 0
          L1[j+1]  <- L1[j] + f1*A[j]  + rnorm(1,0,1)
          L1IPW[j+1] <- L1[j]
          L2[j+1]  <- L2[j] + f2*A[j]  + rnorm(1,0,1)
          L2IPW[j+1] <- L2[j]
          L3[j+1]  <- L3[j] + f3*A[j]  + rnorm(1,0,1)
          L3IPW[j+1] <- L3[j]
          L4[j+1]  <- L4[j] + f4*A[j]  + rnorm(1,0,1)
          L4IPW[j+1] <- L4[j]
          L5[j+1]  <- L5[j] + f5*A[j]  + rnorm(1,0,1)
          L5IPW[j+1] <- L5[j]
        }#end if
        
      }#end for j
      
      Dv <- c(Dv,D)
      L1v <- c(L1v,L1)
      L2v <- c(L2v,L2)
      L3v <- c(L3v,L3)
      L4v <- c(L4v,L4)
      L5v <- c(L5v,L5)
      L1vipw <- c(L1vipw,L1IPW)
      L2vipw <- c(L2vipw,L2IPW)
      L3vipw <- c(L3vipw,L3IPW)
      L4vipw <- c(L4vipw,L4IPW)
      L5vipw <- c(L5vipw,L5IPW)
      Avipw <- c(Avipw,AIPW)
      Av <- c(Av,A)
      
    }#end for
    
    
    #I remove the na``
    df_long   <- na.omit(data.frame(df_long, X1 = L1v, X2 = L2v, X3 = L3v, X4 = L4v, X5 = L5v, TreatmentIPW = Avipw, Treatment = Av))
    
    #I sum the treatment and the covariate over time
    aggT <- aggregate(df_long$Treatment,FUN=sum,by=list(df_long$ID))
    aggL1 <- aggregate(df_long$X1,FUN=sum,by=list(df_long$ID))
    aggL2 <- aggregate(df_long$X2,FUN=sum,by=list(df_long$ID))
    aggL3 <- aggregate(df_long$X3,FUN=sum,by=list(df_long$ID))
    aggL4 <- aggregate(df_long$X4,FUN=sum,by=list(df_long$ID))
    aggL5 <- aggregate(df_long$X5,FUN=sum,by=list(df_long$ID))
    
    df_long   <- data.frame(df_long,sumT = rep(aggT$x, df$count-1))
    df_long   <- data.frame(df_long,sumL1 = rep(aggL1$x, df$count-1))
    df_long   <- data.frame(df_long,sumL2 = rep(aggL2$x, df$count-1))
    df_long   <- data.frame(df_long,sumL3 = rep(aggL3$x, df$count-1))
    df_long   <- data.frame(df_long,sumL4 = rep(aggL4$x, df$count-1))
    df_long   <- data.frame(df_long,sumL5 = rep(aggL5$x, df$count-1))
    
    #I order the dataset by time
    dflong_time <- df_long[order(df_long$time),]

    #I sumulate the outcome
    aas <- dflong_time$sumT[dflong_time$time==1]
    bbs1 <- dflong_time$sumL1[dflong_time$time==1]
    bbs2 <- dflong_time$sumL2[dflong_time$time==1]
    bbs3 <- dflong_time$sumL3[dflong_time$time==1]
    bbs4 <- dflong_time$sumL4[dflong_time$time==1]
    bbs5 <- dflong_time$sumL5[dflong_time$time==1]
    

    aalp <- -tef0*mean(aas) - xef1*mean(bbs1) - xef2*mean(bbs2) - xef3*mean(bbs3) - .05*mean(bbs1*bbs2) - .05*mean(bbs1*bbs3) - .05*mean(bbs2*bbs3) - .05*mean(bbs1*bbs4)  - .05*mean(bbs2*bbs4) - .05*mean(bbs3*bbs4)
    - .05*mean(bbs1*bbs5)  - .05*mean(bbs2*bbs5) - .05*mean(bbs3*bbs5) - .05*mean(bbs4*bbs5)
    
    y <- aalp + tef0*aas + xef1*bbs1 + xef2*bbs2 + xef3*bbs3 + .05*bbs1*bbs2 + .05*bbs1*bbs3 + .05*bbs2*bbs3 + .05*bbs1*bbs4 + .05*bbs2*bbs4 + .05*bbs3*bbs4 + .05*bbs1*bbs5 + .05*bbs2*bbs5 + .05*bbs3*bbs5 + .05*bbs4*bbs5   + rnorm(n,0,5)

    #data <- data.frame(df_long,y,aas,bbs)
    #data <- data.frame(dflong_time,y,aas,bbs1,bbs2,bbs3,wrong_bbs1,wrong_bbs2,wrong_bbs3)
    data <- data.frame(dflong_time,y,aas,bbs1,bbs2,bbs3,bbs4,bbs5)
    data
}


n <- 500


a <- .5
b <- .1
b1 <- .05
b2 <- .08
b3 <- -.03
b4 <- -.03
b5 <- -.03
e <- .5
f <- .2
f1 <- .1
f2 <- .1
f3 <- .1
f4 <- .1
f5 <- .1
tef0 <- .8
xef1 <- .5
xef2 <- .5
xef3 <- .5
xef4 <- .5
xef5 <- .5

# xef1 <- .4
# xef2 <- .4
# xef3 <- .4


ite <- 100

seqqtt <- seq(7,10)
#seqqn <- c(100,1000)

nn <- 100

tt <- 9
REP <- tt
repp <- REP+1

  for(i in 1:ite){
    
    anal.data <- generate(nn,tt,tef0,xef1,xef2,xef3,xef4,xef5,a,e,b1,b2,b3,b4,b5,f1,f2,f3,f4,f5)
    #write.table (anal.data, file=paste("~/KOW/data/Longdata3X_L_L_211117_nn",nn,"_repp",REP,"_",i,".csv",sep=""), row.names=F)
    write.table (anal.data, file=paste("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_211117_VAR5",tt,"_repp",REP,"_",i,".csv",sep=""), row.names=F)
    #print(i)
    
  }
