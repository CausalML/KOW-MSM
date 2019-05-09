
DF <- 3
LWD <- 3
CEX <- 1.5


par(mfrow=c(2,2))
setwd("/Users/micsan/Desktop/check/KOW - Natale/data/")

XUB <- 500

######################################################################

# BIAS LL 1

######################################################################

    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly1_LL.RData")
    load("simu_141217_LL_linpoly1_LL.RData")
    data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
    write.table(data1,"one")
    
    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly1_nn_CBPS.RData")
    load("simu_141217_LL_linpoly1_nn_CBPS.RData")
    data2 <- cbind(MBCBPS[5],MERCBPS[5])
    
    write.table(data2,"two")
    
    one <- read.table("one")
    two <- read.table("two")
    
    #BIAS
    ratioB  <- one$MBIPW^2/one$MB^2
    ratioBs <- one$MBSIPW^2/one$MB^2
    ratioBc <- two$V1^2/one$MB^2
    
    #seqqL[1] <- 1
    
  
    par(mar=c(4,4,0,4))
    layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,2,1,2))
    plot.new()
    text(0.5,0.5,"Linear - Correct",cex=2,font=2)
    plot((seqqL[1:length(ratioB)]),ratioB,type="l",lwd=LWD,ylim = c(0,10),xlim=c(0,XUB),xlab="",ylab="")
    mtext(text=expression( ~ lambda[t]),side=1,line=3)
    mtext(text="Ratio Bias",side=2,line=3)
    #mtext("Correct", side=3, line=3, cex=1)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioBs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=2)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioBc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=3)
    abline(h=1)
    
    
    plot.new()
    text(0.5,0.5,"Linear - Overspecified",cex=2,font=2)


######################################################################
    
# BIAS LL 2
    
######################################################################    
    
    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly2_LL.RData")
    load("simu_141217_LL_linpoly2_LL.RData")
    data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
    write.table(data1,"one")
    
    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly2_nn_CBPS.RData")
    load("simu_141217_LL_linpoly2_nn_CBPS.RData")
    
    data2 <- cbind(MBCBPS[5],MERCBPS[5])
    write.table(data2,"two")
    
    one <- read.table("one")
    two <- read.table("two")
    
    #BIAS
    ratioB  <- one$MBIPW_2^2/one$MB^2
    ratioBs <- one$MBSIPW_2^2/one$MB^2
    ratioBc <- two$V1^2/one$MB^2
    
    
    
    
    plot((seqqL[1:length(ratioB)]),ratioB,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,10), xlim=c(0,XUB))
    mtext(text=expression( ~ lambda[t]),side=1,line=3)
    mtext(text="Ratio Bias",side=2,line=3)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioBs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=2)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioBc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=3)
    abline(h=1)


######################################################################

# MSE LL 1

######################################################################




    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly1_LL.RData")
    load("simu_141217_LL_linpoly1_LL.RData")
    data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
    write.table(data1,"one")
    
    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly1_nn_CBPS.RData")
    load("simu_141217_LL_linpoly1_nn_CBPS.RData")
    data2 <- cbind(MBCBPS[5],MERCBPS[5])
    
    write.table(data2,"two")
    
    one <- read.table("one")
    two <- read.table("two")
    
    
    #MSE
    ratioMSE  <- one$MERRIPW/one$MERR
    ratioMSEs <- one$MERRSIPW/one$MERR
    ratioMSEc <- two$V2/one$MERR
    
    plot((seqqL[1:length(ratioB)]),ratioMSE,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,6.5),xlim=c(0,XUB))
    mtext(text=expression( ~ lambda[t]),side=1,line=3)
    mtext(text="Ratio MSE",side=2,line=3)
    
    abline(h=1)
    abline(v=29,lwd=LWD) #mean values of lambda across simulations (computed separatly)
    lines((seqqL[1:length(ratioB)]),ratioMSEs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=2)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioMSEc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=3)
    abline(h=1)


######################################################################

# MSE LL 2

######################################################################


    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly2_LL.RData")
    load("simu_141217_LL_linpoly2_LL.RData")
    
    
    data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
    write.table(data1,"one")
    
    #load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_LL_linpoly2_nn_CBPS.RData")
    load("simu_141217_LL_linpoly2_nn_CBPS.RData")
    
    data2 <- cbind(MBCBPS[5],MERCBPS[5])
    
    write.table(data2,"two")
    
    one <- read.table("one")
    two <- read.table("two")
    
    #MSE
    ratioMSE  <- one$MERRIPW_2/one$MERR
    ratioMSEs <- one$MERRSIPW_2/one$MERR
    ratioMSEc <- two$V2/one$MERR
    
    plot((seqqL[1:length(ratioB)]),ratioMSE,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,4),xlim=c(0,XUB))
    mtext(text=expression( ~ lambda[t]),side=1,line=3)
    mtext(text="Ratio MSE",side=2,line=3)
    abline(v=29,lwd=LWD) #mean values of lambda across simulations (computed separatly)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioMSEs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=2)
    abline(h=1)
    lines((seqqL[1:length(ratioB)]),ratioMSEc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=3)
    abline(h=1)

















par(mfrow=c(2,2))

######################################################################

# BIAS NLNL 1

######################################################################

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly1_NLNL.RData")
load("simu_141217_NLNL_linpoly1_NLNL.RData")


data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly1_nn_CBPS.RData")
load("simu_141217_NLNL_linpoly1_nn_CBPS.RData")

data2 <- cbind(MBCBPS[5],MERCBPS[5])

write.table(data2,"two")

one <- read.table("one")
two <- read.table("two")

#BIAS
ratioB  <- one$MBIPW^2/one$MB^2
ratioBs <- one$MBSIPW^2/one$MB^2
ratioBc <- two$V1^2/one$MB^2


#seqqL[1] <- 1

par(mar=c(4,4,0,4))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,2,1,2))
plot.new()
text(0.5,0.5,"Nonlinear - Misspecified",cex=2,font=2)
plot((seqqL[1:length(ratioB)]),ratioB,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,5))
mtext(text=expression( ~ lambda[t]),side=1,line=3)
mtext(text="Ratio Bias",side=2,line=3)
abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioBs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=2)
abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioBc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=3)
abline(h=1)



plot.new()
text(0.5,0.5,"Nonlinear - Correct",cex=2,font=2)




######################################################################

# BIAS NLNL 2

######################################################################


#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly2_NLNL.RData")
load("simu_141217_NLNL_linpoly2_NLNL.RData")


data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data1,"one")

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly2_nn_CBPS.RData")
load("simu_141217_NLNL_linpoly2_nn_CBPS.RData")

data2 <- cbind(MBCBPS[5],MERCBPS[5])

write.table(data2,"two")


load("nn_r2000_a.RData")
data3 <- cbind(MBIPW_2,MBSIPW_2,MERRIPW_2,MERRSIPW_2)
write.table(data3,"three")


one <- read.table("one")
two <- read.table("two")
three <- read.table("three")

#BIAS
ratioB  <- three$MBIPW_2[5]^2/one$MB^2
ratioBs <- three$MBSIPW_2[5]^2/one$MB^2
ratioBc <- two$V1^2/one$MB^2


seqqL[1] <- 1

plot((seqqL[1:length(ratioB)]),ratioB,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,30))
mtext(text=expression( ~ lambda[t]),side=1,line=3)
mtext(text="Ratio Bias",side=2,line=3)

abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioBs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=2)
abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioBc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",type="l",lwd=LWD,lty=3)
# abline(h=1)




######################################################################

# MSE NLNL 1

######################################################################


#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly1_NLNL.RData")
load("simu_141217_NLNL_linpoly1_NLNL.RData")


data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly1_nn_CBPS.RData")
load("simu_141217_NLNL_linpoly1_nn_CBPS.RData")

data2 <- cbind(MBCBPS[5],MERCBPS[5])

write.table(data2,"two")

one <- read.table("one")
two <- read.table("two")


#MSE
ratioMSE  <- one$MERRIPW/one$MERR
ratioMSEs <- one$MERRSIPW/one$MERR
ratioMSEc <- two$V2/one$MERR

plot((seqqL[1:length(ratioB)]),ratioMSE,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,3))
mtext(text=expression( ~ lambda[t]),side=1,line=3)
mtext(text="Ratio MSE",side=2,line=3)

abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioMSEs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=2)
abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioMSEc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=3)
# abline(h=1)
abline(v=548,lwd=LWD) #mean values of lambda across simulations (computed separatly)






######################################################################

# MSE NLNL 2

######################################################################

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly2_NLNL.RData")
load("simu_141217_NLNL_linpoly2_NLNL.RData")

data1 <- cbind(MB,MERR)
write.table(data1,"one")

#load("/run/user/1000/gvfs/sftp:host=192.168.120.50/home/likewise-open/IMM/micsan/KOW/data/simu_141217_NLNL_linpoly2_nn_CBPS.RData")
load("simu_141217_NLNL_linpoly2_nn_CBPS.RData")

data2 <- cbind(MBCBPS[5],MERCBPS[5])
write.table(data2,"two")


load("nn_r2000_a.RData")
data3 <- cbind(MBIPW_2,MBSIPW_2,MERRIPW_2,MERRSIPW_2)
write.table(data3,"three")


one <- read.table("one")
two <- read.table("two")
three <- read.table("three")


#MSE
ratioMSE  <- three$MERRIPW_2[5]/one$MERR
ratioMSEs <- three$MERRSIPW_2[5]/one$MERR
ratioMSEc <- two$V2/one$MERR

plot((seqqL[1:length(ratioB)]),ratioMSE,xlab="",ylab="",type="l",lwd=LWD,ylim = c(0,6))
mtext(text=expression( ~ lambda[t]),side=1,line=3)
mtext(text="Ratio MSE",side=2,line=3)

abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioMSEs,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=2)
abline(h=1)
lines((seqqL[1:length(ratioB)]),ratioMSEc,xlab=expression(paste(lambda)),main="KOM vs Stable IPW",ylab="MSE SIPW / MSE KOM",type="l",lwd=LWD,lty=3)
# abline(h=1)
abline(v=82,lwd=LWD) #mean values of lambda across simulations (computed separatly)

