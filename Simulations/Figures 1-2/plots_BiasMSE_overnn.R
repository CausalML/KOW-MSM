
library(splines)


setwd("/Users/micsan/Desktop/KOW - Natale/data/")
#setwd("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/")


DF <- 3
LWD <- 3
CEX <- 1.5


par(mfrow=c(2,2))


######################################################################

# BIAS LL 1

######################################################################


##################################################################################################################################################################
#LInear Correct
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_LL_linpoly1_nn.RData")
load("simu_291117_LL_linpoly1_nn_2.RData")    

data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

load("simu_291117_LL_linpoly1_nn_3b.RData")  

data2 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data2,"two")

load("simu_291117_LL_linpoly1_nn_3.RData")  

data3 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data3,"three")

load("simu_141217_LL_linpoly1_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one,two,three)
data2 <- rbind(four)

seqqn <- seq(100,1000,length.out = 10)

MB <- data$MB
MBIPW <- data$MBIPW
MBSIPW <- data$MBSIPW
MBCBPS <- data2$MBCBPS
MERCBPS <- data2$MERCBPS

uylim <- max(max(MB,MBIPW,MBSIPW,MBCBPS))
lylim <- min(min(MB,MBIPW,MBSIPW,MBCBPS))


par(mar=c(4,4,0,4))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,2,1,2))
plot.new()
text(0.5,0.5,"Linear - Correct",cex=2,font=2)
plot( smooth.spline(seqqn[1:length(MB)],MB,df=DF),ylim=c(lylim,uylim),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="Bias",side=2,line=3)
lines(smooth.spline(seqqn[1:length(MB)],MBIPW,df=DF),col=1,lwd=LWD,lty=2)
lines(smooth.spline(seqqn[1:length(MB)],MBSIPW,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MBCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)


plot.new()
text(0.5,0.5,"Linear - Overspecified",cex=2,font=2)


######################################################################

# BIAS LL 2

######################################################################



##################################################################################################################################################################
#LInear Overspecified
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_LL_linpoly2_nn.RData")
load("simu_291117_LL_linpoly2_nn_2.RData")   

data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data1,"one")

load("simu_291117_LL_linpoly2_nn_2b.RData")  

data2 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data2,"two")

load("simu_141217_LL_linpoly2_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
four <- read.table("four")


data <- rbind(one,two)
data2 <- rbind(four)

MB <- data$MB
MBIPW_2 <- data$MBIPW_2
MBSIPW_2 <- data$MBSIPW_2
MBCBPS <- data2$MBCBPS
MERCBPS <- data2$MERCBPS

seqqn <- seq(100,1000,length.out = 10)
MB
MBIPW_2
MBSIPW_2

uylim <- max(max(MB,MBIPW_2,MBSIPW_2,MBCBPS))
lylim <- min(min(MB,MBIPW_2,MBSIPW_2,MBCBPS))

plot( smooth.spline(seqqn[1:length(MB)],MB,df=DF),ylim=c(lylim,uylim),type="l",lwd=LWD,  xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="Bias",side=2,line=3)
lines(smooth.spline(seqqn[1:length(MB)],MBIPW_2,df=DF),col=1,lwd=LWD,lty=2)
lines(smooth.spline(seqqn[1:length(MB)],MBSIPW_2,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MBCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)


######################################################################

# MSE LL 1

######################################################################



load("simu_291117_LL_linpoly1_nn_2.RData")    

data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

load("simu_291117_LL_linpoly1_nn_3b.RData")  

data2 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data2,"two")

load("simu_291117_LL_linpoly1_nn_3.RData")  

data3 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data3,"three")

load("simu_141217_LL_linpoly1_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one,two,three)
data2 <- rbind(four)

seqqn <- seq(100,1000,length.out = 10)



MERR <- data$MERR
MERRIPW <- data$MERRIPW
MERRSIPW <- data$MERRSIPW
MERCBPS <- data2$MERCBPS

uylim <- max(max(MERR,MERRIPW,MERRSIPW,MERCBPS))
lylim <- min(min(MERR,MERRIPW,MERRSIPW,MERCBPS))

plot(smooth.spline(seqqn[1:length(MERR)],MERR,df=DF),ylim=c(0,uylim),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="MSE",side=2,line=3)
lines( smooth.spline(seqqn[1:length(MERR)],MERRIPW,df=DF),col=1,lwd=LWD,lty=2)
lines( smooth.spline(seqqn[1:length(MERR)],MERRSIPW,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MERCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)



######################################################################

# MSE LL 2

######################################################################



##################################################################################################################################################################
#LInear Overspecified
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_LL_linpoly2_nn.RData")
load("simu_291117_LL_linpoly2_nn_2.RData")   

data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data1,"one")

load("simu_291117_LL_linpoly2_nn_2b.RData")  

data2 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data2,"two")

load("simu_141217_LL_linpoly2_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
four <- read.table("four")


data <- rbind(one,two)
data2 <- rbind(four)

MB <- data$MB
MBIPW_2 <- data$MBIPW_2
MBSIPW_2 <- data$MBSIPW_2
MBCBPS <- data2$MBCBPS
MERCBPS <- data2$MERCBPS

seqqn <- seq(100,1000,length.out = 10)



MERR <- data$MERR
MERRIPW_2 <- data$MERRIPW_2
MERRSIPW_2 <- data$MERRSIPW_2

uylim <- max(max(MERR,MERRIPW_2,MERRSIPW_2))
lylim <- min(min(MERR,MERRIPW_2,MERRSIPW_2))

plot(smooth.spline(seqqn[1:length(MERR)],MERR,df=DF),ylim=c(0,uylim),type="l",lwd=LWD,  xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="MSE",side=2,line=3)
lines( smooth.spline(seqqn[1:length(MERR)],MERRIPW_2,df=DF),col=1,lwd=LWD,lty=2)
lines( smooth.spline(seqqn[1:length(MERR)],MERRSIPW_2,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MERCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)

# 
# MCOVM
# MCOVIPW
# MCOVSIPW
# 
# uylim <- max(max(MCOVM,MCOVIPW,MCOVSIPW))
# lylim <- min(min(MCOVM,MCOVIPW,MCOVSIPW))
# 
# plot(seqqn[1:length(MB)],MCOVM,ylim=c(lylim,1),type="l")
# lines(seqqn[1:length(MB)],MCOVIPW,col=2)
# lines(seqqn[1:length(MB)],MCOVSIPW,col=4)
# abline(h=0.95)

















##NON LInear Misspecified
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_NLNL_linpoly1_nn.RData")
#load("simu_291117_NLNL_linpoly1_nn_2.RData")
load("simu_291117_NLNL_linpoly1_nn_3b.RData")  

data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

load("simu_291117_NLNL_linpoly1_nn_4b.RData")  

data2 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data2,"two")

load("simu_291117_NLNL_linpoly1_nn_3.RData")
data3 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data3,"three")

load("simu_141217_NLNL_linpoly1_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one,two,three)
data2 <- rbind(four)

seqqn <- seq(100,1000,length.out = 10)

par(mfrow=c(2,2))


MB <- data$MB
MBIPW <- data$MBIPW
MBSIPW <- data$MBSIPW
MBCBPS <- data2$MBCBPS
MERCBPS <- data2$MERCBPS

uylim <- max(max(MB,MBIPW,MBSIPW,MBCBPS))
lylim <- min(min(MB,MBIPW,MBSIPW,MBCBPS))


par(mar=c(4,4,0,4))
layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,2,1,2))
plot.new()
text(0.5,0.5,"Nonlinear - Misspecified",cex=2,font=2)
plot( smooth.spline(seqqn[1:length(MB)],MB,df=DF),ylim=c(lylim,4),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="Bias",side=2,line=3)
lines(smooth.spline(seqqn[1:length(MB)],MBIPW,df=DF),col=1,lwd=LWD,lty=2)
lines(smooth.spline(seqqn[1:length(MB)],MBSIPW,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MBCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)




plot.new()
text(0.5,0.5,"Nonlinear - Correct",cex=2,font=2)



#Non LInear Correct
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_NLNL_linpoly2_nn_2.RData")
#load("simu_291117_NLNL_linpoly2_nn_2.RData")
load("simu_291117_NLNL_linpoly2_nn_3.RData")
data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data1,"one")

load("nn_r2000_a.RData")
data3 <- cbind(MBIPW_2,MBSIPW_2,MERRIPW_2,MERRSIPW_2)
write.table(data3,"three")


load("simu_141217_NLNL_linpoly2_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)
write.table(data4,"four")

one <- read.table("one")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one)
data2 <- rbind(four)

MB <- data$MB
MBIPW_2 <- three$MBIPW_2
MBSIPW_2 <- three$MBSIPW_2
MBCBPS <- data2$MBCBPS
MERCBPS <- data2$MERCBPS

uylim <- max(max(MB,MBIPW_2,MBSIPW_2,MBCBPS))
lylim <- min(min(MB,MBIPW_2,MBSIPW_2,MBCBPS))

plot( smooth.spline(seqqn[1:length(MB)],MB,df=DF),ylim=c(lylim,uylim),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="Bias",side=2,line=3)
lines(smooth.spline(seqqn[1:length(MBIPW_2)],MBIPW_2,df=DF),col=1,lwd=LWD,lty=2)
lines(smooth.spline(seqqn[1:length(MBIPW_2)],MBSIPW_2,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MBCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)




##NON LInear Misspecified
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_NLNL_linpoly1_nn.RData")
#load("simu_291117_NLNL_linpoly1_nn_2.RData")
load("simu_291117_NLNL_linpoly1_nn_3b.RData")  

data1 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data1,"one")

load("simu_291117_NLNL_linpoly1_nn_4b.RData")  

data2 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data2,"two")

load("simu_291117_NLNL_linpoly1_nn_3.RData")
data3 <- cbind(MB,MBIPW,MBSIPW,MERR,MERRIPW,MERRSIPW)
write.table(data3,"three")

load("simu_141217_NLNL_linpoly1_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)

write.table(data4,"four")

one <- read.table("one")
two <- read.table("two")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one,two,three)
data2 <- rbind(four)

seqqn <- seq(100,1000,length.out = 10)





MERR <- data$MERR
MERRIPW <- data$MERRIPW
MERRSIPW <- data$MERRSIPW

uylim <- max(max(0.5,MERR,MERRIPW,MERRSIPW))
lylim <- min(min(MERR,MERRIPW,MERRSIPW))

plot(smooth.spline(seqqn[1:length(MERR)],MERR,df=DF),ylim=c(0,uylim),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="MSE",side=2,line=3)
lines( smooth.spline(seqqn[1:length(MERR)],MERRIPW,df=DF),col=1,lwd=LWD,lty=2)
lines( smooth.spline(seqqn[1:length(MERR)],MERRSIPW,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MERCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)




#Non LInear Correct
#load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/simu_291117_NLNL_linpoly2_nn_2.RData")
#load("simu_291117_NLNL_linpoly2_nn_2.RData")
load("simu_291117_NLNL_linpoly2_nn_3.RData")
data1 <- cbind(MB,MBIPW_2,MBSIPW_2,MERR,MERRIPW_2,MERRSIPW_2)
write.table(data1,"one")

load("nn_r2000_a.RData")
data3 <- cbind(MBIPW_2,MBSIPW_2,MERRIPW_2,MERRSIPW_2)
write.table(data3,"three")


load("simu_141217_NLNL_linpoly2_nn_CBPS.RData")
data4 <- cbind(MBCBPS,MERCBPS)
write.table(data4,"four")

one <- read.table("one")
three <- read.table("three")
four <- read.table("four")

data <- rbind(one)
data2 <- rbind(four)

MERR <- data$MERR
MERRIPW_2 <- three$MERRIPW_2
MERRSIPW_2 <- three$MERRSIPW_2


uylim <- max(max(MERR,MERRIPW_2,MERRSIPW_2))
lylim <- min(min(MERR,MERRIPW_2,MERRSIPW_2))

plot(smooth.spline(seqqn[1:length(MERR)],MERR,df=DF),ylim=c(0,uylim),type="l",lwd=LWD, xlab = "", ylab="",cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
mtext(text="Sample size",side=1,line=3)
mtext(text="MSE",side=2,line=3)
lines( smooth.spline(seqqn[1:length(MERRIPW_2)],MERRIPW_2,df=DF),col=1,lwd=LWD,lty=2)
lines( smooth.spline(seqqn[1:length(MERRIPW_2)],MERRSIPW_2,df=DF),col=1,lwd=LWD,lty=3)
lines(smooth.spline(seqqn[1:length(MBCBPS)],MERCBPS,df=DF),col=1,lwd=LWD,lty=4)
abline(h=0)

