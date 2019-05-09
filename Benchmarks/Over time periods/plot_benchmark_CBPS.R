

setwd("/Users/micsan/Desktop/data2/")
load("Longdata3X_L_L_TIME_2_6_CBPS.RData")

  df <- data.frame(MA_TIME_TEMPA)
  
  write.csv(df,file="oneA")
  df <- data.frame(MA_TIME_TEMP)

  write.csv(df,file="one")

load("/Users/micsan/Desktop/data2/Longdata3X_L_L_TIME_2_6.RData")

  df <- data.frame(MA_TIME_TEMP)

  write.csv(df,file="oneKOW")


load("Longdata3X_L_L_TIME_6_10_CBPS_onlyapprox.RData")

  df <- data.frame(MA_TIME_TEMPA)
  write.csv(df,file="twoA")

  
load("/Users/micsan/Desktop/data2/Longdata3X_L_L_TIME_7_10.RData")  

  df <- data.frame(MA_TIME_TEMP)
  write.csv(df,file="twoKOW")

  
one <- read.csv("one")
one <- one[,2:6]
  
dd <- cbind(one)

oneA <- read.csv("oneA")
oneA <- oneA[,2:6]
twoA <- read.csv("twoA")
twoA <- twoA[2:5]

ddA <- cbind(oneA,twoA)
colnames(ddA) <- c(2:10)

oneKOW <- read.csv("oneKOW")
oneKOW <- oneKOW[,2:6]
twoKOW <- read.csv("twoKOW")
twoKOW <- twoKOW[2:5]

ddKOW <- cbind(oneKOW,twoKOW)
colnames(ddKOW) <- c(2:10)


ds <- data.frame(t(seq(1,9)))
colnames(ds) <- c("2","3","4","5","6","7","8","9","10")

plot(seq(1,9),colMeans(ddA),ylim=c(0,45),type="l",lty=2,lwd=3,ylab="Mean time in seconds",xlab="Time periods",xaxt="n")
axis(1,at=1:9,c("2","3","4","5","6","7","8","9","10"))
lines(colMeans(ddKOW),lwd=3)
lines(colMeans(dd),lty=3,lwd=3)