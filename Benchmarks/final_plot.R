
par(mfrow=c(2,2))

CEX <- 1.5

################
# OVER TIME
################


    setwd("/Volumes/biostatistik$/STAFF/Michele Santacatterina/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over time periods/data2/")
    #setwd("~/KI/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over time periods/data2/")
    
    
    load("Longdata3X_L_L_TIME_2_10.RData")
    
    df <- MA_TIME_TEMP
    
    write.csv(df,"oneKOW")
    
    
    load("Longdata3X_L_L_TIME_2_6_CBPS.RData")
    
    df <- MA_TIME_TEMP
    
    write.csv(df,"oneF")
    
    
    load("Longdata3X_L_L_TIME_2_10_CBPS_Approx.RData")
    
    dfA <- MA_TIME_TEMPA
    
    write.csv(dfA,"oneA")
    
    oneK <- read.csv("oneKOW")[,2:10]
    oneF <- read.csv("oneF")[,2:6]
    oneA <- read.csv("oneA")[,2:10]
    
    mK <- c(colMeans(oneK))
    mF <- c(colMeans(oneF))
    mA <- c(colMeans(oneA))
    
    par(mar=c(4,4,0,4))
    layout(matrix(c(1,2,3,4,1,5,3,6),ncol=2),heights=c(1,2,1,2))
    plot.new()
    plot(seq(1,9),mK,type="l",ylim=c(0,55),lty=1,lwd=3,ylab="Mean computational time",xlab="Time periods",xaxt="n",
         cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
    axis(1,at=1:9,c("2","3","4","5","6","7","8","9","10"))
    lines(colMeans(ddKOW),lwd=3)
    lines(colMeans(dd),lty=3,lwd=3)
    lines(mF,lty=2,lwd=3)
    lines(mA,lty=3,lwd=3)
    
    plot.new()
    text(0.5,0.5," ",cex=2,font=2)
    plot(seq(1,9),mK,type="l",ylim=c(0,55),lty=1,lwd=3,ylab="Mean computational time",xlab="Time periods",xaxt="n",
         cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX)
    
################
# OVER COVARIATES
################
    
    
    setwd("/Volumes/biostatistik$/STAFF/Michele Santacatterina/Research/projects/Cornell/Longitudinal/data")
    #setwd("~/KI/projects/Cornell/Longitudinal/data/")
    
    
    load("Longdata3X_L_L_TIME_REP5_VAR3_8.RData") 
    
    df <- MA_TIME_TEMP
    
    write.csv(df,"oneKOW")
    
    load("Longdata3X_L_L_TIME_REP5_VAR3_8_CBPS.RData")
    
    df <- MA_TIME_TEMP
    dfA <- MA_TIME_TEMPA
    
    write.csv(df,"oneF")
    write.csv(dfA,"oneA")
    
    oneK <- read.csv("oneKOW")[,2:7]
    oneF <- read.csv("oneF")[,2:7]
    oneA <- read.csv("oneA")[,2:7]
    
    mK <- c(colMeans(oneK))
    mF <- c(colMeans(oneF))
    mA <- c(colMeans(oneA))
    
   
    plot(seq(1,6),mK,type="l",ylim=c(0,max(22)),lty=1,lwd=3,ylab="Mean computational time",xlab=expression(paste("Number of covariates, ", X[t])),xaxt="n",
         cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX,main="")
    axis(1,at=1:6,c("3","4","5","6","7","8"))
    lines(mF,lty=2,lwd=3)
    lines(mA,lty=3,lwd=3)
    
    plot(seq(1,6),mK,type="l",ylim=c(0,max(22)),lty=1,lwd=3,ylab="Mean computational time",xlab=expression(paste("Number of covariates, ", X[t])),xaxt="n",
         cex.lab=CEX, cex.axis=CEX, cex.main=CEX, cex.sub=CEX,cex=CEX,main="")
    axis(1,at=1:6,c("3","4","5","6","7","8"))
    lines(mF,lty=2,lwd=3)
    lines(mA,lty=3,lwd=3)
    
    
    
