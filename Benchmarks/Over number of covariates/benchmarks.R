

#########################################################################################################
#########################################################################################################v
#########################################################################################################v
#KOW vs CBPS FULL

close(matlab)
source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over number of covariates/simu_3Xs_291117_LL_linPoly1_nn_prova_bVAR.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_8.RData") 

MA_TIME_TEMP
MA_TIME_GPML
MA_TIME_GUROBI

TIME_TEMP
TIME_GPML
TIME_GUROBI

boxplot(MA_TIME_TEMP)



source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over number of covariates/simu_3Xs_291117_LL_linPoly1_nn_prova_CBPS_VAR.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_8_CBPS.RData")

MA_TIME_TEMP

boxplot(MA_TIME_TEMP)

TIME_TEMP

plot(colMeans(MA_TIME_TEMP),type="l",ylim=c(0,max(MA_TIME_TEMP)))





#########################################################################################################v
#########################################################################################################v
#########################################################################################################v
#KOW vs CBPS FULL

close(matlab)
source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over number of covariates/simu_3Xs_291117_LL_linPoly1_nn_prova_bVAR.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_6.RData") 

MA_TIME_TEMP
MA_TIME_GPML
MA_TIME_GUROBI

TIME_TEMP
TIME_GPML
TIME_GUROBI

boxplot(MA_TIME_TEMP)



source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/Over number of covariates/simu_3Xs_291117_LL_linPoly1_nn_prova_CBPS_VAR.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_6_CBPS.RData")

MA_TIME_TEMP

boxplot(MA_TIME_TEMP)

TIME_TEMP
TIME_TEMPA

plot(colMeans(MA_TIME_TEMP),type="l",ylim=c(0,max(MA_TIME_TEMP)))
lines(colMeans(MA_TIME_TEMPA),lty=2)



# 10.041 10.385 10.944 10.671
# 2.3794  4.7188  7.0668 10.6472


load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_8.RData") 

df <- MA_TIME_TEMP

write.csv(df,"oneKOW")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_REP5_VAR3_8_CBPS.RData")

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

plot(mK,type="l",ylim=c(0,max(22)))
lines(mF,lty=2)
lines(mA,lty=3)





