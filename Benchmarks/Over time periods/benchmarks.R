


#LINEAR

close(matlab)
source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/simu_3Xs_291117_LL_linPoly1_nn_prova_bTIME.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6.RData")

MA_TIME_TEMP
MA_TIME_GPML
MA_TIME_GUROBI

TIME_TEMP
TIME_GPML
TIME_GUROBI

boxplot(MA_TIME_TEMP)


source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/simu_3Xs_291117_LL_linPoly1_nn_prova_CBPS_TIME.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_L_L_TIME_2_6_CBPS.RData")

MA_TIME_TEMP
MA_TIME_TEMPA

boxplot(MA_TIME_TEMP)
boxplot(MA_TIME_TEMPA)

TIME_TEMP
TIME_TEMPA

plot(TIME_TEMP,type="l")




#NON-LINEAR

close(matlab)
source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/simu_3Xs_291117_NLNL_linPoly2_nn_prova_TIME.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_NL_NL_TIME_2_6.RData")

MA_TIME_TEMP
MA_TIME_GPML
MA_TIME_GUROBI

boxplot(MA_TIME_TEMP)

TIME_TEMP
TIME_GPML
TIME_GUROBI



source("/home/micsan/research/Research/projects/Cornell/Longitudinal/Rcode/Final for submission/Benchmarks/simu_3Xs_291117_NLNL_linPoly2_nn_prova_CBPS_TIME.R")

load("/home/micsan/research/Research/projects/Cornell/Longitudinal/data/Longdata3X_NL_NL_TIME_2_6_CBPS.RData")

MA_TIME_TEMP
MA_TIME_TEMP[which(MA_TIME_TEMP>200)] <- mean(MA_TIME_TEMP[,4])
MA_TIME_TEMPA

df <- rbind(MA_TIME_TEMP,MA_TIME_TEMPA)
df <- cbind(df,c(rep(1,100),rep(2,100)))

boxplot(MA_TIME_TEMP,add=T)
boxplot(MA_TIME_TEMPA)


TIME_TEMP
TIME_TEMPA

plot(seqqtt[1:length(TIME_TEMP)],TIME_TEMP,type="l")
lines(TIME_TEMPA,col=2)
