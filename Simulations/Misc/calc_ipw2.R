


####################################################################
#
# LINEAR
#

calc_ipw_1 <- function(data){
  
  a1 <-as.integer(data$Treatment[data$time==1]==1)
  a0 <-as.integer(data$Treatment[data$time==1]==0)
  
  r11 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==1))
  r10 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==0))
  r01 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==1))
  r00 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==0))
  
  r111 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r110 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r100 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  r011 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r001 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r101 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r010 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r000 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  
  ####IPTW
  #1/pa0|x0
  
  Z1 <- (data$X1[data$time==1])
  Z2 <- (data$X2[data$time==1])
  Z3 <- (data$X3[data$time==1])
  
  ZZ1 <- (data$X1[data$time==2])
  ZZ2 <- (data$X2[data$time==2])
  ZZ3 <- (data$X3[data$time==2])
  
  ZZZ1 <- (data$X1[data$time==3])
  ZZZ2 <- (data$X2[data$time==3])
  ZZZ3 <- (data$X3[data$time==3])
  
  
  ppr0 <- glm(data$Treatment[data$time==1]~
                Z1
              +Z2
              +Z3
              
              ,family="binomial")$fitted
  ipw0<-rep(0,n)
  ipw0[which(a1==1)] <- 1/ppr0[which(a1==1)]
  ipw0[which(a0==1)] <- 1/(1-ppr0[which(a0==1)])
  ipw1 <-1
  
  ppr1 <- glm(data$Treatment[data$time==2]~
                data$Treatment[data$time==1]
              +Z1
              +Z2
              +Z3
              +ZZ1
              +ZZ2
              +ZZ3
              +data$Treatment[data$time==1]:Z1
              +data$Treatment[data$time==1]:Z2
              +data$Treatment[data$time==1]:Z3
              
              +data$Treatment[data$time==1]:ZZ1
              +data$Treatment[data$time==1]:ZZ2
              +data$Treatment[data$time==1]:ZZ3
              
              ,family="binomial")$fitted
  
  ipw1<-rep(0,n)
  ipw1[which(r11==1)] <- 1/ppr1[which(r11==1)]
  ipw1[which(r10==1)] <- 1/(1-ppr1[which(r10==1)])
  ipw1[which(r01==1)] <- 1/ppr1[which(r01==1)]
  ipw1[which(r00==1)] <- 1/(1-ppr1[which(r00==1)])
  
  ppr2 <- glm(data$Treatment[data$time==3]~
                data$Treatment[data$time==2]
              +data$Treatment[data$time==1]
              +Z1
              +Z2
              +Z3
              +ZZ1
              +ZZ2
              +ZZ3
              +ZZZ1
              +ZZZ2
              +ZZZ3
              
              +data$Treatment[data$time==1]:Z1
              +data$Treatment[data$time==1]:Z2
              +data$Treatment[data$time==1]:Z3
              
              +data$Treatment[data$time==1]:ZZ1
              +data$Treatment[data$time==1]:ZZ2
              +data$Treatment[data$time==1]:ZZ3
              
              +data$Treatment[data$time==1]:ZZZ1
              +data$Treatment[data$time==1]:ZZZ2
              +data$Treatment[data$time==1]:ZZZ3
              
              ,family="binomial")$fitted
  
  ipw2<-rep(0,n)
  ipw2[which(r111==1)] <- 1/ppr2[which(r111==1)]
  ipw2[which(r110==1)] <- 1/(1-ppr2[which(r110==1)])
  ipw2[which(r100==1)] <- 1/(1-ppr2[which(r100==1)])
  ipw2[which(r011==1)] <- 1/(ppr2[which(r011==1)])
  ipw2[which(r001==1)] <- 1/ppr2[which(r001==1)]
  ipw2[which(r101==1)] <- 1/(ppr2[which(r101==1)])
  ipw2[which(r010==1)] <- 1/(1-ppr2[which(r010==1)])
  ipw2[which(r000==1)] <- 1/(1-ppr2[which(r000==1)])
  
  
  
  
  
  ########################################################STABLE IPTW
  #1/pa0|x0
  
  sppr1 <- glm(data$Treatment[data$time==2]~
                 data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw1<-rep(0,n)
  sipw1[which(r11==1)] <- sppr1[which(r11==1)]/ppr1[which(r11==1)]
  sipw1[which(r10==1)] <- (1-sppr1[which(r10==1)])/(1-ppr1[which(r10==1)])
  sipw1[which(r01==1)] <- sppr1[which(r01==1)]/ppr1[which(r01==1)]
  sipw1[which(r00==1)] <- (1-sppr1[which(r00==1)])/(1-ppr1[which(r00==1)])
  
  sppr2 <- glm(data$Treatment[data$time==3]~
                 data$Treatment[data$time==2]
               +data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw2<-rep(0,n)
  sipw2[which(r111==1)] <- sppr2[which(r111==1)]/ppr2[which(r111==1)]
  sipw2[which(r110==1)] <- (1-sppr2[which(r110==1)])/(1-ppr2[which(r110==1)])
  sipw2[which(r100==1)] <- (1-sppr2[which(r100==1)])/(1-ppr2[which(r100==1)])
  sipw2[which(r011==1)] <- (sppr2[which(r011==1)])/(ppr2[which(r011==1)])
  sipw2[which(r001==1)] <- sppr2[which(r001==1)]/ppr2[which(r001==1)]
  sipw2[which(r101==1)] <- (sppr2[which(r101==1)])/(ppr2[which(r101==1)])
  sipw2[which(r010==1)] <- (1-sppr2[which(r010==1)])/(1-ppr2[which(r010==1)])
  sipw2[which(r000==1)] <- (1-sppr2[which(r000==1)])/(1-ppr2[which(r000==1)])
  
  
  ipw <- ipw0*ipw1*ipw2
  sipw <- ipw0*sipw1*sipw2
  
  return(list(ipw_1=ipw,sipw_1=sipw))
  
}



####################################################################
#
# QUADRATIC
#

calc_ipw_2 <- function(data){
  
  a1 <-as.integer(data$Treatment[data$time==1]==1)
  a0 <-as.integer(data$Treatment[data$time==1]==0)
  
  r11 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==1))
  r10 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==0))
  r01 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==1))
  r00 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==0))
  
  r111 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r110 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r100 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  r011 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r001 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r101 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r010 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r000 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  
  ####IPTW
  #1/pa0|x0
  
  Z1 <- (data$X1[data$time==1])
  Z2 <- (data$X2[data$time==1])
  Z3 <- (data$X3[data$time==1])
  
  ZZ1 <- (data$X1[data$time==2])
  ZZ2 <- (data$X2[data$time==2])
  ZZ3 <- (data$X3[data$time==2])
  
  ZZZ1 <- (data$X1[data$time==3])
  ZZZ2 <- (data$X2[data$time==3])
  ZZZ3 <- (data$X3[data$time==3])
  
  
  two_Z1Z2 <- 2*(Z1)*(Z2)
  two_Z1Z3 <- 2*(Z1)*(Z3)
  two_Z2Z3 <- 2*(Z2)*(Z3)
  
  two_Z1ZZ1 <- 2*(Z1)*(ZZ1)
  two_Z1ZZ2 <- 2*(Z1)*(ZZ2)
  two_Z1ZZ3 <- 2*(Z1)*(ZZ3)
  
  two_Z2ZZ1 <- 2*(Z2)*(ZZ1)
  two_Z2ZZ2 <- 2*(Z2)*(ZZ2)
  two_Z2ZZ3 <- 2*(Z2)*(ZZ3)
  
  two_Z3ZZ1 <- 2*(Z3)*(ZZ1)
  two_Z3ZZ2 <- 2*(Z3)*(ZZ2)
  two_Z3ZZ3 <- 2*(Z3)*(ZZ3)
  
  two_Z1Z3 <- 2*(Z1)*(Z3)
  two_Z2Z3 <- 2*(Z2)*(Z3)
  
  two_ZZ1ZZ2 <- 2*(ZZ1)*(ZZ2)
  two_ZZ1ZZ3 <- 2*(ZZ1)*(ZZ3)
  two_ZZ2ZZ3 <- 2*(ZZ2)*(ZZ3)
  
  two_ZZZ1ZZZ2 <- 2*(ZZZ1)*(ZZZ2)
  two_ZZZ1ZZZ3 <- 2*(ZZZ1)*(ZZZ3)
  two_ZZZ2ZZZ3 <- 2*(ZZZ2)*(ZZZ3)
  
  two_ZZ1ZZZ1 <- 2*(ZZ1)*(ZZZ1)
  two_ZZ1ZZZ2 <- 2*(ZZ1)*(ZZZ2)
  two_ZZ1ZZZ3 <- 2*(ZZ1)*(ZZZ3)
  
  two_ZZ2ZZZ1 <- 2*(ZZ2)*(ZZZ1)
  two_ZZ2ZZZ2 <- 2*(ZZ2)*(ZZZ2)
  two_ZZ2ZZZ3 <- 2*(ZZ2)*(ZZZ3)
  
  two_ZZ3ZZZ1 <- 2*(ZZ3)*(ZZZ1)
  two_ZZ3ZZZ2 <- 2*(ZZ3)*(ZZZ2)
  two_ZZ3ZZZ3 <- 2*(ZZ3)*(ZZZ3)
  
  two_Z1ZZZ1 <- 2*(Z1)*(ZZZ1)
  two_Z1ZZZ2 <- 2*(Z1)*(ZZZ2)
  two_Z1ZZZ3 <- 2*(Z1)*(ZZZ3)
  
  two_Z2ZZZ1 <- 2*(Z2)*(ZZZ1)
  two_Z2ZZZ2 <- 2*(Z2)*(ZZZ2)
  two_Z2ZZZ3 <- 2*(Z2)*(ZZZ3)
  
  two_Z3ZZZ1 <- 2*(Z3)*(ZZZ1)
  two_Z3ZZZ2 <- 2*(Z3)*(ZZZ2)
  two_Z3ZZZ3 <- 2*(Z3)*(ZZZ3)
  
  
  ppr0 <- glm(data$Treatment[data$time==1]~
                I(Z1^2)
              +I(Z2^2)
              +I(Z3^2)
              +Z1:Z2
              +Z1:Z3
              +Z2:Z3
              
              ,family="binomial")$fitted
  ipw0<-rep(0,n)
  ipw0[which(a1==1)] <- 1/ppr0[which(a1==1)]
  ipw0[which(a0==1)] <- 1/(1-ppr0[which(a0==1)])
  ipw1 <-1
  
  ppr1 <- glm(data$Treatment[data$time==2]~
                I(Z1^2)
              +I(Z2^2)
              +I(Z3^2)
              +I(ZZ1^2)
              +I(ZZ2^2)
              +I(ZZ3^2)
              +Z1
              +Z1
              +Z2
              +ZZ1
              +ZZ2
              +ZZ3
                
              +  data$Treatment[data$time==1]:I(Z1^2)
              + data$Treatment[data$time==1]:I(Z2^2)
              + data$Treatment[data$time==1]:I(Z3^2)
              
              + data$Treatment[data$time==1]:Z1
              + data$Treatment[data$time==1]:Z2
              + data$Treatment[data$time==1]:Z3
              
              + data$Treatment[data$time==1]:Z1:Z2
              + data$Treatment[data$time==1]:Z1:Z3
              + data$Treatment[data$time==1]:Z2:Z3
              
              + data$Treatment[data$time==1]:I(ZZ1^2)
              + data$Treatment[data$time==1]:I(ZZ2^2)
              + data$Treatment[data$time==1]:I(ZZ3^2)
              
              + data$Treatment[data$time==1]:ZZ1
              + data$Treatment[data$time==1]:ZZ2
              + data$Treatment[data$time==1]:ZZ3
              
              + data$Treatment[data$time==1]:ZZ1:ZZ2
              + data$Treatment[data$time==1]:ZZ1:ZZ3
              + data$Treatment[data$time==1]:ZZ2:ZZ3
              
              + data$Treatment[data$time==1]:Z1:ZZ1
              + data$Treatment[data$time==1]:Z1:ZZ2
              + data$Treatment[data$time==1]:Z1:ZZ3
              
              + data$Treatment[data$time==1]:Z2:ZZ1
              + data$Treatment[data$time==1]:Z2:ZZ2
              + data$Treatment[data$time==1]:Z2:ZZ3
              
              + data$Treatment[data$time==1]:Z3:ZZ1
              + data$Treatment[data$time==1]:Z3:ZZ2
              + data$Treatment[data$time==1]:Z3:ZZ3
              
              ,family="binomial")$fitted
  
  ipw1<-rep(0,n)
  ipw1[which(r11==1)] <- 1/ppr1[which(r11==1)]
  ipw1[which(r10==1)] <- 1/(1-ppr1[which(r10==1)])
  ipw1[which(r01==1)] <- 1/ppr1[which(r01==1)]
  ipw1[which(r00==1)] <- 1/(1-ppr1[which(r00==1)])
  
  ppr2 <- glm(data$Treatment[data$time==3]~
                data$Treatment[data$time==1]
              + data$Treatment[data$time==2]
              +  I(ZZZ1^2)
              + I(ZZZ1^2)
              + I(ZZZ1^2)
              +Z1
              +Z2
              +Z3
              +ZZ1
              +ZZ2
              +ZZ3
              +ZZZ1
              +ZZZ2
              +ZZZ3
              
              + data$Treatment[data$time==1]:Z1
              + data$Treatment[data$time==1]:Z2
              + data$Treatment[data$time==1]:Z3
              
              + data$Treatment[data$time==1]:ZZ1
              + data$Treatment[data$time==1]:ZZ2
              + data$Treatment[data$time==1]:ZZ3
              
              + data$Treatment[data$time==1]:ZZZ1
              + data$Treatment[data$time==1]:ZZZ2
              + data$Treatment[data$time==1]:ZZZ3
              
              +  data$Treatment[data$time==1]:I(ZZZ1^2)
              + data$Treatment[data$time==1]:I(ZZZ2^2)
              + data$Treatment[data$time==1]:I(ZZZ3^2)
              + data$Treatment[data$time==1]:ZZZ1:ZZZ2
              + data$Treatment[data$time==1]:ZZZ1:ZZZ3
              + data$Treatment[data$time==1]:ZZZ2:ZZZ3
              
              +  data$Treatment[data$time==1]:I(Z1^2)
              + data$Treatment[data$time==1]:I(Z2^2)
              + data$Treatment[data$time==1]:I(Z3^2)
              + data$Treatment[data$time==1]:Z1:Z2
              + data$Treatment[data$time==1]:Z1:Z3
              + data$Treatment[data$time==1]:Z2:Z3
              
              + data$Treatment[data$time==1]:I(ZZ1^2)
              + data$Treatment[data$time==1]:I(ZZ2^2)
              + data$Treatment[data$time==1]:I(ZZ3^2)
              + data$Treatment[data$time==1]:ZZ1:ZZ2
              + data$Treatment[data$time==1]:ZZ1:ZZ3
              + data$Treatment[data$time==1]:ZZ2:ZZ3
              
              + data$Treatment[data$time==1]:Z1:ZZ1
              + data$Treatment[data$time==1]:Z1:ZZ2
              + data$Treatment[data$time==1]:Z1:ZZ3
              
              + data$Treatment[data$time==1]:Z2:ZZ1
              + data$Treatment[data$time==1]:Z2:ZZ2
              + data$Treatment[data$time==1]:Z2:ZZ3
              
              + data$Treatment[data$time==1]:Z3:ZZ1
              + data$Treatment[data$time==1]:Z3:ZZ2
              + data$Treatment[data$time==1]:Z3:ZZ3
              
              + data$Treatment[data$time==1]:Z1:ZZZ1
              + data$Treatment[data$time==1]:Z1:ZZZ2
              + data$Treatment[data$time==1]:Z1:ZZZ3
              
              + data$Treatment[data$time==1]:Z2:ZZZ1
              + data$Treatment[data$time==1]:Z2:ZZZ2
              + data$Treatment[data$time==1]:Z2:ZZZ3
              
              + data$Treatment[data$time==1]:Z3:ZZZ1
              + data$Treatment[data$time==1]:Z3:ZZZ2
              + data$Treatment[data$time==1]:Z3:ZZZ3
              
              + data$Treatment[data$time==1]:ZZ1:ZZZ1
              + data$Treatment[data$time==1]:ZZ1:ZZZ2
              + data$Treatment[data$time==1]:ZZ1:ZZZ3
              
              + data$Treatment[data$time==1]:ZZ2:ZZZ1
              + data$Treatment[data$time==1]:ZZ2:ZZZ2
              + data$Treatment[data$time==1]:ZZ2:ZZZ3
              
              + data$Treatment[data$time==1]:ZZ3:ZZZ1
              + data$Treatment[data$time==1]:ZZ3:ZZZ2
              + data$Treatment[data$time==1]:ZZ3:ZZZ3
              
              
              + data$Treatment[data$time==2]:I(ZZZ1^2)
              + data$Treatment[data$time==2]:I(ZZZ2^2)
              + data$Treatment[data$time==2]:I(ZZZ3^2)
              + data$Treatment[data$time==2]:ZZZ1:ZZZ2
              + data$Treatment[data$time==2]:ZZZ1:ZZZ3
              + data$Treatment[data$time==2]:ZZZ2:ZZZ3
              
              + data$Treatment[data$time==2]:I(ZZZ1^2)
              + data$Treatment[data$time==2]:I(ZZZ2^2)
              + data$Treatment[data$time==2]:I(ZZZ3^2)
              + data$Treatment[data$time==2]:ZZZ1:ZZZ2
              + data$Treatment[data$time==2]:ZZZ1:ZZZ3
              + data$Treatment[data$time==2]:ZZZ2:ZZZ3
              
              +  data$Treatment[data$time==2]:I(Z1^2)
              + data$Treatment[data$time==2]:I(Z2^2)
              + data$Treatment[data$time==2]:I(Z3^2)
              + data$Treatment[data$time==2]:Z1:Z2
              + data$Treatment[data$time==2]:Z1:Z3
              + data$Treatment[data$time==2]:Z2:Z3
              
              + data$Treatment[data$time==2]:I(ZZ1^2)
              + data$Treatment[data$time==2]:I(ZZ2^2)
              + data$Treatment[data$time==2]:I(ZZ3^2)
              + data$Treatment[data$time==2]:ZZ1:ZZ2
              + data$Treatment[data$time==2]:ZZ1:ZZ3
              + data$Treatment[data$time==2]:ZZ2:ZZ3
              
              + data$Treatment[data$time==2]:Z1:ZZ1
              + data$Treatment[data$time==2]:Z1:ZZ2
              + data$Treatment[data$time==2]:Z1:ZZ3
              
              + data$Treatment[data$time==2]:Z2:ZZ1
              + data$Treatment[data$time==2]:Z2:ZZ2
              + data$Treatment[data$time==2]:Z2:ZZ3
              
              + data$Treatment[data$time==2]:Z3:ZZ1
              + data$Treatment[data$time==2]:Z3:ZZ2
              + data$Treatment[data$time==2]:Z3:ZZ3
              
              + data$Treatment[data$time==2]:Z1:ZZZ1
              + data$Treatment[data$time==2]:Z1:ZZZ2
              + data$Treatment[data$time==2]:Z1:ZZZ3
              
              + data$Treatment[data$time==2]:Z2:ZZZ1
              + data$Treatment[data$time==2]:Z2:ZZZ2
              + data$Treatment[data$time==2]:Z2:ZZZ3
              
              + data$Treatment[data$time==2]:Z3:ZZZ1
              + data$Treatment[data$time==2]:Z3:ZZZ2
              + data$Treatment[data$time==2]:Z3:ZZZ3
              
              + data$Treatment[data$time==2]:ZZ1:ZZZ1
              + data$Treatment[data$time==2]:ZZ1:ZZZ2
              + data$Treatment[data$time==2]:ZZ1:ZZZ3
              
              + data$Treatment[data$time==2]:ZZ2:ZZZ1
              + data$Treatment[data$time==2]:ZZ2:ZZZ2
              + data$Treatment[data$time==2]:ZZ2:ZZZ3
              
              + data$Treatment[data$time==2]:ZZ3:ZZZ1
              + data$Treatment[data$time==2]:ZZ3:ZZZ2
              + data$Treatment[data$time==2]:ZZ3:ZZZ3
              
              ,family="binomial")$fitted
  
  ipw2<-rep(0,n)
  ipw2[which(r111==1)] <- 1/ppr2[which(r111==1)]
  ipw2[which(r110==1)] <- 1/(1-ppr2[which(r110==1)])
  ipw2[which(r100==1)] <- 1/(1-ppr2[which(r100==1)])
  ipw2[which(r011==1)] <- 1/(ppr2[which(r011==1)])
  ipw2[which(r001==1)] <- 1/ppr2[which(r001==1)]
  ipw2[which(r101==1)] <- 1/(ppr2[which(r101==1)])
  ipw2[which(r010==1)] <- 1/(1-ppr2[which(r010==1)])
  ipw2[which(r000==1)] <- 1/(1-ppr2[which(r000==1)])
  
  
  
  
  
  ########################################################STABLE IPTW
  #1/pa0|x0
  
  sppr1 <- glm(data$Treatment[data$time==2]~
                 data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw1<-rep(0,n)
  sipw1[which(r11==1)] <- sppr1[which(r11==1)]/ppr1[which(r11==1)]
  sipw1[which(r10==1)] <- (1-sppr1[which(r10==1)])/(1-ppr1[which(r10==1)])
  sipw1[which(r01==1)] <- sppr1[which(r01==1)]/ppr1[which(r01==1)]
  sipw1[which(r00==1)] <- (1-sppr1[which(r00==1)])/(1-ppr1[which(r00==1)])
  
  sppr2 <- glm(data$Treatment[data$time==3]~
                 data$Treatment[data$time==2]
               +data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw2<-rep(0,n)
  sipw2[which(r111==1)] <- sppr2[which(r111==1)]/ppr2[which(r111==1)]
  sipw2[which(r110==1)] <- (1-sppr2[which(r110==1)])/(1-ppr2[which(r110==1)])
  sipw2[which(r100==1)] <- (1-sppr2[which(r100==1)])/(1-ppr2[which(r100==1)])
  sipw2[which(r011==1)] <- (sppr2[which(r011==1)])/(ppr2[which(r011==1)])
  sipw2[which(r001==1)] <- sppr2[which(r001==1)]/ppr2[which(r001==1)]
  sipw2[which(r101==1)] <- (sppr2[which(r101==1)])/(ppr2[which(r101==1)])
  sipw2[which(r010==1)] <- (1-sppr2[which(r010==1)])/(1-ppr2[which(r010==1)])
  sipw2[which(r000==1)] <- (1-sppr2[which(r000==1)])/(1-ppr2[which(r000==1)])
  
  
  ipw <- ipw0*ipw1*ipw2
  sipw <- ipw0*sipw1*sipw2
  
  return(list(ipw_2=ipw,sipw_2=sipw))
  
}



####################################################################
#
# CUBIC
#

calc_ipw_3 <- function(data){
  
  a1 <-as.integer(data$Treatment[data$time==1]==1)
  a0 <-as.integer(data$Treatment[data$time==1]==0)
  
  r11 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==1))
  r10 <-as.integer((data$Treatment[data$time==(1)]==1)&(data$Treatment[data$time==2]==0))
  r01 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==1))
  r00 <-as.integer((data$Treatment[data$time==(1)]==0)&(data$Treatment[data$time==2]==0))
  
  r111 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r110 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r100 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  r011 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==1))
  r001 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r101 <-as.integer((data$Treatment[data$time==(1)]==1)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==1))
  r010 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==1&(data$Treatment[data$time==3]==0))
  r000 <-as.integer((data$Treatment[data$time==(1)]==0)&data$Treatment[data$time==(2)]==0&(data$Treatment[data$time==3]==0))
  
  ####IPTW
  #1/pa0|x0
  
  Z1 <- (data$X1[data$time==1])
  Z2 <- (data$X2[data$time==1])
  Z3 <- (data$X3[data$time==1])
  
  ZZ1 <- (data$X1[data$time==2])
  ZZ2 <- (data$X2[data$time==2])
  ZZ3 <- (data$X3[data$time==2])
  
  ZZZ1 <- (data$X1[data$time==3])
  ZZZ2 <- (data$X2[data$time==3])
  ZZZ3 <- (data$X3[data$time==3])
  
  
  two_Z1Z2 <- 2*(Z1)*(Z2)
  two_Z1Z3 <- 2*(Z1)*(Z3)
  two_Z2Z3 <- 2*(Z2)*(Z3)
  
  two_ZZ1ZZ2 <- 2*(ZZ1)*(ZZ2)
  two_ZZ1ZZ3 <- 2*(ZZ1)*(ZZ3)
  two_ZZ2ZZ3 <- 2*(ZZ2)*(ZZ3)
  
  two_ZZZ1ZZZ2 <- 2*(ZZZ1)*(ZZZ2)
  two_ZZZ1ZZZ3 <- 2*(ZZZ1)*(ZZZ3)
  two_ZZZ2ZZZ3 <- 2*(ZZZ2)*(ZZZ3)
  
  ppr0 <- glm(data$Treatment[data$time==1]~
                I(Z1^2)
              +I(Z2^2)
              +I(Z3^2)
              +two_Z1Z2
              +two_Z1Z3
              +two_Z2Z3
              
              ,family="binomial")$fitted
  ipw0<-rep(0,n)
  ipw0[which(a1==1)] <- 1/ppr0[which(a1==1)]
  ipw0[which(a0==1)] <- 1/(1-ppr0[which(a0==1)])
  ipw1 <-1
  
  ppr1 <- glm(data$Treatment[data$time==2]~
                data$Treatment[data$time==1]
              + I(Z1^2)
              +I(Z2^2)
              +I(Z3^2)
              +two_Z1Z2
              +two_Z1Z3
              +two_Z2Z3
              
              + I(ZZ1^2)
              +I(ZZ2^2)
              +I(ZZ3^2)
              +two_ZZ1ZZ2
              +two_ZZ1ZZ3
              +two_ZZ2ZZ3
              
              ,family="binomial")$fitted
  
  ipw1<-rep(0,n)
  ipw1[which(r11==1)] <- 1/ppr1[which(r11==1)]
  ipw1[which(r10==1)] <- 1/(1-ppr1[which(r10==1)])
  ipw1[which(r01==1)] <- 1/ppr1[which(r01==1)]
  ipw1[which(r00==1)] <- 1/(1-ppr1[which(r00==1)])
  
  ppr2 <- glm(data$Treatment[data$time==3]~
                data$Treatment[data$time==2]
              +data$Treatment[data$time==1]
              + I(Z1^2)
              +I(Z2^2)
              +I(Z3^2)
              +two_Z1Z2
              +two_Z1Z3
              +two_Z2Z3
              
              + I(ZZ1^2)
              +I(ZZ2^2)
              +I(ZZ3^2)
              +two_ZZ1ZZ2
              +two_ZZ1ZZ3
              +two_ZZ2ZZ3
              
              + I(ZZZ1^2)
              +I(ZZZ2^2)
              +I(ZZZ3^2)
              +two_ZZZ1ZZZ2
              +two_ZZZ1ZZZ3
              +two_ZZZ2ZZZ3
              
              ,family="binomial")$fitted
  
  ipw2<-rep(0,n)
  ipw2[which(r111==1)] <- 1/ppr2[which(r111==1)]
  ipw2[which(r110==1)] <- 1/(1-ppr2[which(r110==1)])
  ipw2[which(r100==1)] <- 1/(1-ppr2[which(r100==1)])
  ipw2[which(r011==1)] <- 1/(ppr2[which(r011==1)])
  ipw2[which(r001==1)] <- 1/ppr2[which(r001==1)]
  ipw2[which(r101==1)] <- 1/(ppr2[which(r101==1)])
  ipw2[which(r010==1)] <- 1/(1-ppr2[which(r010==1)])
  ipw2[which(r000==1)] <- 1/(1-ppr2[which(r000==1)])
  
  
  
  
  
  ########################################################STABLE IPTW
  #1/pa0|x0
  
  sppr1 <- glm(data$Treatment[data$time==2]~
                 data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw1<-rep(0,n)
  sipw1[which(r11==1)] <- sppr1[which(r11==1)]/ppr1[which(r11==1)]
  sipw1[which(r10==1)] <- (1-sppr1[which(r10==1)])/(1-ppr1[which(r10==1)])
  sipw1[which(r01==1)] <- sppr1[which(r01==1)]/ppr1[which(r01==1)]
  sipw1[which(r00==1)] <- (1-sppr1[which(r00==1)])/(1-ppr1[which(r00==1)])
  
  sppr2 <- glm(data$Treatment[data$time==3]~
                 data$Treatment[data$time==2]
               +data$Treatment[data$time==1]
               ,family="binomial")$fitted
  
  sipw2<-rep(0,n)
  sipw2[which(r111==1)] <- sppr2[which(r111==1)]/ppr2[which(r111==1)]
  sipw2[which(r110==1)] <- (1-sppr2[which(r110==1)])/(1-ppr2[which(r110==1)])
  sipw2[which(r100==1)] <- (1-sppr2[which(r100==1)])/(1-ppr2[which(r100==1)])
  sipw2[which(r011==1)] <- (sppr2[which(r011==1)])/(ppr2[which(r011==1)])
  sipw2[which(r001==1)] <- sppr2[which(r001==1)]/ppr2[which(r001==1)]
  sipw2[which(r101==1)] <- (sppr2[which(r101==1)])/(ppr2[which(r101==1)])
  sipw2[which(r010==1)] <- (1-sppr2[which(r010==1)])/(1-ppr2[which(r010==1)])
  sipw2[which(r000==1)] <- (1-sppr2[which(r000==1)])/(1-ppr2[which(r000==1)])
  
  
  ipw <- ipw0*ipw1*ipw2
  sipw <- ipw0*sipw1*sipw2
  
  return(list(ipw=ipw,sipw=sipw))
  
}