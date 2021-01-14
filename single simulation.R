# The following code is the simulation process
# To shown the process, need to set a value firstly. 
# simulation based on R version 3.5.3

###### setting

C_ini= 0 #the amplitude of the seasonal oscillation, range from 0 to 1, C=0 means no seasonality.
sp= 4    #choose case species #1.anglerfish 2.lizardfish 3.yellow croaker 4.greenling
k= 1     #choose a seed to make the result repeatable
setwd("D:/R") #set working direction


###### required functions  ######
library(TropFishR);library(pracma);library(nlstools);library(ggplot2)
source('operating model.R')
#soVBGF
sogf<-function(t,Linf,K,t0,C,ts){
  Linf*(1-exp(-K*(t-t0)-(C*K/(2*pi))*(sin(2*pi*(t-ts))-sin(2*pi*(t0-ts)))))
}
#VBGF
vbgf<-function(t,Linf,K,t0){
  Linf*(1-exp(-K*(t-t0)))
}
#length-weight relationship
lengthweightfun<-function(L){
  a*L^b      }
#mature curve
maturefun<-function(l){1/(1+exp(-rm*(l-lm)))}
#selectivity curve
selectfun<-function(l){1/(1+exp(-r*(l-l50)))}
#yield per recruit for VBGF
ypr_vbgf<-function(f){
  YPR=0
  for(i in 1:as.numeric(length(l_ypr))){
    YPR1<-0
    for(j in 1:(i-1)){
      if(i==1){YPR1=0}
      else{YPR1<-((s_ypr[j]*f+M)*derT_ypr[j]+YPR1)}
    }
    YPR<-(w_ypr[i]*s_ypr[i]*f*(1-exp(-(s_ypr[i]*f+M)*derT_ypr[i]))*
            exp(-YPR1)/(s_ypr[i]*f+M)+YPR)
  }
  YPR
}
#spwaning stock biomass per recruit for VBGF
sbpr_vbgf<-function(f){
  SBPR=0
  for(i in 1:as.numeric(length(l_ypr))){
    YPR1<-0
    for(j in 1:(i-1)){
      if(i==1){YPR1=0}
      else{YPR1<-((s_ypr[j]*f+M)*derT_ypr[j]+YPR1)}
    }
    SBPR<-(w_ypr[i]*ma[i]*exp(-YPR1)+SBPR)
  }
  SBPR
}
#yield per recruit for soVBGF
ypr_sogf<-function(f){
  YPR=0
  for(i in 1:as.numeric(length(l_ypr))){
    YPR1<-0
    for(j in 1:(i-1)){
      if(i==1){YPR1=0}
      else{YPR1<-((s_ypr[j]*f+M)*derT_ypr[j]+YPR1)}
    }
    YPR<-(w_ypr[i]*s_ypr[i]*f*(1-exp(-(s_ypr[i]*f+M)*derT_ypr[i]))*
            exp(-YPR1)/(s_ypr[i]*f+M)+YPR)
  }
  YPR
}
#spwaning stock biomass per recruit for soVBGF
sbpr_sogf<-function(f){
  SBPR=0
  for(i in 1:as.numeric(length(l_ypr))){
    YPR1<-0
    for(j in 1:(i-1)){
      if(i==1){YPR1=0}
      else{YPR1<-((s_ypr[j]*f+M)*derT_ypr[j]+YPR1)}
    }
    SBPR<-(w_ypr[i]*ma[i]*exp(-YPR1)+SBPR)
  }
  SBPR
}


# create storage
S1<-data.frame(C_ini=rep(NA,100),linf_I=rep(NA,100),K_I=rep(NA,100),
               t0_I=rep(NA,100),f0.1_I=rep(NA,100),fmax_I=rep(NA,100),
               f40_I=rep(NA,100),f20_I=rep(NA,100),linf_II=rep(NA,100),
               K_II=rep(NA,100),t0_II=rep(NA,100),C_II=rep(NA,100),
               ts_II=rep(NA,100),f0.1_ii=rep(NA,100),fmax_ii=rep(NA,100),
               f40_ii=rep(NA,100),f20_ii=rep(NA,100),f0.1_II=rep(NA,100),
               fmax_II=rep(NA,100),f40_II=rep(NA,100),f20_II=rep(NA,100))

S2<-data.frame(C_ini=rep(NA,100),linf_III=rep(NA,100),K_III=rep(NA,100),
               t0_III=rep(NA,100),f0.1_III=rep(NA,100),fmax_III=rep(NA,100),
               f40_III=rep(NA,100),f20_III=rep(NA,100),linf_IV=rep(NA,100),
               K_IV=rep(NA,100),ta_IV=rep(NA,100),C_IV=rep(NA,100),
               ts_IV=rep(NA,100),f0.1_iv=rep(NA,100),fmax_iv=rep(NA,100),
               f40_iv=rep(NA,100),f20_iv=rep(NA,100),f0.1_IV=rep(NA,100),
               fmax_IV=rep(NA,100),f40_IV=rep(NA,100),f20_IV=rep(NA,100))



#####  input parameter  #####

biodata=read.csv("biodata.csv", header = TRUE,row.names = 1) 
S1[k,1]=C_ini
S2[k,1]=C_ini
C=C_ini
#length-weight
a=biodata[1,sp]/1000;b=biodata[2,sp]
#growth
Linf=biodata[3,sp]
K=biodata[4,sp]
t0=biodata[5,sp]
lr=biodata[11,sp]#length at recruitment
##mortality
M<-biodata[6,sp]
Z=biodata[12,sp]
##maturity
lm<-biodata[9,sp]
rm<-biodata[10,sp]

wmat=-(log(1/3)-log(3))/rm
##selectivity
r<-biodata[8,sp]
l50<-biodata[7,sp]

wqs=-(log(1/3)-log(3))/r
ts<-1/12  # 7 means growth fast
if(ts<0){ts=ts+1}else{ts=ts}


######## operating model #######

set.seed(k)
res=lfqSim(tincr=1/12,K.mu = K,Linf.mu = Linf,ts=ts,C=C,LWa = a,LWb = b,Lmat = lm,wmat = wmat,
           repro_wt = c(0,0,0,0,0,0,0,0,0,0,1,0),M=M,harvest_rate = Z-M,L50 = l50,wqs=wqs,
           bin.size = 1,timemin.date = as.Date("1990-01-01"),fished_t = seq(15+1/12,30,1/12),
           progressBar=F,timemax = 30)
a=biodata[1,sp]

####### estimation  

nls.control(maxiter = 200,warnOnly = T)

#  generate length frequency data 
catch=cbind(res$lfqbin$catch[,"27.8333333333333"],res$lfqbin$catch[,"28.0833333333333"],
            res$lfqbin$catch[,"28.3333333333333"],res$lfqbin$catch[,"28.5833333333333"])
date <- c("2017-10-30","2018-1-30","2018-4-30","2018-7-30")
date <- as.Date(date)
midlength=res$lfqbin$midLengths
elefan_data<-list(dates=date,midLengths=midlength,catch=catch)

#########  scenario III  ############

# growth estimation
elefan_vbgf=ELEFAN_GA(elefan_data,MA=5,run = 13)
plot(elefan_vbgf,main="scenario III",date.format = "%m",date.axis = "modern",xlab = "month")
Linf=elefan_vbgf$par$Linf
K=elefan_vbgf$par$K
ta=elefan_vbgf$par$t_anchor

# BRPs estimation
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
derT_ypr<-log((Linf-l_ypr)/(Linf-l_ypr-d))/K
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)

#f0.1
ypr0.1<-0.1*fderiv(ypr_vbgf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_vbgf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_vbgf(f))
fmax=f[which.max(y$ypr_vbgf.f.)]
#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_vbgf(f))
x.abs<-abs(x$sbpr_vbgf.f.-0.4*sbpr_vbgf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_vbgf.f.-0.2*sbpr_vbgf(0))
f20<-f[which.min(y.abs)]


S2[k,2]=Linf  # save 
S2[k,3]=K
S2[k,4]=ta
S2[k,5]=f0.1
S2[k,6]=fmax
S2[k,7]=f40
S2[k,8]=f20


####### scenario iv ##########

#growth estimation
elefan_sogf=ELEFAN_GA(elefan_data,MA=5,seasonalised = T,run = 13)
plot(elefan_sogf,main="scenario IV",date.format = "%m",date.axis = "modern",xlab = "month")

Linf=elefan_sogf$par$Linf
K=elefan_sogf$par$K
ta=elefan_sogf$par$t_anchor
C=elefan_sogf$par$C
ts=elefan_sogf$par$ts-ta
if(ts<0){ts=ts+1}else{ts=ts}
t0=0

#BRPs
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
derT_ypr<-log((Linf-l_ypr)/(Linf-l_ypr-d))/K
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)
#f0.1
ypr0.1<-0.1*fderiv(ypr_vbgf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_vbgf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_vbgf(f))
fmax=f[which.max(y$ypr_vbgf.f.)]
#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_vbgf(f))
x.abs<-abs(x$sbpr_vbgf.f.-0.4*sbpr_vbgf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_vbgf.f.-0.2*sbpr_vbgf(0))
f20<-f[which.min(y.abs)]


S2[k,9]=Linf  # save 
S2[k,10]=K
S2[k,11]=ta
S2[k,12]=C
S2[k,13]=ts
S2[k,14]=f0.1
S2[k,15]=fmax
S2[k,16]=f40
S2[k,17]=f20

####  scenario IV ######
# scenario IV don't need growth estimation because the growth parameters is the same as scenario iv
# BRPs
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)
T<-rep(0,(length(l_ypr)+1))
for (i in 1:length(l_ypr)) {
  T[i]=fzero(function(t) sogf(t,Linf,K,t0,C,ts)-l_ypr[i],c(-5,100))$x
}
T[length(T)]=fzero(function(t) sogf(t,Linf,K,t0,C,ts)-(l_ypr[i]+d),c(-5,100))$x
for(i in 1:length(l_ypr)){derT_ypr[i]<-T[i+1]-T[i]}# calculate dert

#f0.1
ypr0.1<-0.1*fderiv(ypr_sogf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_sogf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_sogf(f))
fmax=f[which.max(y$ypr_sogf.f.)]

#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_sogf(f))
x.abs<-abs(x$sbpr_sogf.f.-0.4*sbpr_sogf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_sogf.f.-0.2*sbpr_sogf(0))
f20<-f[which.min(y.abs)]


S2[k,18]=f0.1 #save 
S2[k,19]=fmax
S2[k,20]=f40
S2[k,21]=f20


#####  generate length-at-age data
## sample 400 inds
length_age<-rbind(res$LA$`27.8333333333333`,res$LA$`28.0833333333333`,res$LA$`28.3333333333333`,
                  res$LA$`28.5833333333333`)
set.seed(k)
samp = sample(nrow(length_age),400)
dataS1 = length_age[samp,]
dataS1 = as.data.frame(dataS1)
colnames(dataS1) = c("lt","t")

###### scenario ii ###########
#growth
set.seed(k)
output_sogf<-nls(
  lt~sogf(t,Linf,K,t0,C,ts),
  data = dataS1,
  start = list(Linf=38,K=0.36,t0=-0.25,C=0.5,ts=0.25),
  algorithm = "port",
  lower=list(Linf=0,K=0,t0=-Inf,C=0,ts=0),
  upper=list(Linf=Inf,K=Inf,t0=Inf,C=1,ts=1)
)
Linf<-as.numeric(coef(output_sogf)[1])
K<-as.numeric(coef(output_sogf)[2])
t0<-as.numeric(coef(output_sogf)[3])
C<-as.numeric(coef(output_sogf)[4])
ts<-as.numeric(coef(output_sogf)[5])
ggplot(dataS1,aes(x=t,y=lt))+
  geom_point()+
  geom_line(aes(x=t,y=sogf(t=t,Linf = Linf,K=K,t0=t0,C=C,ts=ts)),col="red")+
  theme_bw()+labs(x="age",y="length",title="scenario II")

#BRPs
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
derT_ypr<-log((Linf-l_ypr)/(Linf-l_ypr-d))/K
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)

#f0.1
ypr0.1<-0.1*fderiv(ypr_vbgf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_vbgf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_vbgf(f))
fmax=f[which.max(y$ypr_vbgf.f.)]
#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_vbgf(f))
x.abs<-abs(x$sbpr_vbgf.f.-0.4*sbpr_vbgf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_vbgf.f.-0.2*sbpr_vbgf(0))
f20<-f[which.min(y.abs)]


S1[k,9]=Linf  # save 
S1[k,10]=K
S1[k,11]=t0
S1[k,12]=C
S1[k,13]=ts
S1[k,14]=f0.1
S1[k,15]=fmax
S1[k,16]=f40
S1[k,17]=f20

##### scenario II ###########
# scenario II also don't need growth estimation because the growth parameters is the same as scenario ii
#BRPs
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)
T<-rep(0,(length(l_ypr)+1))
for (i in 1:length(l_ypr)) {
  T[i]=fzero(function(t) sogf(t,Linf,K,t0,C,ts)-l_ypr[i],c(-5,100))$x
}
T[length(T)]=fzero(function(t) sogf(t,Linf,K,t0,C,ts)-(l_ypr[i]+d),c(-5,100))$x
for(i in 1:length(l_ypr)){derT_ypr[i]<-T[i+1]-T[i]}# calculate dert

#f0.1
ypr0.1<-0.1*fderiv(ypr_sogf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_sogf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_sogf(f))
fmax=f[which.max(y$ypr_sogf.f.)]
#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_sogf(f))
x.abs<-abs(x$sbpr_sogf.f.-0.4*sbpr_sogf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_sogf.f.-0.2*sbpr_sogf(0))
f20<-f[which.min(y.abs)]


S1[k,18]=f0.1  #save 
S1[k,19]=fmax
S1[k,20]=f40
S1[k,21]=f20

####### scenario I ############
#growth
set.seed(k)
output_vbgf<-nls(
  lt~vbgf(t,Linf,K,t0),
  data = dataS1,
  start = list(Linf=38,K=0.36,t0=-0.25),
  algorithm = "port",
  lower=list(Linf=0,K=0,t0=-Inf),
  upper=list(Linf=Inf,K=Inf,t0=Inf)
)
Linf<-as.numeric(coef(output_vbgf)[1])#output
K<-as.numeric(coef(output_vbgf)[2])
t0<-as.numeric(coef(output_vbgf)[3])

ggplot(dataS1,aes(x=t,y=lt))+
  geom_point()+
  geom_line(aes(x=t,y=vbgf(t=t,Linf = Linf,K=K,t0=t0)),col="red")+
  theme_bw()+labs(x="age",y="length",title="scenario I")
#BRPs
d=1 #interval=1cm
l_ypr<-seq(lr,(floor(Linf)-3*d),d)
derT_ypr<-log((Linf-l_ypr)/(Linf-l_ypr-d))/K
w_ypr<-lengthweightfun(l_ypr)
s_ypr<-selectfun(l_ypr)

#f0.1
ypr0.1<-0.1*fderiv(ypr_vbgf,0,1)
f<-seq(0,3,0.001)
yprderiv<-fderiv(ypr_vbgf,f,1)
x=data.frame(f,yprderiv)
x.abs<-abs(yprderiv-ypr0.1)
f0.1=f[which.min(x.abs)]
#fmax
y=data.frame(f,ypr_vbgf(f))
fmax=f[which.max(y$ypr_vbgf.f.)]
#f40%
ma<-maturefun(l_ypr)
x<-data.frame(f,sbpr_vbgf(f))
x.abs<-abs(x$sbpr_vbgf.f.-0.4*sbpr_vbgf(0))
f40<-f[which.min(x.abs)]
#f20%
y.abs<-abs(x$sbpr_vbgf.f.-0.2*sbpr_vbgf(0))
f20<-f[which.min(y.abs)]


S1[k,2]=Linf  # save 
S1[k,3]=K
S1[k,4]=t0
S1[k,5]=f0.1
S1[k,6]=fmax
S1[k,7]=f40
S1[k,8]=f20

