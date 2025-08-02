




library(MASS)
library(lavaan)
library(simsem)
library(semPlot)
library(OpenMx)
library(tidySEM)
library(sjPlot)
library(foreach)
library(doParallel)



population.model<-'
f1=~0.7*x1+0.7*x2+0.75*x3+0.8*x4+0.8*x5
f2=~0.7*x6+0.7*x7+0.75*x8+0.8*x9+0.8*x10
f3=~0.7*x11+0.7*x12+0.75*x13+0.8*x14+0.8*x15
f1~~1*f1
f2~~1*f2
f3~~1*f3
f1~~0.3*f2
f1~~0.4*f3
f2~~0.5*f3
x1~~1*x1
x2~~1*x2
x3~~1*x3
x4~~1*x4
x5~~1*x5
x6~~1*x6
x7~~1*x7
x8~~1*x8
x9~~1*x9
x10~~1*x10
x11~~1*x11
x12~~1*x12
x13~~1*x13
x14~~1*x14
x15~~1*x15
'


my.model<-'
f1=~0.7*x1+x2+x3+x4+x5
f2=~0.7*x6+x7+x8+x9+x10
f3=~0.7*x11+x12+x13+x14+x15
'


setwd('C:/Users/xin512/Desktop/R')
#load("data.RData")

u1<-list(mean=0,sd=1)      ##  or df=1
u2<-list(mean=0,sd=1)        ##  or df=1
u3<-list(mean=0,sd=1)       ##  or df=1
v<-list(mean=0,sd=1)
distname<-c(rep("norm",3))   ##  or chisq 
distname1<-c(rep("norm",15))
facdist<-bindDist(distname,u1,u2,u3)
errordist<-bindDist(distname1,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v)

##

m<-3
R<-1000


########  

lav_t_check <- function(object, verbose = FALSE) {
  
  #stopifnot(inherits(object, "lavaan"))
  lavpartable    <- object@ParTable
  lavmodel       <- object@Model
  lavdata        <- object@Data
  
  var.ov.ok <- var.lv.ok <- result.ok <- TRUE
  
  # 1a. check for negative variances ov
  var.idx <- which(lavpartable$op == "~~" &
                     lavpartable$lhs %in% lavNames(object, "ov") &
                     lavpartable$lhs == lavpartable$rhs)
  if(length(var.idx) > 0L && any(lavpartable$est[var.idx] < 0.0)) {
    result.ok <- var.ov.ok <- FALSE
    #warning("lavaan WARNING: some estimated ov variances are negative")
  }
  
  # 1b. check for negative variances lv
  var.idx <- which(lavpartable$op == "~~" &
                     lavpartable$lhs %in% lavNames(object, "lv") &
                     lavpartable$lhs == lavpartable$rhs)
  if(length(var.idx) > 0L && any(lavpartable$est[var.idx] < 0.0)) {
    result.ok <- var.lv.ok <- FALSE
    #warning("lavaan WARNING: some estimated lv variances are negative")
  }
  
  result.ok
}

###




Test.all<-function(N,m,R) {
  
  n<-R
  n1<-0
  n2<-0
  n3<-0
  n4<-0
  n5<-0
  n6<-0
  n7<-0
  n8<-0
  
  a1<-rep(NA,n)
  a2<-rep(NA,n)
  a3<-rep(NA,n)
  a4<-rep(NA,n)
  a5<-rep(NA,n)
  a6<-rep(NA,n)
  a7<-rep(NA,n)
  a8<-rep(NA,n)
  
  u<-0
  
  S<-data.frame(matrix(nrow = 25,ncol = 1))
  names(S)<-c("N")
  
  for (i in 1:n) {
    data<-generate(population.model,n=N,facDist = facdist,errorDist = errordist)
    fit<-cfa(my.model,data =data,estimator="MLM")
    fit1<-cfa(my.model,data =data,ridge = TRUE, ridge.constant = 1e-3,estimator="MLM")  ##   ridge
    # 
    
    if(fit@optim$converged){
      if(lav_t_check(fit)){
        
        df<-fitmeasures(fit,"df")
        fml_value<-fitMeasures(fit,"fmin") 
        fml_value1<-fitMeasures(fit1,"fmin")
        fml<-2*fml_value ### 
        fml_r<-2*fml_value1
        
        ###  Bartlett 
         c_b<-1-((2*ncol(data)+11)/6-2*m/3)/nrow(data)   ###   N = n
        #c_b <- (1 - (2*ncol(data) + 4*m + 11)/(6*nrow(data))) * nrow(data) / (nrow(data) - 1)  ## N-1= n
         ######################
        ## lnλ=N/2*[ln|Σ|-ln|S|+tr(SΣ⁻¹)-p]=N*fml=n*fml=N/2*T;  fml_value=T/2;
        ## Tml=N*fml_value
        ###############
        
        ###
        c_p<-ncol(data)/nrow(data)       ##c=p/n
        mu<-ncol(data)-ncol(data)*(1-1/c_p)*log(1-c_p)-0.5*log(1-c_p)
        sd_c<-sqrt(-2*log(1-c_p)-2*c_p)
        
        mu_y0<-2.015910*(m/ncol(data))+1.291412*(ncol(data)/nrow(data))-0.278377*m+0.036066*ncol(data)-2.393643   ## g
        mu_y1<-1.531018*(ncol(data)/nrow(data))-0.237301*m+0.029336*ncol(data)-2.035224                           ## g1
        
        ###  
        chi.sc<-fitmeasures(fit,"chisq.scaled")
        chi<-fitmeasures(fit,"chisq")
        chi.sc1<-fitmeasures(fit1,"chisq.scaled")
        
        
        a1[i]<-T1<-c_b*chi.sc                              ##     Trmlb
        a2[i]<-T2<-c_b*chi                                 ##     Tmlb  
        
        a3[i]<-T3<-(fml-mu)/sd_c                           #      T_F
        a4[i]<-T4<-(fml-mu-mu_y0*sd_c)/sd_c                ###    T_FC
        a5[i]<-T5<-(fml_r-mu-mu_y0*sd_c)/sd_c              ###    T_FCr  
        
        a6[i]<-T6<-(chi.sc/nrow(data)-mu)/sd_c             ##     T_CsF
        a7[i]<-T7<-(chi.sc/nrow(data)-mu-mu_y1*sd_c)/sd_c  ##     T_CsFC
        a8[i]<-T8<-(chi.sc1/nrow(data)-mu-mu_y1*sd_c)/sd_c ##     T_CsFCr
        
        
        if((1-pchisq(T1,df))<=0.05){n1<-n1+1}
        if((1-pchisq(T2,df))<=0.05){n2<-n2+1}
        
        if((1-pnorm(T3,0,1))<=0.05){n3<-n3+1}
        if((1-pnorm(T4,0,1))<=0.05){n4<-n4+1}
        if((1-pnorm(T5,0,1))<=0.05){n5<-n5+1}
        
        if((1-pnorm(T6,0,1))<=0.05){n6<-n6+1}
        if((1-pnorm(T7,0,1))<=0.05){n7<-n7+1}
        if((1-pnorm(T8,0,1))<=0.05){n8<-n8+1}
        
        u<-u+1
      }
    }
  }
  
  erra1<-n1/u
  erra2<-n2/u
  erra3<-n3/u
  erra4<-n4/u
  erra5<-n5/u
  
  erra6<-n6/u
  erra7<-n7/u
  erra8<-n8/u
  
  p<-u/n
  
  S[1]<-c(mean(a1,na.rm=TRUE),sd(a1,na.rm=TRUE),erra1,
          mean(a2,na.rm=TRUE),sd(a2,na.rm=TRUE),erra2,
          mean(a3,na.rm=TRUE),sd(a3,na.rm=TRUE),erra3,
          mean(a4,na.rm=TRUE),sd(a4,na.rm=TRUE),erra4,
          mean(a5,na.rm=TRUE),sd(a5,na.rm=TRUE),erra5,
          mean(a6,na.rm=TRUE),sd(a6,na.rm=TRUE),erra6,
          mean(a7,na.rm=TRUE),sd(a7,na.rm=TRUE),erra7,
          mean(a8,na.rm=TRUE),sd(a8,na.rm=TRUE),erra8,p)
  
  
  return(S)
}




###

N<-c(50,80,100,200,300,400,500,800,1000,2000)
L<-data.frame(matrix(nrow = 25,ncol = length(N)))
cl<-detectCores(logical = F)-1 
registerDoParallel(makeCluster(cl)) 
L<-foreach(i=1:length(N),.packages = c("simsem","lavaan"), .combine=cbind) %dopar% Test.all(N[i],m,R)
stopImplicitCluster()

names(L)<-c("50","80","100","200","300","400","500","800","1000","2000")
tab_df(L,digits = 3)








