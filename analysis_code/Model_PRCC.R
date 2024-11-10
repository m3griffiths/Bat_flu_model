library(epiR)
library(tidyverse)

#step 1 - deterministic model

deterministic_model<-function(params, t_max, init_val){
  
  alpha<-params["alpha"]
  r_bd<-params["b"]
  gamma<-params["gamma"]
  delta<-params["delta"]
  beta0<-params["beta0"]
  beta1<-params["beta1"]
  beta2<-params["beta2"]
  period1<-params["period1"]
  period2<-params["period2"]

  
  s<-c()
  e<-c()
  i<-c()
  r<-c()
  
  
  s[1]<-init_val[1]
  e[1]<-init_val[2]
  i[1]<-init_val[3]
  r[1]<-init_val[4]
  
  
  for (t in 2:(t_max*30)){
    

      beta<-((1/(beta0/4*(1+beta1*sin(2*t*pi/period1))*
                   (1+beta2*sin(2*t*pi/period2))+1)))

    
    s[t] <- s[t-1] + 
      (1/r_bd)*(s[t-1]+e[t-1]+i[t-1]+r[t-1]) + 
      (1/gamma)*r[t-1] - 
      beta*
      s[t-1]*i[t-1] - 
      (1/r_bd)*s[t-1] 
    
    
    e[t] <- e[t-1] + 
      beta*
      s[t-1]*i[t-1] -
      (1/delta)*e[t-1] - 
      (1/r_bd)*e[t-1]
    
    i[t] <- i[t-1] + 
      (1/delta)*e[t-1] - 
      (1/alpha)*i[t-1]- 
      (1/r_bd)*i[t-1]
    
    r[t] <- r[t-1] + 
      (1/alpha)*i[t-1] - 
      (1/gamma)*r[t-1] - 
      (1/r_bd)*r[t-1]
    
    
  }
  
  time<-seq(1:(30*t_max))
  return(cbind(time,s,e,i,r))
}



init_val<-c(0.45,0.1,0,0.45)


params<- c(b=rnorm(n=1, mean=365*3.5, sd=365),
           alpha=rnorm(n=1, mean=7, sd=2),
           delta=rnorm(n=1, mean=2, sd=0.5),
           gamma=rnorm(n=1, mean=200, sd=50), 
           beta0=rnorm(n=1, mean=8, sd=2),
           beta1=rbeta(n=1, shape1=2.625, shape2=2.625),
           beta2=rbeta(n=1, shape1=2.625, shape2=2.625),
           period1=365,
           period2=runif(n=1, min=1200, max=1700)
)





out<-as.data.frame(deterministic_model(params, 144, init_val))

ggplot(data=out, aes(x=(time/365)+2007, y=r/(s+e+i+r)))+geom_line(col="green")+
  geom_line(aes(x=(time/365)+2007, y=e/(s+e+i+r)), col="red")+
  geom_line(aes(x=(time/365)+2007, y=i/(s+e+i+r)), col="blue")+
  ylim(0,1)+
  theme_bw()

#step 2 - run prcc

prcc_seir_df<-data.frame(matrix(ncol=4320, nrow=0))
prcc_seir_in_df<-data.frame(matrix(ncol=9, nrow=0))
colnames(prcc_seir_in_df)<-names(params)

for(n in 1:500){
  
  params<- c(b=rnorm(n=1, mean=365*3.5, sd=365/1.5),
             alpha=rnorm(n=1, mean=7, sd=1.5),
             delta=rgamma(n=1, shape=6, rate=2),
             gamma=rnorm(n=1, mean=200, sd=50), 
             beta0=rnorm(n=1, mean=8, sd=2),
             beta1=rbeta(n=1, shape1=2.625, shape2=2.625),
             beta2=rbeta(n=1, shape1=2.625, shape2=2.625),
             period1=365,
             period2=runif(n=1, min=1200, max=1700)
  )
  
  out<-as.data.frame(deterministic_model(params, 144, init_val))
  
  prcc_seir_in_df<-rbind(prcc_seir_in_df, t(as.data.frame(params)))
  prcc_seir_df<-rbind(prcc_seir_df, data.frame(t(out$r)))

}

#arrange into an appropriate dataframe
prcc_in_df<-cbind(prcc_seir_in_df, prcc_seir_df)


prcc_out<-epi.prcc(prcc_in_df[,c(1:7,9,911)], sided.test = 2, conf.level = 0.95)

ggplot(data=prcc_out)+geom_point(aes(x=(est), y=var))+theme_bw()+
  geom_vline(aes(xintercept = 0), linetype="dashed")+
  geom_errorbar(aes(xmax=upper, xmin=lower, y=var))

prcc_out

prcc_out_list<-list()

for(n in 1:4320){
  
  prcc_out<-epi.prcc(prcc_in_df[,c(1:7,9,9+n)], sided.test = 2, conf.level = 0.95)
  prcc_out_list<-append(prcc_out_list, list(prcc_out))
  
}

prcc_out<-prcc_out_list[[320]]
prcc_mean<-data.frame(var=prcc_out$var, est=prcc_out$est, lower=prcc_out$est-prcc_out$lower, upper=prcc_out$upper-prcc_out$est)


for(n in 321:4320){
  for(r in 1:nrow(prcc_mean)){
    prcc_mean[r,2]<-prcc_mean[r,2]+abs(prcc_out_list[[n]][r,2])
    prcc_mean[r,3]<-prcc_mean[r,3]+(prcc_out_list[[n]][r,2]-prcc_out_list[[n]][r,3])
    prcc_mean[r,4]<-prcc_mean[r,4]+(prcc_out_list[[n]][r,4]-prcc_out_list[[n]][r,2])
  }
}


ggplot(data=prcc_mean)+geom_point(aes(x=(est)/4000, y=var))+theme_bw()+
  geom_vline(aes(xintercept = 0), linetype="dashed")+
  geom_errorbar(aes(xmax=((est)/4000)+(upper/4000), xmin=((est)/4000)-(lower/4000), y=var))+
  labs(x="Partial rank correlation coefficient", y="Parameters")

