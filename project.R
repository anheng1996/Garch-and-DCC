####Garch(1,1) Calculation of vol
require(gdata)
library(foreign)
interest=0.0169

stock_price<-read.xls("~/R/FRM/Components.xlsx")
return_matrix<-matrix(nrow = 250,ncol = 30)

garch_coef<-matrix(nrow = 4,ncol = 30)

for (i in 1:30) {
  return_matrix[,i]<-diff(log(stock_price[,i]))
  fun1<-function(x) {
    sigmasqhat=rep(0,length(return_matrix[,i]))
    sigmasqhat[1]=x[4]^2
    if (x[1]+x[2]>=1 ||x[1]<0||x[2]<0 ||x[3]<0.001 ||x[4]<0.001) {
      NeglogLH=9999
    } else {
      for (j in 1:(length(return_matrix[,i])-1)) {
        sigmasqhat[j+1]=(1-x[1]-x[2])*x[3]^2+x[1]*return_matrix[j,i]^2+x[2]*sigmasqhat[j]
      }
      f<-(1/(sqrt(2*pi*sigmasqhat)))*exp(-0.5*return_matrix[,i]^2/sigmasqhat)
      NeglogLH=-sum(log(f))
    }
    
    return(NeglogLH)
  }
  guess<-c(0.1,0.8,0.01,0.01)
  output_1=optim(guess,fun1)
  
  garch_coef[1,i]=output_1$par[1]  #alpha
  garch_coef[2,i]=output_1$par[2]  #beta
  garch_coef[3,i]=output_1$par[3]  #sigma
  garch_coef[4,i]=output_1$par[4]  #sigma1
}

####Calculation of sigma 1 to sigma 250
stock_vol<-matrix(nrow = 250, ncol = 30)
for (i in 1:30) {
  stock_vol[1,i]=garch_coef[4,i]
  for (j in 2:250) {
    stock_vol[j,i]=sqrt((1-garch_coef[1,i]-garch_coef[2,i])*garch_coef[3,i]^2+garch_coef[1,i]*return_matrix[j-1,i]^2+garch_coef[2,i]*stock_vol[(j-1),i]^2)
  }
}

std_matrix<-matrix(0,nrow = 30,ncol = 30)   ####the vols at t+1 are on the diagnoal
for (i in 1:30) {
  vari=(1-garch_coef[1,i]-garch_coef[2,i])*garch_coef[3,i]^2+garch_coef[1,i]*return_matrix[250,i]^2+garch_coef[2,i]*stock_vol[250,i]^2
  std_matrix[i,i]=sqrt(vari)
}

####standard return
standard_return<-return_matrix/stock_vol

####DCC(1,1)  Calculation of correlation
z_matrix<-array(0,dim = c(30,30,250))

for (i in 1:250) {
  for (j in 1:30) {
    for (k in j:30) {
      z_matrix[j,k,i]=standard_return[i,j]*standard_return[i,k]
    }
  }
}

fun2<-function(x) {
  if (x[1]+x[2]>=1 ||x[1]<0||x[2]<0) {
    NeglogLH=9999999999
  } else {
    q_matrix<-array(0,dim = c(30,30,250))
    for (i in 1:30) {
      q_matrix[i,i,1]=1
    }
    for (i in 1:29){
      for (j in (i+1):30) {
        q_matrix[i,j,1]=mean(z_matrix[i,j,1:250])
      }
    }
    for (i in 2:250) {
      for (j in 1:30) {
        for (k in j:30) {
          q_matrix[j,k,i]=q_matrix[j,k,1]+x[1]*(z_matrix[j,k,(i-1)]-q_matrix[j,k,1])+x[2]*(q_matrix[j,k,(i-1)]-q_matrix[j,k,1])
        }
      }
    }
    r_matrix<-array(0,dim = c(30,30,250))
    for (i in 1:250) {
      for (j in 1:29) {
        for (k in (j+1):30) {
          r_matrix[j,k,i]=q_matrix[j,k,i]/sqrt(q_matrix[j,j,i]*q_matrix[k,k,i])
        }
      }
    }
    logf_matrix<-array(0,dim = c(30,30,250))
    for (i in 1:250) {
      for (j in 1:29) {
        for (k in (j+1):30) {
          logf_matrix[j,k,i]=-log(sqrt(1-r_matrix[j,k,i]^2))-0.5*(z_matrix[j,j,i]+z_matrix[k,k,i]-2*r_matrix[j,k,i]*z_matrix[j,k,i])/(1-r_matrix[j,k,i]^2)
        }
      }
    }
    NeglogLH=-sum(logf_matrix)
  }
  return(NeglogLH)
}

guess2<-c(0.2,0.7)
output_2=optim(guess2,fun2)
alpha_rho=output_2$par[1]
beta_rho=output_2$par[2]

####simulation of rho
q_matrix<-array(0,dim = c(30,30,250))
for (i in 1:30) {
  q_matrix[i,i,1]=1
}
for (i in 1:29){
  for (j in (i+1):30) {
    q_matrix[i,j,1]=mean(z_matrix[i,j,1:250])
  }
}
for (i in 2:250) {
  for (j in 1:30) {
    for (k in j:30) {
      q_matrix[j,k,i]=q_matrix[j,k,1]+alpha_rho*(z_matrix[j,k,(i-1)]-q_matrix[j,k,1])+beta_rho*(q_matrix[j,k,(i-1)]-q_matrix[j,k,1])
    }
  }
}

r_matrix<-matrix(nrow = 30,ncol = 30)  ####correlation matrix at t+1
for (i in 1:30) {
  r_matrix[i,i]=1
}

for (i in 1:29){
  for (j in (i+1):30) {
    qij=q_matrix[i,j,1]+alpha_rho*(z_matrix[i,j,250]-q_matrix[i,j,1])+beta_rho*(q_matrix[i,j,250]-q_matrix[i,j,1])
    qii=1+alpha_rho*(z_matrix[i,i,250]-1)+beta_rho*(q_matrix[i,i,250]-1)
    qjj=1+alpha_rho*(z_matrix[j,j,250]-1)+beta_rho*(q_matrix[j,j,250]-1)
    r_matrix[i,j]=qij/sqrt(qii*qjj)
    r_matrix[j,i]=r_matrix[i,j]
  }
}

cov_matrix<-std_matrix%*%r_matrix%*%std_matrix

####weight
weight<-read.xls("~/R/FRM/weights.xlsx")[,2]

####dividend yield
yield<-read.xls("~/R/FRM/div yield.xlsx")[,2]

####Calculation about DJI
dji_price<-read.xls("~/R/FRM/DJI.xlsx")[,1]
return_dji<-diff(log(dji_price))
dow_divisor=sum(stock_price[251,1:30])/dji_price[251]

fun3<-function(x) {
  sigmasqhat=rep(0,length(return_dji))
  sigmasqhat[1]=x[4]^2
  if (x[1]+x[2]>=1 ||x[1]<0||x[2]<0 ||x[3]<0.001 ||x[4]<0.001) {
    NeglogLH=9999
  } else {
    for (j in 1:(length(return_dji)-1)) {
      sigmasqhat[j+1]=(1-x[1]-x[2])*x[3]^2+x[1]*return_dji[j]^2+x[2]*sigmasqhat[j]
    }
  }
  f<-(1/(sqrt(2*pi*sigmasqhat)))*exp(-0.5*return_dji^2/sigmasqhat)
  NeglogLH=-sum(log(f))
  return(NeglogLH)
}
guess3<-c(0.1,0.8,0.01,0.01)
output_3=optim(guess3,fun3)

alpha_dji=output_3$par[1]  #alpha
beta_dji=output_3$par[2]  #beta
sigma_dji=output_3$par[3]  #sigma
sigma1_dji=output_3$par[4]  #sigma1

dji_vol<-rep(0,251)
dji_vol[1]=sigma1_dji
for (j in 1:249) {
  dji_vol[j+1]=sqrt((1-alpha_dji-beta_dji)*sigma_dji^2+alpha_dji*return_dji[j]^2+beta_dji*dji_vol[j]^2)
}
dji_vol[251]=sqrt((1-alpha_dji-beta_dji)*sigma_dji^2+alpha_dji*return_dji[250]^2+beta_dji*dji_vol[250]^2)

####vega neutral
library(fOptions)
tau=0.25
vega_strike_dji=-100*GBSGreeks(Selection = "vega",TypeFlag = "c",S=dji_price[251]/100,X=dji_price[251]/100,
                          Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252))-
  100*GBSGreeks(Selection = "vega",TypeFlag = "p",S=dji_price[251]/100,X=dji_price[251]/100,
            Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252))

vega_strike_stock<-matrix(nrow = 30,ncol = 3) ##vega of call is in the 1st col, vega of put is in the 2nd col, vega of strike is in the 3rd row)
for (i in 1:30){
  vega_strike_stock[i,1]=GBSGreeks(Selection = "vega",TypeFlag = "c",S=stock_price[251,i],X=stock_price[251,i],
                                   Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252))
  vega_strike_stock[i,2]=GBSGreeks(Selection = "vega",TypeFlag = "p",S=stock_price[251,i],X=stock_price[251,i],
                                   Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252))
  vega_strike_stock[i,3]=vega_strike_stock[i,1]+vega_strike_stock[i,2]
}

total_vega_stock=0
for (i in 1:30) {
  total_vega_stock=total_vega_stock+weight[i]*vega_strike_stock[i,3]
}

k=-vega_strike_dji/total_vega_stock
kw=k*weight

####option prices
call_price<-rep(0,31)
put_price<-rep(0,31)
library(fOptions)

for (i in 1:30) {
  call_price[i]=GBSOption(TypeFlag = "c",S=stock_price[251,i],X=stock_price[251,i],
                          Time=tau,r=interest,b=interest-yield[i],sigma =stock_vol[250,i]*sqrt(252))@price
}
call_price[31]=100*GBSOption(TypeFlag = "c",S=dji_price[251]/100,X=dji_price[251]/100,
                             Time=tau,r=interest,b=interest-yield[31],sigma =dji_vol[250]*sqrt(252))@price

for (i in 1:30) {
  put_price[i]=GBSOption(TypeFlag = "p",S=stock_price[251,i],X=stock_price[251,i],
                         Time=tau,r=interest,b=interest-yield[i],sigma =stock_vol[250,i]*sqrt(252))@price
}
put_price[31]=100*GBSOption(TypeFlag = "p",S=dji_price[251]/100,X=dji_price[251]/100,
                            Time=tau,r=interest,b=interest-yield[31],sigma =dji_vol[250]*sqrt(252))@price

####today's value of my portfolio
value0=as.numeric(-call_price[31]-put_price[31])
for (i in 1:30) {
  value0=value0+as.numeric(kw[i]*(call_price[i]+put_price[i])) 
}

####Calculation of VaR
library(MASS)
set.seed(123)
n_sim=10000  ####number of simulation
return_sim<-mvrnorm(n=n_sim,mu=rep(0,len=30),Sigma=cov_matrix)

stock_sim<-matrix(nrow = n_sim,ncol = 30)
dji_sim<-rep(0,n_sim)

for (i in 1:30) {
  for (j in 1:n_sim) {
    stock_sim[j,i]=exp(return_sim[j,i])*stock_price[251,i]
  }
}

for (i in 1:n_sim) {
  dji_sim[i]=sum(stock_sim[i,1:30])/dow_divisor
}

call_stock_sim<-matrix(nrow = n_sim,ncol = 30)
for (i in 1:n_sim) {
  for (j in 1:30) {
    call_stock_sim[i,j]=GBSOption(TypeFlag = "c",S=stock_sim[i,j],X=stock_price[251,j],
                        Time=tau-1/252,r=interest,b=interest-yield[j],sigma =std_matrix[j,j]*sqrt(252))@price
  }
}

put_stock_sim<-matrix(nrow = n_sim,ncol = 30)
for (i in 1:n_sim) {
  for (j in 1:30) {
    put_stock_sim[i,j]=GBSOption(TypeFlag = "p",S=stock_sim[i,j],X=stock_price[251,j],
                                  Time=tau-1/252,r=interest,b=interest-yield[j],sigma =std_matrix[j,j]*sqrt(252))@price
  }
}

call_dji_sim<-rep(0,n_sim)
for (i in 1:n_sim) {
  call_dji_sim[i]=100*GBSOption(TypeFlag = "c",S=dji_sim[i]/100,X=dji_price[251]/100,
                Time=tau-1/252,r=interest,b=interest-yield[31],sigma = dji_vol[251]*sqrt(252))@price
}

put_dji_sim<-rep(0,n_sim)
for (i in 1:n_sim) {
  put_dji_sim[i]=100*GBSOption(TypeFlag = "p",S=dji_sim[i]/100,X=dji_price[251]/100,
                                Time=tau-1/252,r=interest,b=interest-yield[31],sigma = dji_vol[251]*sqrt(252))@price
}

value1<-rep(0,n_sim)   ####possible portfolio value at t+1
for (i in 1:n_sim) {
  value1[i]=-call_dji_sim[i]-put_dji_sim[i]+sum(kw*call_stock_sim[i,])+sum(kw*put_stock_sim[i,])
}

return_portfolio=value1-rep(value0,n_sim)
VaR_0.95=-quantile(return_portfolio,0.05)   ####95% VaR
VaR_0.99=-quantile(return_portfolio,0.01)   ####95% VaR

####Expeced Shortfall
loss=-return_portfolio
ES_0.95=mean(loss[loss>=VaR_0.95])
ES_0.99=mean(loss[loss>=VaR_0.99])

####Calculation of greeks
#vega,delta
vega=delta=0

#gamma
gamma_matrix<-matrix(nrow = 31,ncol = 3)
for (i in 1:30) {
  gamma_matrix[i,1]=GBSGreeks(Selection = "Gamma",TypeFlag = "c",S=stock_price[251,i],X=stock_price[251,i],
                              Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  gamma_matrix[i,2]=GBSGreeks(Selection = "Gamma",TypeFlag = "p",S=stock_price[251,i],X=stock_price[251,i],
                              Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  gamma_matrix[i,3]=kw[i]*(gamma_matrix[i,1]+gamma_matrix[i,2])
}

gamma_matrix[31,1]=100*GBSGreeks(Selection = "Gamma",TypeFlag = "c",S=dji_price[251]/100,X=dji_price[251]/100,
                                 Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
gamma_matrix[31,2]=100*GBSGreeks(Selection = "Gamma",TypeFlag = "p",S=dji_price[251]/100,X=dji_price[251]/100,
                                 Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
gamma_matrix[31,3]=-(gamma_matrix[31,1]+gamma_matrix[31,2])

gamma=sum(gamma_matrix[,3])

#theta
theta_matrix<-matrix(nrow = 31,ncol = 3)
for (i in 1:30) {
  theta_matrix[i,1]=GBSGreeks(Selection = "theta",TypeFlag = "c",S=stock_price[251,i],X=stock_price[251,i],
                              Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  theta_matrix[i,2]=GBSGreeks(Selection = "theta",TypeFlag = "p",S=stock_price[251,i],X=stock_price[251,i],
                              Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  theta_matrix[i,3]=kw[i]*(theta_matrix[i,1]+theta_matrix[i,2])
}

theta_matrix[31,1]=100*GBSGreeks(Selection = "theta",TypeFlag = "c",S=dji_price[251]/100,X=dji_price[251]/100,
                                 Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
theta_matrix[31,2]=100*GBSGreeks(Selection = "theta",TypeFlag = "p",S=dji_price[251]/100,X=dji_price[251]/100,
                                 Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
theta_matrix[31,3]=-(theta_matrix[31,1]+theta_matrix[31,2])

theta=sum(theta_matrix[,3])

#rho
rho_matrix<-matrix(nrow = 31,ncol = 3)
for (i in 1:30) {
  rho_matrix[i,1]=GBSGreeks(Selection = "rho",TypeFlag = "c",S=stock_price[251,i],X=stock_price[251,i],
                            Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  rho_matrix[i,2]=GBSGreeks(Selection = "rho",TypeFlag = "p",S=stock_price[251,i],X=stock_price[251,i],
                            Time=tau,r=interest,b=interest-yield[i],sigma = stock_vol[250,i]*sqrt(252) )
  rho_matrix[i,3]=kw[i]*(rho_matrix[i,1]+rho_matrix[i,2])
}

rho_matrix[31,1]=100*GBSGreeks(Selection = "rho",TypeFlag = "c",S=dji_price[251]/100,X=dji_price[251]/100,
                               Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
rho_matrix[31,2]=100*GBSGreeks(Selection = "rho",TypeFlag = "p",S=dji_price[251]/100,X=dji_price[251]/100,
                               Time=tau,r=interest,b=interest-yield[31],sigma = dji_vol[250]*sqrt(252) )
rho_matrix[31,3]=-(rho_matrix[31,1]+rho_matrix[31,2])

rho=sum(rho_matrix[,3])