library("PMA")
library("ggplot2")
library("cowplot")
library("gridExtra")

################################################################################
# generate simulated data 1
data1 = function(seed0=1){
  n = 50
  p = 80
  q = 100
  
  u = c(rep(1,10),  rep(-1,10), rep(0,p-20))
  v = c(rep(-1,15), rep(1,15),  rep(0,q-30))
  w = c(rep(1,25),  rep(0,n-25))
  
  SNR = 0.1
  gamma1 = sqrt(norm(u%*%t(v),"F")^2/(SNR*n*p))
  gamma2 = sqrt(norm(w%*%t(v),"F")^2/(SNR*n*q))
  
  set.seed(seed0)
  X = w%*%t(u) + gamma1*matrix(rnorm(n*p),ncol=p)
  Y = w%*%t(v) + gamma2*matrix(rnorm(n*q),ncol=q)
  
  ku = 20
  kv = 30
  kw = 25
  
  data = list(X=X,Y=Y,u=u,v=v,w=w, ku=ku, kv=kv, kw=kw)
}

################################################################################
# generate simulated data 2
data2 = function(seed0=1){
  n = 50*2
  p = 800
  q = 1000
  
  u = c(rep(1,100),  rep(-1,100), rep(0,p-200))
  v = c(rep(-1,150), rep(1,150),  rep(0,q-300))
  w = c(rep(1,50),  rep(0,n-50))
  
  SNR = 0.1
  gamma1 = sqrt(norm(u%*%t(v),"F")^2/(SNR*n*p))
  gamma2 = sqrt(norm(w%*%t(v),"F")^2/(SNR*n*q))
  
  set.seed(1*seed0)
  X = w%*%t(u) + gamma1*matrix(rnorm(n*p),ncol=p)
  Y = w%*%t(v) + gamma2*matrix(rnorm(n*q),ncol=q)
  
  ku = 200
  kv = 300
  kw = 50
  
  data = list(X=X,Y=Y,u=u,v=v,w=w, ku=ku, kv=kv, kw=kw)
}

################################################################################
# generate simulated data 3
data3 = function(seed0=1){
  n = 500
  p = 8000
  q = 10000
  
  u = c(rep(1,1000),  rep(-1,1000), rep(0,p-2000))
  v = c(rep(-1,1500), rep(1,1500),  rep(0,q-3000))
  w = c(rep(1,250),  rep(0,n-250))
  
  
  SNR = 0.1
  gamma1 = sqrt(norm(u%*%t(v),"F")^2/(SNR*n*p))
  gamma2 = sqrt(norm(w%*%t(v),"F")^2/(SNR*n*q))
  
  set.seed(1*seed0)
  X = w%*%t(u) + gamma1*matrix(rnorm(n*p),ncol=p)
  Y = w%*%t(v) + gamma2*matrix(rnorm(n*q),ncol=q)
  
  ku = 2000
  kv = 3000
  kw = 250
  
  data = list(X=X,Y=Y,u=u,v=v,w=w, ku=ku, kv=kv, kw=kw)
}

################################################################################
# testing on the simulated data 
simulation = function(flag=1){
  rep_num = 10 # Repeat the experiment for 10 times
  table_list = list()
  for(jj in 1:rep_num){
    print(jj)
    ############################################################################
    # obtain simulated data
    if(flag==1){
      data = data1(jj)
      X = data$X; Y = data$Y
      ku = data$ku; kv = data$kv; kw = data$kw; 
    }
    
    if(flag==2){
      data = data2(jj)
      X = data$X; Y = data$Y
      ku = data$ku; kv = data$kv; kw = data$kw; 
    }
    
    if(flag==3){
      data = data3(jj)
      X = data$X; Y = data$Y
      ku = data$ku; kv = data$kv; kw = data$kw; 
    }
    
    ############################################################################
    # 记录程序运行的时间+++++
    time_vec = c()
    
    # 1-PLS
    library("PMA") # PLS with no sparsity
    ptm  = proc.time()
    out1 = CCA(X, Y, typex="standard", typez="standard", K=1, niter=100,
               penaltyx = 1, penaltyz = 1)
    tim  = proc.time() - ptm # 用户, CPU time
    time_vec[1] = tim[1]
    
    # 2-PMD
    library("PMA") # sCCA with L1-norm constrain
    ptm  = proc.time()
    out2 = CCA(X, Y, typex="standard", typez="standard", K=1, niter=100,
               penaltyx = 0.4, penaltyz = 0.467)
    tim  = proc.time() - ptm # 用户, CPU time
    time_vec[2] = tim[1]
    
    # 3-L0-sPLS
    source('Function/sPLS.R')
    ptm  = proc.time()
    out3 = L0_sPLS(X, Y, ku, kv, niter=100, err=10^(-5), nstart=5)
    tim  = proc.time() - ptm # 用户, CPU time
    time_vec[3] = tim[1]
    
    # 4-L2-wsPLS
    source('Function/SWCCA_CJE18.R')
    ptm  = proc.time()
    out4 = SWCCA_CJE18(X, Y, ku, kv, kw, seed0=5)
    tim  = proc.time() - ptm # 用户, CPU time
    time_vec[4] = tim[1]
    
    # 5-l_infty_wsPLS
    source('Function/PALM_wsPLS.R')
    ptm  = proc.time()
    out5 = PALM_wsPLS(X, Y, ku, kv, kw, Lc=0.1, nstart=5, seed0=2)
    tim = proc.time() - ptm # 用户, CPU time
    time_vec[5] = tim[1]
    
    # For comparison, we assume that w = 1 for these methods
    out1$w = out2$w = out3$w = rep(1,length(data$w))
    
    ############################################################################
    source("Function/Fun_TPR_TNR_ACC.R")
    Res_table = Fun_TPR_TNR_ACC_table(data, out1, out2, out3, out4, out5)
    
    Res_table$time = time_vec
    
  
    table_list[[jj]] = Res_table
    if(jj==1){
      avg_table = Res_table
    }else{
      avg_table = avg_table + Res_table
    }
  }
  
  # Final results 
  table_avg = avg_table/rep_num
  
  return(list(table_avg=table_avg, table_list=table_list))
}

################################################################################
################################################################################
# testing on the simulated data
Res_data1 = simulation(flag=1)
Res_data2 = simulation(flag=2)
Res_data3 = simulation(flag=3)

save.image("F2_results_summary.RData")