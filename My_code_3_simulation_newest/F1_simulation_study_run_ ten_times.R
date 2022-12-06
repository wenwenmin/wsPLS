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
    
    if(jj==1){
      avg_table = Res_table
    }else{
      avg_table = avg_table + Res_table
    }
  }
  
  # Final results 
  avg_table_data1 = avg_table/rep_num
  
  return(avg_table_data1)
}

################################################################################
################################################################################
# testing on the simulated data
avg_table_data1 = simulation(flag=1)
print(avg_table_data1)
# simulated data1
# ACC_all   ACC_u ACC_v ACC_w   TPR_all TPR_u     TPR_v TPR_w   TNR_all     TNR_u     TNR_v TNR_w  time
# PLS           0.3260870 0.25000 0.300 0.500 1.0000000 1.000 1.0000000 1.000 0.0000000 0.0000000 0.0000000 0.000 0.003
# PMD           0.7995652 0.91625 0.856 0.500 0.8800000 0.815 0.8233333 1.000 0.7606452 0.9500000 0.8700000 0.000 0.021
# L0sPLS        0.8695652 0.98000 0.966 0.500 0.9666667 0.960 0.9433333 1.000 0.8225806 0.9866667 0.9757143 0.000 0.018
# L2wsPLS       0.7095652 0.72750 0.692 0.716 0.5546667 0.455 0.4866667 0.716 0.7845161 0.8183333 0.7800000 0.716 0.025
# L_infty_wsPLS 0.9791304 0.98000 0.972 0.992 0.9680000 0.960 0.9533333 0.992 0.9845161 0.9866667 0.9800000 0.992 0.051

avg_table_data2 = simulation(flag=2)
print(avg_table_data2)
# simulated data2
# ACC_all   ACC_u  ACC_v ACC_w   TPR_all  TPR_u     TPR_v TPR_w   TNR_all  TNR_u     TNR_v TNR_w  time
# PLS           0.2894737 0.25000 0.3000 0.500 1.0000000 1.0000 1.0000000 1.000 0.0000000 0.0000 0.0000000 0.000 0.065
# PMD           0.8681579 0.79550 0.9631 0.500 0.8161818 0.6245 0.9133333 1.000 0.8893333 0.8525 0.9844286 0.000 0.087
# L0sPLS        0.9256842 0.88975 0.9970 0.500 0.9170909 0.7795 0.9950000 1.000 0.9291852 0.9265 0.9978571 0.000 0.109
# L2wsPLS       0.6169474 0.62725 0.6050 0.654 0.3383636 0.2545 0.3416667 0.654 0.7304444 0.7515 0.7178571 0.654 0.026
# L_infty_wsPLS 0.9527368 0.89125 0.9972 1.000 0.9183636 0.7825 0.9953333 1.000 0.9667407 0.9275 0.9980000 1.000 0.114


avg_table_data3 = simulation(flag=3)
print(avg_table_data3)

# TPR_all   TPR_u     TPR_v  TPR_w   TNR_all     TNR_u     TNR_v  TNR_w   ACC_all    ACC_u   ACC_v  ACC_w   time
# PLS           1.0000000 1.00000 1.0000000 1.0000 0.0000000 0.0000000 0.0000000 0.0000 0.2837838 0.250000 0.30000 0.5000 11.396
# PMD           0.8921714 0.82065 0.9308667 1.0000 0.9596000 0.9524500 1.0000000 0.0000 0.9404649 0.919500 0.97926 0.5000 11.504
# L0sPLS        0.9823429 0.95365 1.0000000 1.0000 0.9741358 0.9845500 1.0000000 0.0000 0.9764649 0.976825 1.00000 0.5000 25.141
# L2wsPLS       0.8600000 0.81300 0.8804333 0.9908 0.9445283 0.9376667 0.9487571 0.9908 0.9205405 0.906500 0.92826 0.9908  0.361
# L_infty_wsPLS 0.9823048 0.95355 1.0000000 1.0000 0.9929887 0.9845167 1.0000000 1.0000 0.9899568 0.976775 1.00000 1.0000  5.774

save.image("F1_results_summary.RData")