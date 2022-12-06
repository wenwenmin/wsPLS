################################################################################
################################################################################
Fun_TPR_TNR_ACC_1 = function(uu0, uu1){
  # # example 1
  # # True Label
  # uu0 = c(1,   1,  1, 1,   1,   0, 0, 0, 0, 0)
  # #Predicted Label
  # uu1 = c(0.9, 08, 0, 0.7, 0.9, 0, 0, 0, 0, 0.1)
  # 
  # # example 2
  # uu0 = data$u
  # uu1 = out2$u
  
  # The number of positive samples  
  P  = length(which(uu0!=0))
  
  # The number of (uu0!=0 and uu1 !=0)
  TP = length(which(uu1[uu0!=0]!=0))
  
  N  = length(which(uu0==0))
  
  TN = length(which(uu1[uu0==0]==0))
  
  TPR = TP/P
  TNR = TN/N
  ACC = (TP+TN)/(P+N) 
  
  return(list(TPR=TPR,TNR=TNR,ACC=ACC))
}

Fun_TPR_TNR_ACC_2 = function(data, out1){
  
  TPR_u   = Fun_TPR_TNR_ACC_1(data$u, out1$u)$TPR
  TPR_v   = Fun_TPR_TNR_ACC_1(data$v, out1$v)$TPR
  TPR_w   = Fun_TPR_TNR_ACC_1(data$w, out1$w)$TPR
  TPR_all = Fun_TPR_TNR_ACC_1(c(data$u, data$v, data$w), c(out1$u, out1$v, out1$w))$TPR
  TPR     = c(TPR_all=TPR_all, TPR_u=TPR_u, TPR_v=TPR_v, TPR_w=TPR_w)
  
  TNR_u   = Fun_TPR_TNR_ACC_1(data$u, out1$u)$TNR
  TNR_v   = Fun_TPR_TNR_ACC_1(data$v, out1$v)$TNR
  TNR_w   = Fun_TPR_TNR_ACC_1(data$w, out1$w)$TNR
  TNR_all = Fun_TPR_TNR_ACC_1(c(data$u, data$v, data$w), c(out1$u, out1$v, out1$w))$TNR
  TNR     = c(TNR_all=TNR_all, TNR_u=TNR_u, TNR_v=TNR_v, TNR_w=TNR_w)
  
  ACC_u   = Fun_TPR_TNR_ACC_1(data$u, out1$u)$ACC
  ACC_v   = Fun_TPR_TNR_ACC_1(data$v, out1$v)$ACC
  ACC_w   = Fun_TPR_TNR_ACC_1(data$w, out1$w)$ACC
  ACC_all = Fun_TPR_TNR_ACC_1(c(data$u, data$v, data$w), c(out1$u, out1$v, out1$w))$ACC
  ACC     = c(ACC_all=ACC_all, ACC_u=ACC_u, ACC_v=ACC_v, ACC_w=ACC_w)
  
  return(list(TPR=TPR,TNR=TNR,ACC=ACC))
}

Fun_TPR_TNR_ACC_table = function(data, out1, out2, out3, out4, out5){
  
  TPR_table = data.frame(matrix(0,5,4))
  TPR_table[1,] = Fun_TPR_TNR_ACC_2(data, out1)$TPR # 1-PLS
  TPR_table[2,] = Fun_TPR_TNR_ACC_2(data, out2)$TPR # 2-PMD
  TPR_table[3,] = Fun_TPR_TNR_ACC_2(data, out3)$TPR # 3-L0-sPLS
  TPR_table[4,] = Fun_TPR_TNR_ACC_2(data, out4)$TPR # 4-L2-wsPLS
  TPR_table[5,] = Fun_TPR_TNR_ACC_2(data, out5)$TPR # l_infty_wsPLS
  
  TNR_table = data.frame(matrix(0,5,4))
  TNR_table[1,] = Fun_TPR_TNR_ACC_2(data, out1)$TNR # 1-PLS
  TNR_table[2,] = Fun_TPR_TNR_ACC_2(data, out2)$TNR # 2-PMD
  TNR_table[3,] = Fun_TPR_TNR_ACC_2(data, out3)$TNR # 3-L0-sPLS
  TNR_table[4,] = Fun_TPR_TNR_ACC_2(data, out4)$TNR # 4-L2-wsPLS
  TNR_table[5,] = Fun_TPR_TNR_ACC_2(data, out5)$TNR # l_infty_wsPLS
  
  ACC_table = data.frame(matrix(0,5,4))
  ACC_table[1,] = Fun_TPR_TNR_ACC_2(data, out1)$ACC # 1-PLS
  ACC_table[2,] = Fun_TPR_TNR_ACC_2(data, out2)$ACC # 2-PMD
  ACC_table[3,] = Fun_TPR_TNR_ACC_2(data, out3)$ACC # 3-L0-sPLS
  ACC_table[4,] = Fun_TPR_TNR_ACC_2(data, out4)$ACC # 4-L2-wsPLS
  ACC_table[5,] = Fun_TPR_TNR_ACC_2(data, out5)$ACC # l_infty_wsPLS
  
  Res_table = cbind(cbind(ACC_table,TPR_table),TNR_table)
  
  
  row.names(Res_table) = c("PLS", "PMD", "L0sPLS", "L2wsPLS", "L_infty_wsPLS")
  
  colnames(Res_table)  = c("ACC_all", "ACC_u", "ACC_v", "ACC_w",
                           "TPR_all", "TPR_u", "TPR_v", "TPR_w",
                           "TNR_all", "TNR_u", "TNR_v", "TNR_w")
  return(Res_table)
}

