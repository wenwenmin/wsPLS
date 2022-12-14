mwsPLS_scheme2 = function(X1, X2, X3, k1, k2, k3, kn, Lc=0.1, niter=100, err=10^(-3), nstart=5){
  # ----------------------------------------------------------------------------
  # Input data
  # X \in R^{n \times p1} (n:samples, p1:variables)
  # Y \in R^{n \times p2} (n:samples, p2:variables)
  # Z \in R^{n \times p3} (n:samples, p3:variables)
  # Lc: Lipschitz constant
  # Sparse level: {ku, kv, kw, kn}
  # ----------------------------------------------------------------------------
  # The optimization problems solved are as follows
  # \begin{equation}\label{equ-15}
  # \begin{aligned}
  # & \underset{\bm{u}_i, \bm{w}}{\text{minimize}} &&- \bm{w}^{\rm{T}}\big[\bigodot\limits_{i=1}^M(\bm{X}_i\bm{u}_i)\big] \\
  # & \text{subject to} &&\|\bm{u}_i\| = 1,\|\bm{u}_i\|_0\leq k_i~\forall  i\\
  # &                   &&\|\bm{w}\|_0\leq k_w, \|\bm{w}\|_\infty \leq 1, \bm{w}\geq 0,\\
  # \end{aligned}
  # \end{equation}
  # ----------------------------------------------------------------------------
  
  # Initialize optimal objective function
  opt_obj = 0
  
  # Initialize the objective function value of each initial point
  obj_nstart = rep(0,nstart); names(obj_nstart) = paste("init",1:nstart,sep = "_")
  
  for(init in 1:nstart){
    # Initialize u v and w
    set.seed(init)
    u1 = u1_0 = matrix(rnorm(ncol(X1)),ncol=1); u1 = u1/norm(u1,'E')
    u2 = u2_0 = matrix(rnorm(ncol(X2)),ncol=1); u2 = u2/norm(u2,'E')
    u3 = u3_0 = matrix(rnorm(ncol(X3)),ncol=1); u3 = u3/norm(u3,'E')
    w  = w0   = matrix(rep(1,nrow(X1)),ncol=1)  #Each element is one
    
    objs = c()
    # Iterative algorithm
    for (i in 1:niter){
      # update u1
      z1 = crossprod(X1, w*(X2%*%u2)*(X3%*%u3))
      u1 = PALM_sproject(z1, u1_0, Lc, k1)
      
      # update u2
      z2 = crossprod(X2, w*(X1%*%u1)*(X3%*%u3))
      u2 = PALM_sproject(z2, u2_0, Lc, k2)
      
      # update u3
      z3 = crossprod(X3, w*(X1%*%u1)*(X2%*%u2))
      u3 = PALM_sproject(z3, u3_0, Lc, k3)
      
      # update w
      z4 = (X1%*%u1)*(X2%*%u2)*(X3%*%u3)
      w  = PALM_update_w_l0_equal_k(z4, w0, Lc, kn)
      
      # Record objective function value
      obj = sum((X1%*%u1)*(X2%*%u2)*(X3%*%u3)*w)
      objs = c(objs,obj)
      
      # Algorithm termination condition
      if(i>=4&&abs(objs[i]-objs[i-1])<= err){
        break}else{
          u1_0=u1;u2_0=u2;u3_0=u3;w_0=w}
    }
    
    # --------------------------------------------------------------------------
    # Record objective every initialization
    obj_nstart[init] = obj
    # Keep optimal solution
    if(obj>opt_obj){
      opt_obj = obj
      opt_out = list(u1      = u1, 
                     u2      = u2,
                     u3      = u3,
                     w       = w, 
                     objs    = objs, 
                     opt_obj = opt_obj)
      
    }
  }
  
  # ----------------------------------------------------------------------------
  # Output
  out = list(u1         = opt_out$u1, 
             u2         = opt_out$u2, 
             u3         = opt_out$u3, 
             w          = opt_out$w, 
             objs       = opt_out$objs, 
             opt_obj    = opt_out$opt_obj,
             obj_nstart = obj_nstart)
  return(out)
}

# ----------------------------------------------------------------------------
# Sparse project function
# s.t. ||w||_0 <= k_w,  0<=w_i<=1
PALM_sproject = function(z, u0, Lc, k){
  #Lc is the Lipschitz constant
  z = u0 + z/Lc
  if((length(z)<k)|(sum(z^2)==0)) return(z)  
  u = abs(z);
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u) 
}
# ----------------------------------------------------------------------------
# Learning w
# s.t. ||w||_0 <= k_w,  w_i \in {0,1}^n
PALM_update_w_l0_equal_k = function(z, w0, Lc, k){
  #Lc is the Lipschitz constant
  z = w0 + z/Lc
  if((length(z)<=k)|(sum(z^2)==0)) {
    w = matrix(rep(1,length(z)),ncol=1)
    return(w) }

  w = matrix(rep(1,length(z)),ncol=1) # ones
  w[-order(z,decreasing=T)[1:k]] = 0

  return(w)
}
# ----------------------------------------------------------------------------
# PALM_update_w_l0_equal_k = function(z, w0, Lc, k){
#   #Lc is the Lipschitz constant
#   z = w0 + z/Lc
#   if((length(z)<=k)|(sum(z^2)==0)) {
#     w = matrix(rep(1,length(z)),ncol=1)
#     return(w) }
#   z[-order(z,decreasing=T)[1:k]] = 0
#   z = ifelse(z>0,1,-1)
#   return(w) 
# }
# ----------------------------------------------------------------------------