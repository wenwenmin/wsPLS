L0_sCCA = function(X, Y, ku, kv, niter=100, err=10^(-5), nstart=5){
  # ----------------------------------------------------------------------------
  # Input
  # X \in R^{n \times p} (n:samples, p:variables)
  # Y \in R^{n \times q} (n:samples, q:variables)
  # ----------------------------------------------------------------------------
  
  XTY = crossprod(X,Y)
  # Initialize optimal objective function
  opt_obj = 0
  
  # Initialize the objective function value of each initial point
  obj_nstart = rep(0,nstart); names(obj_nstart) = paste("init",1:nstart,sep = "_")
  
  for(init in 1:nstart){
    
    # Initialize u v and w
    set.seed(init)
    u0 = matrix(rnorm(ncol(X)),ncol=1);u0 = u0/norm(u0,'E')
    v0 = matrix(rnorm(ncol(Y)),ncol=1);v0 = v0/norm(v0,'E')
    u = u0; v = v0
    
    objs = c()
    # Iterative algorithm to solve u and v
    for (i in 1:niter){
      # update u
      z.u = XTY%*%v
      u = L0_Ball_Project(z.u, ku) 
      # update v
      z.v = crossprod(XTY,u) 
      v = L0_Ball_Project(z.v, kv) 
      # Record objective 
      obj = t(u)%*%XTY%*%v
      objs = c(objs,obj)
      # Algorithm termination condition
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)){
        break}else {
          u0 = u;v0 = v}
    }
   
    # Record objective every initialization
    obj_nstart[init] = obj
    # Keep optimal solution
    if(obj>opt_obj){
      opt_out = list(u=u, v=v, objs=objs, opt_obj=obj)
      opt_obj = obj
    }
  }

  # Output
  out = list(u          = opt_out$u, 
             v          = opt_out$v, 
             objs       = opt_out$objs, 
             opt_obj    = opt_out$opt_obj,
             obj_nstart = obj_nstart)
  return(out)
}

# Sparse project function
L0_Ball_Project = function(z, k){
  if((length(z)<=k)|(sum(z^2)==0)) return(z)  
  u = abs(z);
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u) 
}