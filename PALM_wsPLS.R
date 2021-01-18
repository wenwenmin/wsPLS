PALM_wsPLS = function(X, Y, ku, kv, kw, Lc=0.1, niter=100, err=10^(-5), nstart=5, seed0=1){
  # ----------------------------------------------------------------------------
  # Input
  # X \in R^{n \times p} (n:samples, p:variables)
  # Y \in R^{n \times q} (n:samples, q:variables)
  # Lc: Lipschitz constant
  # Sparse level: {ku, kv, kw}
  # ----------------------------------------------------------------------------
  
  # Initialize optimal objective function
  opt_obj = 0
  
  # Initialize the objective function value of each initial point
  obj_nstart = rep(0,nstart); names(obj_nstart) = paste("init",1:nstart,sep = "_")
  
  for(init in 1:nstart){
    
    # Initialize u v and w
    set.seed(init*seed0)
    u0 = matrix(rnorm(ncol(X)),ncol=1);u0 = u0/norm(u0,'E')
    v0 = matrix(rnorm(ncol(Y)),ncol=1);v0 = v0/norm(v0,'E')
    # Each element is one
    w0 = matrix(rep(1,nrow(X)),ncol=1)
    u = u0; v = v0; w = w0
    
    objs = cors = c()
    # Iterative algorithm
    for (i in 1:niter){
      # update u
      z.u = crossprod(X, w*(Y%*%v))
      u = PALM_sproject(z.u, u0, Lc, ku)
      # update v
      z.v = crossprod(Y, w*(X%*%u))
      v = PALM_sproject(z.v, v0, Lc, kv) 
      # update w
      z.w = (X%*%u)*(Y%*%v)
      w = PALM_update_w(z.w, w0, Lc, kw)
      
      # Record objective function value
      obj = sum((X%*%u)*(Y%*%v)*w)
      objs = c(objs,obj)
      
      cc = cor((X%*%u)[which(w!=0)], (Y%*%v)[which(w!=0)])
      cors = c(cors,cc)
      
      # Algorithm termination condition
      if ((i>=20)&(norm(u-u0,'E')<= err)&(norm(v-v0,'E')<= err)&(norm(w-w0,'E')<= err)){
        break}else{
          u0=u; v0=v; w0=w}
    }
    
    # Record objective every initialization
    obj_nstart[init] = obj
    # Keep optimal solution
    if(obj>opt_obj){
      opt_out = list(u=u, v=v, w=w, objs=objs, cors=cors, opt_obj=obj)
      opt_obj = obj
    }
  }
  
  # Output
  out = list(u          = opt_out$u, 
             v          = opt_out$v, 
             w          = opt_out$w, 
             objs       = opt_out$objs,
             cors       = opt_out$cors,			 
             opt_obj    = opt_out$opt_obj,
             obj_nstart = obj_nstart)
  return(out)
}
# ----------------------------------------------------------------------------
# Sparse project function
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
PALM_update_w = function(z, w0, Lc, k){
  #Lc is the Lipschitz constant
  z = w0 + z/Lc
  if((length(z)<=k)|(sum(z^2)==0)) {
    w = matrix(rep(1,length(z)),ncol=1)
    return(w) }
  z[z<0] = 0
  z[-order(z,decreasing=T)[1:k]] = 0
  w = sign(z)
  return(w) 
}
# ----------------------------------------------------------------------------
