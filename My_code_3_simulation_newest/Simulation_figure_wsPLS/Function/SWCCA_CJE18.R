SWCCA_CJE18 = function(X, Y, ku, kv, kw, seed0=100, niter=100, err=0.0001){
  # ----------------------------------------------------------------------------
  # Input
  # X \in R^{n \times p} (n:samples, p:variables)
  # Y \in R^{n \times q} (n:samples, q:variables)
  # Ref CJE 18 Sparse Weighted Canonical Correlation Analysis
  # ----------------------------------------------------------------------------
  set.seed(seed0)
  u0 = matrix(rnorm(ncol(X)),ncol=1);u0 = u0/norm(u0,'E')
  v0 = matrix(rnorm(ncol(Y)),ncol=1);v0 = v0/norm(v0,'E')
  w0 = matrix(rnorm(nrow(X)),ncol=1);w0 = w0/norm(w0,'E')
  
  u = u0; v = v0; w = w0
  # Iterative algorithm
  for (i in 1:niter){
    z.u = crossprod(X, w*(Y%*%v))
    u = SWCCA_CJE18_project(z.u, ku)
    
    z.v = crossprod(Y, w*(X%*%u))
    v = SWCCA_CJE18_project(z.v, kv) 
    
    z.w = (X%*%u)*(Y%*%v)
    w = SWCCA_CJE18_project(z.w, kw)
    
    # Algorithm termination condition
    if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)&(norm(w - w0,'E')<= err)){break}
    else {
      u0=u; v0=v; w0=w}
  }
  obj = sum((X%*%u)*(Y%*%v)*w)
  return (list(u=u, v=v, w=w, obj=obj))
}

# Sparse project function
SWCCA_CJE18_project = function(z, k){
  if((length(z)<k)|(sum(z^2)==0)) return(z)  
  u = abs(z);
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u) 
}
