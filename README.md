# wsPLS
Sparse Partial Least Squares (sPLS) is a common dimensionality reduction technique for data fusion,
which projects data samples from two views by seeking linear combinations with a small number of variables with the maximum variance.
However, sPLS extracts the combinations between two data sets with all data samples so that it cannot detect latent subsets of samples.
To extend the application of sPLS by identifying a specific subset of samples and remove outliers,
we propose an $\ell_\infty/\ell_0$-norm constrained weighted sparse PLS method for joint sample and feature selection,
where the $\ell_\infty/\ell_0$-norm constrains are used to select a subset of samples.
We prove that the $\ell_\infty/\ell_0$-norm constrains have the Kurdyka-\L{ojasiewicz}~property so that a globally convergent algorithm is developed to solve it.
Moreover, multi-view data with a same set of samples can be available in various real problems.
To this end, we extend the $\ell_\infty/\ell_0$-wsPLS model and propose two multi-view wsPLS models for multi-view data fusion.
We develop an efficient iterative algorithm for each multi-view wsPLS model and show its convergence property.
As well as numerical and biomedical data experiments demonstrate the efficiency of the proposed methods.

<p align="center"> 
<img src="https://github.com/wenwenmin/wsPLS/blob/main/Figure_1_overview_wsPLS.png">
</p>



### R code for $\ell_\infty/\ell_0$-wsPLS
The first example explains how to use the $\ell_\infty/\ell_0$-wsPLS algorithm (please see https://github.com/wenwenmin/wsPLS/blob/main/Simulation1_wsPLS/Example1_simulation_for_wsPLS.R)
Before running the script, please first set the path for "Example1_simulation_for_wsPLS.R",
and then run the following R command in the Console. 

More descriptions about wPLS can be seen in the file "README_R_markdown_file.Rmd"(https://github.com/wenwenmin/wsPLS/blob/main/README_R_markdown_file.Rmd).


``` r
> source('Example1_simulation_for_wsPLS.R') 
```

### R code for mwsPLS
The second example explains how to use the mwsPLS algorithm (See Simulation2_mwsPLS/Example2_simulation_for_mwsPLS.R). 
Before running the script, please first set the path for "Example2_simulation_for_mwsPLS.R",
and then run the following R command in the Console. 

``` r
> source('Example2_simulation_for_mwsPLS.R') 
```


### R function of $\ell_\infty/\ell_0$-wsPLS
```{r cars}
wsPLS = function(X, Y, ku, kv, kw, Lc=0.1, niter=100, err=10^(-5), nstart=5, seed0=1){
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
```
