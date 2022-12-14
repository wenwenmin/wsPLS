# wsPLS
Sparse Partial Least Squares (sPLS) is a common dimensionality reduction technique for data fusion,
which projects data samples from two views by seeking linear combinations with a small number of variables with the maximum variance.
However, sPLS extracts the combinations between two data sets with all data samples so that it cannot detect latent subsets of samples.
To extend the application of sPLS by identifying a specific subset of samples and remove outliers,
we propose an $\ell_\infty/\ell_0$-norm constrained weighted sparse PLS ($\ell_\infty/\ell_0$-wsPLS) method for joint sample and feature selection,
where the $\ell_\infty/\ell_0$-norm constrains are used to select a subset of samples.
We prove that the $\ell_\infty/\ell_0$-norm constrains have the Kurdyka-\L{ojasiewicz}~property so that a globally convergent algorithm is developed to solve it.
Moreover, multi-view data with a same set of samples can be available in various real problems.
To this end, we extend the $\ell_\infty/\ell_0$-wsPLS model and propose two multi-view wsPLS models for multi-view data fusion.
We develop an efficient iterative algorithm for each multi-view wsPLS model and show its convergence property.
As well as numerical and biomedical data experiments demonstrate the efficiency of the proposed methods.

More descriptions about wPLS can be seen in the file "README_R_markdown_file.Rmd"(https://github.com/wenwenmin/wsPLS/blob/main/README_R_markdown_file.Rmd).

<p align="center"> 
<img src="https://github.com/wenwenmin/wsPLS/blob/main/Fig0_wsPLS.png">
</p>

### R function of wsPLS
```{r cars}
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
```

### A simple exaple to explain wsPLS
```{r pressure, echo=T}
library("ComplexHeatmap")
library("circlize")

# ------------------------------------------------------------------------------
# Generate simulation data
n = 16
p = 10*2
q = 12*2

u1 = c(rep( 1, 4),  rep(0,p-4))
v1 = c(rep(-1, 6),  rep(0,q-6))
w1 = c(rep( 1, 8),  rep(0,n-8))

u2 = c(rep( 0, 4), rep(-1, 4), rep(0,p-8))
v2 = c(rep( 0, 6), rep( 1, 6), rep(0,q-12))
w2 = c(rep( 0, 8), rep( 1, 8))

ku = 4
kv = 6
kw = 8

set.seed(1)
X0 = w1%*%t(u1) +  w2%*%t(u2) +  0.2*matrix(rnorm(n*p),ncol=p)
Y0 = w1%*%t(v1) +  w2%*%t(v2) +  0.2*matrix(rnorm(n*q),ncol=q)

set.seed(1)
u_order = sample(1:p, size=p, replace = F)
v_order = sample(1:q, size=q, replace = F)
w_order = sample(1:n, size=n, replace = F)

X = X0[w_order,u_order]
Y = Y0[w_order,v_order] 

# ----------------------------------------------------------------------------
# We use wsPLS algorithm on the above datasets
out = PALM_wsPLS(X, Y, ku, kv, kw, Lc=0.1, nstart=5, seed0=1)

u_order2 = order(-abs(out$u))
v_order2 = order(-abs(out$v))
w_order2 = order(-abs(out$w))

X_order = X[w_order2,u_order2]
Y_order = Y[w_order2,v_order2] 


# ------------------------------------------------------------------------------
# plot two submatrices
ht3 = Heatmap(X_order, 
              column_title = "X_order",
              cluster_rows = F,
              cluster_columns = F, 
              column_names_side = "top",
              row_names_side = "left",
              show_heatmap_legend = FALSE, name="foo1")
ht4 = Heatmap(Y_order, 
              column_title = "Y_order",
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_side = "top",
              row_names_side = "left",
              show_heatmap_legend = FALSE, name="foo2")
ht_list = ht3 + ht4
ht = draw(ht_list, auto_adjust = T)

decorate_heatmap_body("foo1", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  grid.rect(x = unit(0/p, "npc"), width = unit(ku/p, "npc"),
            y = unit((n-kw)/n, "npc"), height = unit(kw/n, "npc"),
            hjust = 0, vjust = 0,
            gp=gpar(fill = "transparent", col = "red", lwd = 4))
})

decorate_heatmap_body("foo2", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
  grid.rect(x = unit(0/q, "npc"), width = unit(kv/q, "npc"),
            y = unit((n-kw)/n, "npc"), height = unit(kw/n, "npc"),
            hjust = 0, vjust = 0,
            gp=gpar(fill = "transparent", col = "red", lwd = 4))
})
```

