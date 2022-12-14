library("PMA")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("latex2exp")

################################################################################
#Set Working Directory
setwd("C:/Users/minwenwen/Desktop/github_wsPLS/Simulation1_wsPLS")
print(getwd())


################################################################################
# Generate simulated data
n = 50
p = 80
q = 100

u = c(rep(1,10),  rep(-1,10), rep(0,p-20))
v = c(rep(-1,15), rep(1,15),  rep(0,q-30))
w = c(rep(1,25),  rep(0,n-25))

SNR = 0.1
gamma1 = sqrt(norm(u%*%t(v),"F")^2/(SNR*n*p))
gamma2 = sqrt(norm(w%*%t(v),"F")^2/(SNR*n*q))

set.seed(1)
X = w%*%t(u) + gamma1*matrix(rnorm(n*p),ncol=p)
Y = w%*%t(v) + gamma2*matrix(rnorm(n*q),ncol=q)

ku = 20
kv = 30
kw = 25
################################################################################
# ------------------------------------------------------------------------------
# PLS
library("PMA") 
out0 = CCA(X, Y, typex="standard", typez="standard", K=1, penaltyx = 1, penaltyz = 1)

# sPLS with L1-norm constrain
library("PMA") 
out1 = CCA(X, Y, typex="standard", typez="standard", K=1, penaltyx = 0.4, penaltyz = 0.467)

# L0_sPLS
source('L0_sPLS.R')
out2 = L0_sPLS(X, Y, ku, kv, niter=100, err=10^(-5), nstart=5)

# L2-wsPLS is solved by SWCCA_CJE18
source('SWCCA_CJE18.R')
out3 = SWCCA_CJE18(X, Y, ku, kv, kw, seed0=1)

# l_infty_wsPLS
source('wsPLS.R')
out4 = wsPLS(X, Y, ku, kv, kw, Lc=0.1, nstart=5, seed0=2)
# ------------------------------------------------------------------------------
################################################################################

source('Functions_plot_figure.R')

# Plot true pattern
u_true_df = data.frame(x = 1:length(u), u = u, group = factor(ifelse(u!=0,1,0)))
v_true_df = data.frame(x = 1:length(v), v = v, group = factor(ifelse(v!=0,1,0)))
w_true_df = data.frame(x = 1:length(w), w = w, group = factor(ifelse(w!=0,1,0)))
fig1 = get_figs(u_true_df,v_true_df,w_true_df,Meth="True")

# PLS
u0 = data.frame(x = 1:length(u), u = out0$u, group = factor(ifelse(u!=0,1,0)))
v0 = data.frame(x = 1:length(v), v = out0$v, group = factor(ifelse(v!=0,1,0)))
fig2 = get_figs_uv(u0,v0, Meth="PLS")

# PMD
u1 = data.frame(x = 1:length(u), u = out1$u, group = factor(ifelse(u!=0,1,0)))
v1 = data.frame(x = 1:length(v), v = out1$v, group = factor(ifelse(v!=0,1,0)))
fig3 = get_figs_uv(u1,v1, Meth="PMD-sPLS")

# L0-sPLS
u2 = data.frame(x = 1:length(u), u = out2$u, group = factor(ifelse(u!=0,1,0)))
v2 = data.frame(x = 1:length(v), v = out2$v, group = factor(ifelse(v!=0,1,0)))
fig4 = get_figs_uv(u2,v2,Meth="l0-sPLS")
# fig4 + labs(title ="XXX",x="NULL",y=NULL)


# L2-wsPLS
u3 = data.frame(x = 1:length(u), u = out3$u, group = factor(ifelse(u!=0,1,0)))
v3 = data.frame(x = 1:length(v), v = out3$v, group = factor(ifelse(v!=0,1,0)))
w3 = data.frame(x = 1:length(w), w = out3$w, group = factor(ifelse(w!=0,1,0)))
fig5 = get_figs(u3,v3,w3,Meth="l2-wsPLS")

# l_infty_wsPLS
u4 = data.frame(x = 1:length(u), u = out4$u, group = factor(ifelse(u!=0,1,0)))
v4 = data.frame(x = 1:length(v), v = out4$v, group = factor(ifelse(v!=0,1,0)))
w4 = data.frame(x = 1:length(w), w = out4$w, group = factor(ifelse(w!=0,1,0)))
fig6 = get_figs(u4,v4,w4,Meth="l_infty_wsPLS")


figs = plot_grid(fig1,fig2, fig3,fig4,fig5,fig6, 
                ncol=1, nrow = 6, align="hv", labels =c("A","B","C","D","E","F"))
ggsave(file="Figs_wsPLS.pdf",width=8.3, height=1.6*6.5)
