library("PMA")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("latex2exp")

# source('L0CCA.R')
# source('L0SBCCA.R')

n = 50
p = 80
q = 100

u = c(rep(1,10),  rep(-1,10), rep(0,p-20))
v = c(rep(-1,15), rep(1,15),  rep(0,q-30))
w = c(rep(1,25),  rep(0,n-25))

# set.seed(1)
# X = w%*%t(u) + matrix(rnorm(n*p),ncol=p)
# Y = w%*%t(v) + matrix(rnorm(n*q),ncol=q)

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
library("PMA") # sCCA with L1-norm constrain
out0 = CCA(X, Y, typex="standard", typez="standard", K=1, penaltyx = 1, penaltyz = 1)

library("PMA") # sCCA with L1-norm constrain
out1 = CCA(X, Y, typex="standard", typez="standard", K=1, penaltyx = 0.4, penaltyz = 0.467)

source('Function/sCCA.R')
out2 = L0_sCCA(X, Y, ku, kv, niter=100, err=10^(-5), nstart=5)

source('Function/SWCCA_CJE18.R')
out3 = SWCCA_CJE18(X, Y, ku, kv, kw, seed0=1)

source('Function/PALM_wsPLS.R')
out4 = PALM_wsPLS(X, Y, ku, kv, kw, Lc=0.1, nstart=5, seed0=2)
# ------------------------------------------------------------------------------
################################################################################

source('F1_plot_fun.R')

# Plot true pattern
u_true_df = data.frame(x = 1:length(u), u = u, group = factor(ifelse(u!=0,1,0)))
v_true_df = data.frame(x = 1:length(v), v = v, group = factor(ifelse(v!=0,1,0)))
w_true_df = data.frame(x = 1:length(w), w = w, group = factor(ifelse(w!=0,1,0)))
fig1 = get_figs(u_true_df,v_true_df,w_true_df,Meth="True")
#ggsave(file="Res_figures/Fig1_true.jpg",width=8.3, height=1.6)

# PLS
u0 = data.frame(x = 1:length(u), u = out0$u, group = factor(ifelse(u!=0,1,0)))
v0 = data.frame(x = 1:length(v), v = out0$v, group = factor(ifelse(v!=0,1,0)))
fig2 = get_figs_uv(u0,v0, Meth="PLS")
#ggsave(file="Res_figures/Fig2_PLS.jpg",width=8.3, height=1.6)

# PMD
u1 = data.frame(x = 1:length(u), u = out1$u, group = factor(ifelse(u!=0,1,0)))
v1 = data.frame(x = 1:length(v), v = out1$v, group = factor(ifelse(v!=0,1,0)))
fig3 = get_figs_uv(u1,v1, Meth="PMD-sPLS")
#ggsave(file="Res_figures/Fig3_PMD-sPLS.jpg",width=8.3, height=1.6)

# L0-sPLS
u2 = data.frame(x = 1:length(u), u = out2$u, group = factor(ifelse(u!=0,1,0)))
v2 = data.frame(x = 1:length(v), v = out2$v, group = factor(ifelse(v!=0,1,0)))
fig4 = get_figs_uv(u2,v2,Meth="l0-sPLS")
# fig4 + labs(title ="XXX",x="NULL",y=NULL)

#ggsave(file="Res_figures/Fig4_L0-sPLS.jpg",width=8.3, height=1.6)

# L2-wsPLS
u3 = data.frame(x = 1:length(u), u = out3$u, group = factor(ifelse(u!=0,1,0)))
v3 = data.frame(x = 1:length(v), v = out3$v, group = factor(ifelse(v!=0,1,0)))
w3 = data.frame(x = 1:length(w), w = out3$w, group = factor(ifelse(w!=0,1,0)))
fig5 = get_figs(u3,v3,w3,Meth="l2-wsPLS")
#ggsave(file="Res_figures/Fig5_L2-wsPLS.jpg",width=8.3, height=1.6)

# l_infty_wsPLS
u4 = data.frame(x = 1:length(u), u = out4$u, group = factor(ifelse(u!=0,1,0)))
v4 = data.frame(x = 1:length(v), v = out4$v, group = factor(ifelse(v!=0,1,0)))
w4 = data.frame(x = 1:length(w), w = out4$w, group = factor(ifelse(w!=0,1,0)))
fig6 = get_figs(u4,v4,w4,Meth="l_infty_wsPLS")
#ggsave(file="Res_figures/Fig6_l_infty_wsPLS.jpg",width=8.3, height=1.6)


figs = plot_grid(fig1,fig2, fig3,fig4,fig5,fig6, 
                ncol=1, nrow = 6, align="hv", labels =c("A","B","C","D","E","F"))
ggsave(file="Figure2_diff_methods+.jpg",width=8.3, height=1.6*6.5)
ggsave(file="Figure2_diff_methods+.pdf",width=8.3, height=1.6*6.5)

# ------------------------------------------------------------------------------
save.image("Result_output.RData")
