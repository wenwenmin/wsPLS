library("ggplot2")
library("cowplot")
library("circlize")
library("gridExtra")

################################################################################
#Set Working Directory
setwd("C:/Users/minwenwen/Desktop/github_wsPLS/Simulation2_mwsPLS")
print(getwd())

################################################################################
# Generate simulated data
n = 50
p = 100
q = 200
t = 300

# Number of the points which are truly nonzero
k1 = 20
k2 = 40
k3 = 60
kn = 25

set.seed(1)
u1 = c(rnorm(20),  rep(0,p-20))
u2 = c(rnorm(40),  rep(0,q-40))
u3 = c(rnorm(60),  rep(0,t-60))
w = c(rep( 1,25), rep( 0,n-25))

set.seed(1)
X1 = w%*%t(u1) + matrix(rnorm(n*p),ncol=p)
X2 = w%*%t(u2) + matrix(rnorm(n*q),ncol=q)
X3 = w%*%t(u3) + matrix(rnorm(n*t),ncol=t) 

################################################################################
source('mwsPLS_scheme1.R')
out1 = mwsPLS_scheme1(X1, X2, X3, k1, k2, k3, kn, Lc=10, niter=100, err=10^(-3), nstart=5)

source('mwsPLS_scheme2.R')
out2 = mwsPLS_scheme2(X1, X2, X3, k1, k2, k3, kn, Lc=100, niter=100, err=10^(-3), nstart=5)


################################################################################
# Plot figures
get_u1_fig = function(u1Dat, Meth="mwsPLS-scheme1 u"){
  colors =  c("black","red")
  f1 = ggplot(data=u1Dat, aes(x=x, y = value, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = Meth,x=NULL,y=NULL)+guides(color=FALSE) + theme_bw()
  
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  
  f1 = f1 + my_theme + scale_x_continuous(breaks=as.numeric(c(0, 20, 60,100)))
  y_breaks = c(round(min(u1Dat$value),1), 0, round(max(u1Dat$value),1))
  f1 = f1 + my_theme + scale_y_continuous(breaks=y_breaks)
}

get_u2_fig = function(u2Dat,Meth="mwsPLS-scheme1 u"){
  colors =  c("black","red")
  f1 = ggplot(data=u2Dat, aes(x=x, y = value, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = Meth,x=NULL,y=NULL)+guides(color=FALSE) + theme_bw()
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  f1 = f1 + my_theme + scale_x_continuous(breaks=as.numeric(c(0, 40,110,200)))
  y_breaks = c(round(min(u2Dat$value)+0.1,1), 0, round(max(u2Dat$value),1))
  f1 = f1 + my_theme + scale_y_continuous(breaks=y_breaks)
}

get_u3_fig = function(u3Dat,Meth="mwsPLS-scheme1 u"){
  colors =  c("black","red")
  f1 = ggplot(data=u3Dat, aes(x=x, y = value, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = Meth,x=NULL,y=NULL)+guides(color=FALSE) + theme_bw()
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  f1 = f1 + my_theme + scale_x_continuous(breaks=as.numeric(c(0, 60,180,300)))
  
  y_breaks = c(round(min(u3Dat$value),1), 0, round(max(u3Dat$value),1))
  f1 = f1 + my_theme + scale_y_continuous(breaks=y_breaks)
}

get_w_fig = function(wDat,Meth="mwsPLS-scheme1 w"){
  colors =  c("black","red")
  f1 = ggplot(data=wDat, aes(x=x, y = value, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = Meth,x=NULL,y=NULL)+guides(color=FALSE) + theme_bw()
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  f1 = f1 + my_theme + scale_x_continuous(breaks=as.numeric(c(0, 25,50)))
  f1 = f1 + my_theme + scale_y_continuous(breaks=as.numeric(c( 0.0, 0.5, 1.0)))
}


get_fig_for_a_method = function(out1,u1,u2,u3,w,Meth0="mwsPLS-scheme1"){
  
  u1Dat = data.frame(x = 1:length(u1), value = out1$u1, group = factor(ifelse(u1!=0,1,0)))
  u2Dat = data.frame(x = 1:length(u2), value = out1$u2, group = factor(ifelse(u2!=0,1,0)))
  u3Dat = data.frame(x = 1:length(u3), value = out1$u3, group = factor(ifelse(u3!=0,1,0)))
  wDat  = data.frame(x = 1:length(w),  value = out1$w, group = factor(ifelse(w!=0,1,0)))
  
  figs = list()
  figs$u1 = get_u1_fig(u1Dat, Meth=paste(Meth0,"u1"))
  figs$u2 = get_u2_fig(u2Dat, Meth=paste(Meth0,"u2"))
  figs$u3 = get_u3_fig(u3Dat, Meth=paste(Meth0,"u3"))
  figs$w  = get_w_fig( wDat,  Meth=paste(Meth0,"w"))
  
  figs = plot_grid(figs$u1,figs$u2, figs$u3,figs$w, ncol=4, nrow = 1, align="hv", labels =c(""))
}

out = list(u1=u1,u2=u2,u3=u3,w=w)
m0_fig = get_fig_for_a_method(out, u1,u2,u3,w,Meth0="True")
m1_fig = get_fig_for_a_method(out1,u1,u2,u3,w,Meth0="mwsPLS-scheme1")
m2_fig = get_fig_for_a_method(out2,u1,u2,u3,w,Meth0="mwsPLS-scheme2")

figs = plot_grid(m0_fig,m1_fig,m2_fig,ncol=1, nrow = 3, align="hv", labels =c("A","B","C"))


################################################################################
ggsave(file="Figs_mwsPLS.jpg",width=11.8, height=7)



