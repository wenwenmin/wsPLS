# scale_y_continuous(breaks=c(0,50,100,130))
# ------------------------------------------------------------------------------
get_figs = function(u_true_df, v_true_df,w_true_df, Meth="True"){
  colors =  c("black","red")
 
  f1_y_breaks = c(-0.2,0,0.2)
  f2_y_breaks = c(round(min(v_true_df$v),1), 0, round(max(v_true_df$v),1))
  f3_y_breaks = c(round(min(w_true_df$w),1), round(min(w_true_df$w)/2+max(w_true_df$w)/2,1),round(max(w_true_df$w),1))
  
  if(Meth=="True"){
    f1_y_breaks = c(-1.0,-0.5, 0.0, 0.5, 1.0)
    f2_y_breaks = c(-1.0,-0.5, 0.0, 0.5, 1.0)
    f3_y_breaks = c( 0.0, 0.5, 1.0)
  }
  
  # plot u
  f1 = ggplot(data=u_true_df, aes(x=x, y = u, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = paste(Meth,"u"),x=NULL,y=NULL)+guides(color=FALSE) + theme_bw() + 
    scale_y_continuous(breaks=as.numeric(f1_y_breaks))
  
  # plot v
  f2 = ggplot(data=v_true_df, aes(x=x, y = v, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = paste(Meth,"v"),x=NULL,y=NULL)+guides(color=FALSE) + theme_bw() + 
    scale_y_continuous(breaks=f2_y_breaks)
  
  # plot w
  f3 = ggplot(data=w_true_df, aes(x=x, y = w, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = paste(Meth,"w"),x=NULL,y=NULL)+guides(color=FALSE) + theme_bw() + 
    scale_y_continuous(breaks=f3_y_breaks)
  
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  
  f1 = f1 + my_theme
  f2 = f2 + my_theme
  f3 = f3 + my_theme
  
  if(Meth=="l_infty_wsPLS"){
    # f1 = f1 + labs(title = expression(paste("l",infinity,"-wsPLS u",sep="")))
    # f2 = f2 + labs(title = expression(paste("l",infinity,"-wsPLS v",sep="")))
    # f3 = f3 + labs(title = expression(paste("l",infinity,"-wsPLS w",sep="")))
    f1 = f1 + labs(title =TeX("$l_\\infty-wsPLS~~u$"))
    f2 = f2 + labs(title =TeX("$l_\\infty-wsPLS~~v$"))
    f3 = f3 + labs(title =TeX("$l_\\infty-wsPLS~~w$"))
  }

  if(Meth=="l2-wsPLS"){
    f1 = f1 + labs(title =TeX("$l_2-wsPLS~~u$"))
    f2 = f2 + labs(title =TeX("$l_2-wsPLS~~v$"))
    f3 = f3 + labs(title =TeX("$l_2-wsPLS~~w$"))
  }
  # PLS
  # PMD-sPLS
  # l0-sPLS
  # l2-wsPLS
  # l_infty_wsPLS
  
  fig = plot_grid(f1, f2, f3, ncol=3, nrow =1, align="hv", labels =c("","",""))
  return(fig)
}

get_figs_uv = function(u_true_df,v_true_df,Meth="True"){
  colors =  c("black","red")
  
  f1_y_breaks = c(-0.2,0,0.2)
  f2_y_breaks = c(round(min(v_true_df$v),1), 0, round(max(v_true_df$v),1))
  
  if(Meth=="True"){
    f1_y_breaks = c(-1.0, 0.0, 1.0)
    f2_y_breaks = c(-1.0, 0.0, 1.0)
  }
  
  # plot u
  f1 = ggplot(data=u_true_df, aes(x=x, y = u, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = paste(Meth,"u"),x=NULL,y=NULL)+guides(color=FALSE) + theme_bw() + 
    scale_y_continuous(breaks=f1_y_breaks)
  
  # plot v
  f2 = ggplot(data=v_true_df, aes(x=x, y = v, color = group)) + geom_point(size=1,shape=16)+scale_color_manual(values = colors)+
    labs(title = paste(Meth,"v"),x=NULL,y=NULL)+guides(color=FALSE) + theme_bw() + 
    scale_y_continuous(breaks=f2_y_breaks)
  
  my_theme =  theme(
    plot.title   = element_text(size = rel(1.5)),
    axis.title.y = element_text(size = rel(1.5)),
    axis.title.x = element_text(size = rel(1.5)),
    axis.text=element_text(size = rel(1.2),colour="black"),
    axis.line  = element_line(size = rel(1.2),colour = "black"),
    axis.ticks=element_line(size=1)) 
  
  f1 = f1 + my_theme
  f2 = f2 + my_theme
  
  if(Meth=="l0-sPLS"){
    f1 = f1 + labs(title =TeX("$l_0-sPLS~~u$"))
    f2 = f2 + labs(title =TeX("$l_0-sPLS~~v$"))
  }
  
  fig = plot_grid(f1, f2, NULL, ncol=3, nrow =1, align="hv", labels =c("",""))
  return(fig)
}