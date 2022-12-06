library(ggplot2)
library(ggpubr)
library(latex2exp)

load("Data_Env_RData/tidyEnvData.RData")

AW_mat = env[c(paste("AW",1:9,sep="_")),]
RW_mat = env[c(paste("RW",1:9,sep="_")),]

diff_mat = AW_mat - RW_mat
################################################################################
envName = colnames(diff_mat)

# Removal of Dissolved Oxygen (DO) due to experimental error 
envName = envName[envName!="Do"]

pvalue = data.frame(name=envName, p=rep(1,length(envName)))  

for(i in 1:length(envName)){
  out = wilcox.test(diff_mat[,envName[i]], alternative = "two.sided")
  pvalue[i,2] = out$p.value
}

pvalue$p.adjust = p.adjust(pvalue$p, "BH")


Diff_res = pvalue
Diff_res = Diff_res[order(Diff_res$p),]
Diff_res$if.sig = 0
Diff_res$if.sig[which(Diff_res$p.adjust<0.05)] = 1
################################################################################
write.csv(Diff_res, file = "F3_Diff_res_plus.csv")


################################################################################
dat = Diff_res
dat$log_adj_p = -log10(dat$p.adjust)
dat$group = "Not Significant"
dat$group[dat$p.adjust<0.05] = "Significant"
#"name","p","p.adjust","if.sig","log_adj_p","group"  
colnames(dat) = c("factor",
                  "p_value","adj_p_value","if.sig","negative_log_adj_p","group") 

dat = dat[order(-dat$p_value),]
dat$factor = factor(dat$factor, levels = dat$factor)
dat$group = factor(dat$group,levels=c("Significant","Not Significant"))

f = ggplot(dat, aes(x=factor, y=negative_log_adj_p, fill=group))+
  geom_bar(stat="identity") 
# Horizontal box plot
f = f +  geom_hline(aes(yintercept=-log10(0.05)),  colour="red", linetype="dashed")+
  coord_flip() + scale_y_continuous(expand= c(0.01,0))

f = f +  theme(axis.line  = element_line(size=rel(1.2)),
               axis.ticks = element_line(size=rel(1.2)),
               axis.title = element_text(size = rel(1.3)),
               axis.text  = element_text(colour = "black", size=rel(1.1)),
               axis.text.x = element_text(angle = 60, hjust = 1),
               legend.position = c(.85, .1),
               legend.title =element_text(size=1, colour="white"), 
               legend.text=element_text(size=rel(1.2), colour="black")) 

f = f + labs(y = TeX("$-log_{10}(adj-p-value)$")) +
  labs(x = "Physicochemical factor") +
  labs(title = "Differential expression analysis of physicochemical factors")
  title = "This is title"

################################################################################
# change name of factors
# from down to up 

# dat$factor
# [1] Mo          Ba          Temperature NO3-        Eh          EC         
# [7] Cl-         Pb          Cu          NO2-        PO          Zn         
# [13] pH          Ca          Na          K           Al          Cd         
# [19] Fe          U           F-          Mg          As          Cr         
# [25] Mn          SO42-
  
f + scale_x_discrete(labels=c(
"Mo",          
"Ba",          
"Temperature", 
TeX("$NO_3^{-}$"), # "NO3-"        
"Eh",          
"EC",         
TeX("$Cl^{-}$"), # "Cl-"         
"Pb",          
"Cu", 
TeX("$NO_2^{-}$"), # "NO2-"       
"PO",          
"Zn",         
"pH",          
"Ca",          
"Na",          
"K",           
"Al",          
"Cd",         
"Fe",          
"U",           
TeX("$F^{-}$"), # "F-"           
"Mg",          
"As",          
"Cr",         
"Mn",
TeX("$SO_4^{2-}$")
))     
  
#f + scale_fill_manual(values=c("#56B4E9",'lightgray')) 
ggsave(file="F3_diff_bar+.svg", units = c("in"), width=7, height=7)
ggsave(file="F3_diff_bar+.pdf", units = c("in"), width=7, height=7)
ggsave(file="F3_diff_bar+.png", units = c("in"), width=7, height=7,dpi = 300)

