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


# ------------------------------------------------------------------------------
# We use wsPLS algorithm on the above datasets
source('PALM_wsPLS.R')
out = PALM_wsPLS(X, Y, ku, kv, kw, Lc=0.1, nstart=5, seed0=1)

u_order2 = order(-abs(out$u))
v_order2 = order(-abs(out$v))
w_order2 = order(-abs(out$w))

X_order = X[w_order2,u_order2]
Y_order = Y[w_order2,v_order2] 


# ------------------------------------------------------------------------------
png(filename = "Fig1_original_data_matrices.png",width=9.15, height=4.5, units = "in",res=600)
ht1 = Heatmap(X, 
              column_title = "X",
              cluster_rows = F,
              cluster_columns = F, 
              column_names_side = "top",
              row_names_side = "left",
              show_heatmap_legend = FALSE, name="foo1")
ht2 = Heatmap(Y, 
              column_title = "Y",
              cluster_rows = F, 
              cluster_columns = F, 
              column_names_side = "top",
              row_names_side = "left",
              show_heatmap_legend = FALSE, name="foo2")
ht_list = ht1 + ht2

ht = draw(ht_list, auto_adjust = T)

decorate_heatmap_body("foo1", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
})

decorate_heatmap_body("foo2", {
  grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
})
dev.off()

# ------------------------------------------------------------------------------
png(filename = "Fig2_order_data_matrices.png",width=9.15, height=4.5, units = "in",res=600)
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
dev.off()

# ------------------------------------------------------------------------------
# png(filename = "Fig3_figures.png",width=9.15*2, height=4.5, units = "in",res=600)
# ht_list = ht1 + ht2 + ht3 + ht4
# ht = draw(ht_list, auto_adjust = T)
# 
# decorate_heatmap_body("foo1", {
#   grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
#   grid.rect(x = unit(0/p, "npc"), width = unit(ku/p, "npc"),
#             y = unit((n-kw)/n, "npc"), height = unit(kw/n, "npc"),
#             hjust = 0, vjust = 0,
#             gp=gpar(fill = "transparent", col = "red", lwd = 4))
# })
# 
# decorate_heatmap_body("foo2", {
#   grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))
#   grid.rect(x = unit(0/q, "npc"), width = unit(kv/q, "npc"),
#             y = unit((n-kw)/n, "npc"), height = unit(kw/n, "npc"),
#             hjust = 0, vjust = 0,
#             gp=gpar(fill = "transparent", col = "red", lwd = 4))
# })
# dev.off()
# dev.off()
# ------------------------------------------------------------------------------