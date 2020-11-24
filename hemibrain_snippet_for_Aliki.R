library(neuprintr)
library(nat)
library(ggpubr)


## connection to neuprint API
conn = neuprint_login(server= "https://neuprint.janelia.org/",
                      token= "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJlbWFpbCI6InphdmFyb3BvdWxvdS5hbGlraUBnbWFpbC5jb20iLCJsZXZlbCI6Im5vYXV0aCIsImltYWdlLXVybCI6Imh0dHBzOi8vbGg2Lmdvb2dsZXVzZXJjb250ZW50LmNvbS8tQWJDLWx0S19mRzgvQUFBQUFBQUFBQUkvQUFBQUFBQUFBQUEvQU1adXVjay1hR1ViNjZ5UnVmVHpiaXhodWxrcVJ3blMxdy9zOTYtYy9waG90by5qcGc_c3o9NTA_c3o9NTAiLCJleHAiOjE3ODM5MDQ2NjJ9.YGVe1LLz_a_zWroPT87oLygD_Rcxom0gRIOMGI8Rz2w",
                      dataset = "hemibrain:v1.1" )


## bodyIDs from Aljoscha's cluster analysis for the two DRAMeTu subgroups 
DRAMeTu_sub1_ID = c(1078097609, 1201227161, 1171548392, 1076729226, 1202250779, 1077756950, 5813016373, 1198127720, 1078102468,
                    1108796315, 1108446169, 5812987715, 1077761540, 1170179426)

DRAMeTu_sub2_ID = c(1139826933, 1077756731, 1078097817, 1198809632, 1076729233, 1078102490, 1107423158, 1078102267, 1139485872,
                    1077761532, 1078102350, 1139148877, 1140172133, 1170861854, 1201559521, 1077424293, 1201572715, 1076729262,
                    1201210029, 5901210109, 1199155111, 1077424204, 1198809646, 1171211215, 5813104106, 1078102400, 1077424267,
                    1139831143, 1139831106, 1200873417, 1232939935, 1140517287, 1108450770, 5812986562, 1078102419, 1169834028,
                    1202250762, 1201905111, 1046039977)

## getting 'meta-data'-data.frame for cells
DRAMeTu_sub1 = neuprint_get_meta(DRAMeTu_sub1_ID)
DRAMeTu_sub2 = neuprint_get_meta(DRAMeTu_sub2_ID)
MC61 = neuprint_search(".*MC61.*")
MC64 = neuprint_search(".*MC64.*")
TuBu = neuprint_search(".*TuBu.*")
TuTu = neuprint_search(".*TuTu.*")
DN1p = neuprint_search(".*DN1p.*")
Ring = neuprint_search(".*ring.*")

#neupritn_get_synapsis
## renaming 'type' and 'name' for the DRAMeTu subtypes
DRAMeTu_sub1$name = "DRAMeTu_sub1"
DRAMeTu_sub1$type = "DRAMeTu"

DRAMeTu_sub2$name = "DRAMeTu_sub2"
DRAMeTu_sub2$type = "DRAMeTu"

## creating data.frame for MC61 with DRAMeTubodyIDs (all.MC61) and with them removed (only.MC61)
all.MC61 = MC61
only.MC61 = MC61[ !(MC61$bodyid %in% DRAMeTu_sub1$bodyid), ]
only.MC61 = only.MC61[ !(only.MC61$bodyid %in% DRAMeTu_sub2$bodyid), ]

## subsetting and rbinding the 'meta-data'-data.frames for further analysis
all.AVP = rbind(DRAMeTu_sub1,DRAMeTu_sub2,only.MC61,MC64,TuBu,TuTu,DN1p)

TuBu1 = subset(all.AVP,name == "TuBu01_ABU_R")
TuBu2 = subset(all.AVP,name == "TuBu02_IBU_R")
TuBu3 = subset(all.AVP,name == "TuBu03_IBU_R")
TuBu4 = subset(all.AVP,name == "TuBu04_IBU_R")
TuBu5 = subset(all.AVP,name == "TuBu05_SBU_R")
TuBu6 = subset(all.AVP,name == "TuBu06_SBU_R")
TuBu7 = subset(all.AVP,name == "TuBu07_SBU_R")
TuBu8 = subset(all.AVP,name == "TuBu08_SBU_R")
TuBu9 = subset(all.AVP,name == "TuBu09_SBU_R")
TuBu10 = subset(all.AVP,name == "TuBu10_SBU_R")
TuTuB_a_L = subset(all.AVP,name == "TuTuB_a(ADL19)_L")
TuTuB_a_R = subset(all.AVP,name == "TuTuB_a_R")
TuTuB_b_L = subset(all.AVP,name == "TuTuB_b(ADL19)_L")
TuTuB_b_R = subset(all.AVP,name == "TuTuB_b_R")
DN1pB = subset(all.AVP,name == "DN1pB_R")

TuTuA = subset(all.AVP,type == "TuTuA")
TuBuL = subset(all.AVP,name == "TuBu_L")
DN1pA = subset(all.AVP,name == "DN1pA_R")

all.POL = rbind(DRAMeTu_sub1,DRAMeTu_sub2,TuBu1,TuBu6,TuTuB_b_L,TuTuB_b_R,DN1pB)
TuBu = rbind(TuBu1,TuBu2,TuBu3,TuBu4,TuBu5,TuBu6,TuBu7,TuBu8,TuBu9,TuBu10)
MeTu = rbind(DRAMeTu_sub1,DRAMeTu_sub2,only.MC61,MC64)

## removing 'TuTuA', 'TuBuL' & 'DN1pA' neurons for further analysis
all.AVP = all.AVP[ !(all.AVP$bodyid %in% TuTuA$bodyid),]
all.AVP = all.AVP[ !(all.AVP$bodyid %in% TuBuL$bodyid),]
all.AVP = all.AVP[ !(all.AVP$bodyid %in% DN1pA$bodyid),]

## plotting average downstream/upstream synapses for all.AVP
p = ggerrorplot(all.AVP, x = "name", y = "downstream", add = "jitter", desc_stat = "mean_sd", color = "name", add.params = list(color = "darkgray"))
p = p + rotate_x_text()
ggpar(p,legend = "none")

p = ggerrorplot(all.AVP, x = "name", y = "upstream", add = "jitter", desc_stat = "mean_sd", color = "name", add.params = list(color = "darkgray"))
p = p + rotate_x_text()
ggpar(p,legend = "none")

## reading in the neurons with 'neuprint_read_neurons()'
only.MC61.n = neuprint_read_neurons(only.MC61$bodyid)
MC64.n = neuprint_read_neurons(MC64$bodyid)
DRAMeTu_sub1.n = neuprint_read_neurons(DRAMeTu_sub1$bodyid)
DRAMeTu_sub2.n = neuprint_read_neurons(DRAMeTu_sub2$bodyid)

DRAMeTu_sub1.n[,'type'] = 'DRAMeTu'
DRAMeTu_sub2.n[,'type'] = 'DRAMeTu'
DRAMeTu_sub1.n[,'name'] = 'DRAMeTu_sub1'
DRAMeTu_sub2.n[,'name'] = 'DRAMeTu_sub2'

TuBu.n = neuprint_read_neurons(TuBu$bodyid)
# TuBu1.n = neuprint_read_neurons(TuBu1$bodyid)
# TuBu2.n = neuprint_read_neurons(TuBu2$bodyid)
# TuBu3.n = neuprint_read_neurons(TuBu3$bodyid)
# TuBu4.n = neuprint_read_neurons(TuBu4$bodyid)
# TuBu5.n = neuprint_read_neurons(TuBu5$bodyid)
# TuBu6.n = neuprint_read_neurons(TuBu6$bodyid)
# TuBu7.n = neuprint_read_neurons(TuBu7$bodyid)
# TuBu8.n = neuprint_read_neurons(TuBu8$bodyid)
# TuBu9.n = neuprint_read_neurons(TuBu9$bodyid)
# TuBu10.n = neuprint_read_neurons(TuBu10$bodyid)


TuTuB_a_L.n = neuprint_read_neurons(TuTuB_a_L$bodyid)
TuTuB_a_R.n = neuprint_read_neurons(TuTuB_a_R$bodyid)
TuTuB_b_L.n = neuprint_read_neurons(TuTuB_b_L$bodyid)
TuTuB_b_R.n = neuprint_read_neurons(TuTuB_b_R$bodyid)

DN1pB.n = neuprint_read_neurons(DN1pB$bodyid)

## reading in relavant ROI's. Lookup all ROI's with neuprint_ROIs()
BU.mesh = neuprint_ROI_mesh(roi = "BU(R)")
EB.mesh = neuprint_ROI_mesh(roi = "EB")
AOTU.mesh = neuprint_ROI_mesh(roi = "AOTU(R)")
ME.mesh = neuprint_ROI_mesh(roi = "ME(R)")

## transforming ROI's into hxsurf 
BU.mesh.surf = as.hxsurf(BU.mesh)
EB.mesh.surf = as.hxsurf(EB.mesh)
AOTU.mesh.surf = as.hxsurf(AOTU.mesh)
ME.mesh.surf = as.hxsurf(ME.mesh)

## pruning different neurons in 'AOTU' or 'BU'
prune_AOTU_only.MC61.n = prune_in_volume(only.MC61.n,surf = AOTU.mesh.surf)
prune_AOTU_MC64.n = prune_in_volume(MC64.n,surf = AOTU.mesh.surf)
prune_AOTU_DRAMeTu_sub1.n = prune_in_volume(DRAMeTu_sub1.n,surf = AOTU.mesh.surf)
prune_AOTU_DRAMeTu_sub2.n = prune_in_volume(DRAMeTu_sub2.n,surf = AOTU.mesh.surf)
prune_AOTU_TuBu.n = prune_in_volume(TuBu.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu1.n = prune_in_volume(TuBu1.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu2.n = prune_in_volume(TuBu2.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu3.n = prune_in_volume(TuBu3.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu4.n = prune_in_volume(TuBu4.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu5.n = prune_in_volume(TuBu5.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu6.n = prune_in_volume(TuBu6.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu7.n = prune_in_volume(TuBu7.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu8.n = prune_in_volume(TuBu8.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu9.n = prune_in_volume(TuBu9.n,surf = AOTU.mesh.surf)
# prune_AOTU_TuBu10.n = prune_in_volume(TuBu10.n,surf = AOTU.mesh.surf)
prune_AOTU_TuTuB_a_L.n = prune_in_volume(TuTuB_a_L.n,surf = AOTU.mesh.surf)
prune_AOTU_TuTuB_a_R.n = prune_in_volume(TuTuB_a_R.n,surf = AOTU.mesh.surf)
prune_AOTU_TuTuB_b_L.n = prune_in_volume(TuTuB_b_L.n,surf = AOTU.mesh.surf)
prune_AOTU_TuTuB_b_R.n = prune_in_volume(TuTuB_b_R.n,surf = AOTU.mesh.surf)
prune_AOTU_DN1pB.n = prune_in_volume(DN1pB.n,surf = AOTU.mesh.surf)

prune_BU_TuBu.n = prune_in_volume(TuBu.n,surf = BU.mesh.surf)
# prune_BU_TuBu1.n = prune_in_volume(TuBu1.n,surf = BU.mesh.surf)
# prune_BU_TuBu2.n = prune_in_volume(TuBu2.n,surf = BU.mesh.surf)
# prune_BU_TuBu3.n = prune_in_volume(TuBu3.n,surf = BU.mesh.surf)
# prune_BU_TuBu4.n = prune_in_volume(TuBu4.n,surf = BU.mesh.surf)
# prune_BU_TuBu5.n = prune_in_volume(TuBu5.n,surf = BU.mesh.surf)
# prune_BU_TuBu6.n = prune_in_volume(TuBu6.n,surf = BU.mesh.surf)
# prune_BU_TuBu7.n = prune_in_volume(TuBu7.n,surf = BU.mesh.surf)
# prune_BU_TuBu8.n = prune_in_volume(TuBu8.n,surf = BU.mesh.surf)
# prune_BU_TuBu9.n = prune_in_volume(TuBu9.n,surf = BU.mesh.surf)
# prune_BU_TuBu10.n = prune_in_volume(TuBu10.n,surf = BU.mesh.surf)


## plotting stuff
##

matrix.BU.view = structure(c(0.33, 0.57, 0.75, 0,
                             -0.71, -0.37, 0.60, 0,
                             0.62, -0.73, 0.28, 0,
                             0, 0, 0, 1), 
                           .Dim = c(4L, 4L))

matrix.AOTU.dorsal.view = structure(c(0.76, 0.41, -0.49, 0,
                                      0.39, -0.91, -0.16, 0,
                                      -0.52, -0.06, -0.85, 0,
                                      0, 0, 0, 1), 
                                    .Dim = c(4L, 4L))

matrix.AOTU.anterior.view = structure(c(0.89, -0.30, -0.33, 0,
                                        0.37, 0.05, 0.93, 0,
                                        -0.27, -0.95, 0.15, 0,
                                        0, 0, 0, 1), 
                                      .Dim = c(4L, 4L))

window.BU.view = c(104L, 1095L, 947L, 1824L)
window.AOTU.dorsal.view = window.BU.view 
window.AOTU.anterior.view = window.BU.view 

zoom.BU.view = 1
zoom.AOTU.dorsal.view = zoom.BU.view
zoom.AOTU.anterior.view = zoom.BU.view

nopen3d(userMatrix = matrix.BU.view,
        zoom = zoom.BU.view, 
        windowRect = window.BU.view)

nopen3d(userMatrix = matrix.AOTU.dorsal.view,
        zoom = zoom.AOTU.dorsal.view, 
        windowRect = window.AOTU.dorsal.view)

nopen3d(userMatrix = matrix.AOTU.anterior.view,
        zoom = zoom.AOTU.anterior.view, 
        windowRect = window.AOTU.anterior.view)

open3d()

plot3d(TuBu.n, col = type)
plot3d(prune_AOTU_TuBu.n, type == 'TuBu01', col = "darkgreen")
plot3d(prune_AOTU_TuBu.n, type == 'TuBu06', col = "darkblue")
plot3d(prune_BU_TuBu.n, col = type)
plot3d(prune_BU_TuBu.n, type == 'TuBu01', col = "darkgreen",WithLine = F,WithNodes = T)
plot3d(prune_BU_TuBu.n, type != 'TuBu01', col = "darkgrey",WithLine = F,WithNodes = T)

# plot3d(prune_BU_TuBu1.n,col = "darkgreen")
# plot3d(prune_BU_TuBu2.n,col = "darkgrey")
# plot3d(prune_BU_TuBu3.n,col = "darkgrey")
# plot3d(prune_BU_TuBu4.n,col = "darkgrey")
# plot3d(prune_BU_TuBu5.n,col = "darkgrey")
# plot3d(prune_BU_TuBu6.n,col = "blue")
# plot3d(prune_BU_TuBu7.n,col = "darkgrey")
# plot3d(prune_BU_TuBu8.n,col = "darkgrey")
# plot3d(prune_BU_TuBu9.n,col = "darkgrey")
# plot3d(prune_BU_TuBu9.n,col = "darkgrey")
# plot3d(prune_BU_TuBu10.n,col = "darkgrey")
# 
# plot3d(prune_AOTU_TuBu1.n,col = "darkgreen")
# plot3d(prune_AOTU_TuBu2.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu3.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu4.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu5.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu6.n,col = "blue")
# plot3d(prune_AOTU_TuBu7.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu8.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu9.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu9.n,col = "darkgrey")
# plot3d(prune_AOTU_TuBu10.n,col = "darkgrey")

plot3d(AOTU.mesh, add = TRUE, alpha = 0.1)
plot3d(BU.mesh, add = TRUE, alpha = 0.1)
#plot3d(EB.mesh, add = TRUE, alpha = 0.1)
#plot3d(ME.mesh, add = TRUE, alpha = 0.1)

# plot3d(DRAMeTu_sub1.n,col = "green", lwd = 2)
# plot3d(DRAMeTu_sub2.n,col = "yellow", lwd = 2)
# plot3d(TuBu1.n,col = "blue", lwd = 2)
# plot3d(TuBu2.n,col = "green", lwd = 2)
# plot3d(TuBu6.n,col = "red", lwd = 2)
# plot3d(MC61.n,col = "red", lwd = 2)

## connectivity matrix
MeTu.TuBu.adj = neuprint_get_adjacency_matrix(outputids = TuBu$bodyid, inputids = MeTu$bodyid)
colnames(MeTu.TuBu.adj) = TuBu$name 
rownames(MeTu.TuBu.adj) = MeTu$name
MeTu.TuBu.adj.comp = MeTu.TuBu.adj
MeTu.TuBu.adj.comp = t(apply(t(MeTu.TuBu.adj.comp), 2, function(x) tapply(x, colnames(MeTu.TuBu.adj.comp), sum, na.rm = TRUE)))
MeTu.TuBu.adj.comp = apply(MeTu.TuBu.adj.comp, 2, function(x) tapply(x, rownames(MeTu.TuBu.adj.comp), sum, na.rm = TRUE))

sampled.matrix <- MeTu.TuBu.adj.comp[sample(nrow(MeTu.TuBu.adj.comp)),]
sampled.matrix <- sampled.matrix[,sample(ncol(sampled.matrix))]

heatmap(MeTu.TuBu.adj, Rowv = NA, Colv = NA)
heatmap(MeTu.TuBu.adj.comp,col = col)

heatmap(sampled.matrix, Rowv = NA, Colv = NA)
heatmap(sampled.matrix)

scaled.MeTu.TuBu.adj = scale(MeTu.TuBu.adj)

d <- dist(scaled.MeTu.TuBu.adj, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1)
sub_grp <- cutree(hc1, k = 4)
rect.hclust(hc1, k = 4, border = 2:5)
table(sub_grp)
des = which("3" == sub_grp)

maybe = names(des)
setequal(maybe,DRAMeTu_stats$bodyid)

heatmap(df5, col = col,Rowv = as.dendrogram(hc1), Colv = NA)
colnames(df5) = all.TuBu$name
df6 = t(apply(t(df5), 2, function(x) tapply(x, colnames(df5), sum, na.rm = TRUE)))
all.POL.adj = neuprint_get_adjacency_matrix(bodyids = all.POL$bodyid)
rownames(all.POL.adj) = colnames(all.POL.adj) = all.POL$name

heatmap(all.POL.adj)


all.POL.adj.comp = all.POL.adj
rownames(all.POL.adj.comp) = colnames(all.POL.adj.comp) = all.POL$name
all.POL.adj.comp = t(apply(t(all.POL.adj.comp), 2, function(x) tapply(x, colnames(all.POL.adj.comp), sum, na.rm = TRUE)))
all.POL.adj.comp = apply(all.POL.adj.comp, 2, function(x) tapply(x, rownames(all.POL.adj.comp), sum, na.rm = TRUE))

heatmap(all.POL.adj.comp, Rowv = NA, Colv = NA, revC = T)

all.POL.adj.comp[all.POL.adj.comp < 10] <- 0
all.POL.adj[all.POL.adj < 3] <- 0

n = network(all.POL.adj,
            matrix.type = "adjacency",
            ignore.eval = FALSE,
            layout = "fruchtermanreingold",
            names.eval = "weight",
            directed = TRUE)

n = ggnetwork(n, cell.jitter = 0.75, arrow.gap = 0.01)

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = vertex.names),
             curvature = 0.05,
             arrow = arrow(length = unit(6, "pt"),
                           type = "closed")) +
  geom_nodes(aes(color = vertex.names, size = 6)) +
  geom_edgetext(aes(label = weight, color = vertex.names), fill = NA) +
  geom_nodelabel_repel(aes(color = vertex.names, label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank()+
  guides(color = FALSE, shape = FALSE, fill = FALSE, size = FALSE, linetype = FALSE) + ylab("") + xlab("")



all.AVP.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = all.AVP, FUN = mean)
all.AVP.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = all.AVP, FUN = sum)
all.AVP.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = all.AVP, FUN = length)


all.AVP.adj = neuprint_get_adjacency_matrix(bodyids = all.AVP$bodyid)
rownames(all.AVP.adj) = colnames(all.AVP.adj) = all.AVP$name

heatmap(all.AVP.adj)


all.AVP.adj.comp = all.AVP.adj
rownames(all.AVP.adj.comp) = colnames(all.AVP.adj.comp) = all.AVP$name
all.AVP.adj.comp = t(apply(t(all.AVP.adj.comp), 2, function(x) tapply(x, colnames(all.AVP.adj.comp), sum, na.rm = TRUE)))
all.AVP.adj.comp = apply(all.AVP.adj.comp, 2, function(x) tapply(x, rownames(all.AVP.adj.comp), sum, na.rm = TRUE))

heatmap(all.AVP.adj.comp)


all.AVP.subset = rbind(DRAMeTu_sub1,DRAMeTu_sub2,MC61,MC64,TuBu,TuTu,DN1p)


all.AVP.subset.adj = neuprint_get_adjacency_matrix(bodyids = all.AVP.subset$bodyid)
rownames(all.AVP.subset.adj) = colnames(all.AVP.subset.adj) = all.AVP.subset$name

heatmap(all.AVP.subset.adj, Rowv = NA, Colv = NA, revC = T)


all.AVP.subset.adj.comp = all.AVP.subset.adj
rownames(all.AVP.subset.adj.comp) = colnames(all.AVP.subset.adj.comp) = all.AVP.subset$name
all.AVP.subset.adj.comp = t(apply(t(all.AVP.subset.adj.comp), 2, function(x) tapply(x, colnames(all.AVP.subset.adj.comp), sum, na.rm = TRUE)))
all.AVP.subset.adj.comp = apply(all.AVP.subset.adj.comp, 2, function(x) tapply(x, rownames(all.AVP.subset.adj.comp), sum, na.rm = TRUE))

heatmap(all.AVP.subset.adj.comp, Rowv = NA, Colv = NA, revC = T)

all.AVP.subset.adj.comp.2 = all.AVP.subset.adj
rownames(all.AVP.subset.adj.comp.2) = colnames(all.AVP.subset.adj.comp.2) = all.AVP.subset$type
all.AVP.subset.adj.comp.2 = t(apply(t(all.AVP.subset.adj.comp.2), 2, function(x) tapply(x, colnames(all.AVP.subset.adj.comp.2), sum, na.rm = TRUE)))
all.AVP.subset.adj.comp.2 = apply(all.AVP.subset.adj.comp.2, 2, function(x) tapply(x, rownames(all.AVP.subset.adj.comp.2), sum, na.rm = TRUE))

heatmap(all.AVP.subset.adj.comp.2, Rowv = NA, Colv = NA, revC = T)

all.AVP.subset.adj.comp.2[all.AVP.subset.adj.comp.2 < 51] <- 0
all.AVP.subset.adj.comp[all.AVP.subset.adj.comp < 51] <- 0


n = network(all.AVP.subset.adj.comp,
            matrix.type = "adjacency",
            ignore.eval = FALSE,
            layout = "drl",
            names.eval = "weight",
            directed = TRUE,
            loops = TRUE)

n = network(all.AVP.subset.adj.comp.2,
            matrix.type = "adjacency",
            ignore.eval = FALSE,
            layout = "fruchtermanreingold",
            names.eval = "weight",
            directed = TRUE)

n = ggnetwork(n, cell.jitter = 0.75, arrow.gap = 0.01)

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = vertex.names),
             curvature = 0.05,
             arrow = arrow(length = unit(6, "pt"),
                           type = "closed")) +
  geom_nodes(aes(color = vertex.names, size = 6)) +
  geom_edgetext(aes(label = weight, color = vertex.names), fill = NA) +
  geom_nodelabel_repel(aes(color = vertex.names, label = vertex.names),
                       fontface = "bold", box.padding = unit(1, "lines")) +
  theme_blank()+
  guides(color = FALSE, shape = FALSE, fill = FALSE, size = FALSE, linetype = FALSE) + ylab("") + xlab("")









MC61.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC61, FUN = mean)
MC61.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC61, FUN = sum)
MC61.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC61, FUN = length)

MC64.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC64, FUN = mean)
MC64.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC64, FUN = sum)
MC64.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = MC64, FUN = length)

TuBu.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuBu, FUN = mean)
TuBu.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuBu, FUN = sum)
TuBu.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuBu, FUN = length)

TuTu.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuTu, FUN = mean)
TuTu.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuTu, FUN = sum)
TuTu.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = TuTu, FUN = length)

DN1p.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = DN1p, FUN = mean)
DN1p.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = DN1p, FUN = sum)
DN1p.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = DN1p, FUN = length)

Ring.agg.mean    = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = Ring, FUN = mean)
Ring.agg.sum     = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = Ring, FUN = sum)
Ring.agg.length  = aggregate(cbind(pre,post,upstream,downstream) ~ name, data = Ring, FUN = length)

DRAMeTu_1_ID = unique(DRA_MeTu_group1$bodyId)
DRAMeTu_2_ID = unique(DRA_MeTu_group2$bodyId)

DRAMeTu_sub1 = neuprint_read_neurons(DRAMeTu_sub1_ID)
DRAMeTu_sub2 = neuprint_read_neurons(DRAMeTu_sub2_ID)

DRAMeTu_sub1_stats = summary(DRAMeTu_sub1,include.attached.dataframe = TRUE)
DRAMeTu_sub1_stats$subtype = "sub1"

DRAMeTu_sub2_stats = summary(DRAMeTu_sub2,include.attached.dataframe = TRUE)
DRAMeTu_sub2_stats$subtype = "sub2"

DRAMeTu_stats = rbind(DRAMeTu_sub1_stats,DRAMeTu_sub2_stats)

test = c(10479161,10737978,15743226,11903992,11908710,11992932,11993030,11993064,11993314,11993352,11993382,11994448)
drametu.fafb=read.neurons.catmaid(test)

drametu.fafb.prune =nlapply(drametu.fafb, subset, X>363000, OmitFailures = T)

drametu.fafb.prune.stats = summary(drametu.fafb.prune,include.attached.dataframe = TRUE)
drametu.fafb.prune.stats$subtype = "FAFB"

n = subset(DRAMeTu_stats, select = c(branchpoints, subtype, endpoints,segments))
m = subset(drametu.fafb.prune.stats, select = c(branchpoints, subtype, endpoints,segments))


prune = rbind(n,m)

p = ggerrorplot(prune, x = "subtype", y = "branchpoints", add = "jitter", desc_stat = "mean_sd", color = "subtype", add.params = list(color = "darkgray"), palette = "jco")
#ggerrorplot(prune, x = "subtype", y = "endpoints", add = "jitter", desc_stat = "mean_sd", color = "subtype", add.params = list(color = "darkgray"), palette = "jco")

ggpar(p,legend = "none")

# ggsave("branchpoints.png",
#        width = 3,
#        height = 3,
#        dpi = 600)            
TuBu1 = neuprint_search("TuBu1.*")
TuBu4 = neuprint_search("TuBu4.*")
TuTu = neuprint_search(".*TuTu.*")
Ring = neuprint_search(".*ring.*")
R4m = neuprint_search(".*R4m.*")
ExR1 = neuprint_search(".*ExR1.*")
DN1pB = neuprint_search(".*DN1pB.*")





R4mR = neuprint_read_neurons(R4m$bodyid[6:10])
ExR1R = neuprint_read_neurons(ExR1$bodyid[5:6])
ring = neuprint_read_neurons(Ring$bodyid[1])
rings = neuprint_read_neurons(Ring$bodyid)
DRAMeTu_sub1 = neuprint_read_neurons(unique(DRA_MeTu_group1$bodyId))
DRAMeTu_sub2 = neuprint_read_neurons(unique(DRA_MeTu_group2$bodyId))
TuBu4_R = neuprint_read_neurons(TuBu4)
DN1pB_R = neuprint_read_neurons(DN1pB)
TuBu1_R = neuprint_read_neurons(TuBu1)
MC61_R = neuprint_read_neurons(MC61$bodyid)
MC64_R = neuprint_read_neurons(MC64$bodyid)



TuTuB = neuprint_search(".*TuTuB.*")
TuTuB_n = neuprint_read_neurons(TuTuB$bodyid)
