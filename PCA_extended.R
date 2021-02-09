library(neuprintr)
library(nat)
library(ggpubr)
library(ggfortify)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(dplyr)


### small functions ###

PCAanalysis_metu <- function(a, namingplot){
  # use the nonzero columns 
  pdf(namingplot)
  length( a[,apply(a, 2, function(x) !all(x==0))] )
  nonz <- a[,apply(a, 2, function(x) !all(x==0))] 
  # PCA scaling the matrix to unit variance - seems to have the same results as standardization which makes sense as PCA looks at variation, thus variance 
  # somepca <- PCA(a, scale.unit = TRUE, graph = FALSE)
  # print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity scaled matrix") )
  # print(fviz_eig(somepca, addlabels = TRUE) )
  # PCA without scaling the matrix
  somepca <- PCA(a, scale.unit = FALSE, graph = FALSE)
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity not scaled matrix") )
  print(fviz_eig(somepca, addlabels = TRUE, title="Scree plot not scaled matrix") )
  # PCA with standardization, i.e. subtracting the mean and unit variance
  somepca <- PCA(scale(nonz, center = TRUE, scale = TRUE), scale.unit = FALSE, graph = FALSE)
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity stand. matrix") )
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity stand. matrix", axes=c(1,3)) )
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity stand. matrix", axes=c(2,3)) )
  print(fviz_eig(somepca, addlabels = TRUE, title="Scree plot stand. matrix") )
  # heatmap without the all zero columns 
  print(heatmap(as.matrix(nonz), scale = "col", Rowv = NA, Colv = NA))
  nonzcol <- apply(a, 2, function(x) !all(x==0))
  dev.off()
  return(nonzcol)
}

#### data loading and preprocessing ####

## connection to neuprint API for dataset = "hemibrain:v1.1" by calling another script so as not to reveal token
source("/home/alikiz/Documents/WernerLab/code/loadhemibrain.R")

## getting some info for our cells of interest (MC61 & MC64 aka “MeTu") including the bodyIDs
MC61.info = neuprint_search("MC61") 
MC64.info = neuprint_search("MC64")
MeTu.info = rbind(MC61.info,MC64.info)

### Lets get some info for the TuBu neurons in the hemibrain data set
TuBu.info = neuprint_search(".*TuBu.*")
TuBu.info = TuBu.info[TuBu.info$name != "TuBu_L",] # removing all TuBus from the left hemisphere, 74 TuBu neurons in the right hemisphere

## This is how the grouping contained in the csv file was build by Emil
## Looking for all postsynaptic cells of MeTus
## TODO: Why the only the postsynaptic?
# MeTu.con = neuprint_connection_table(MeTu.info$bodyid,prepost = "POST")
# MeTu.con$bodyid_name = neuprint_get_neuron_names(MeTu.con$bodyid)
# MeTu.con$partner_name = neuprint_get_neuron_names(MeTu.con$partner)
# 
# ## focusing only on MeTu > TuBu connections
# MeTu.TuBu.con <- MeTu.con[grep("TuBu", MeTu.con$partner_name),]
# 
# ### finding the strongest connecteted tubu_group for each MeTu
# MeTu_group = MeTu.TuBu.con %>%
#   group_by(bodyid, subgroup) %>%
#   summarize(sum_weight = sum(weight)) %>%
#   top_n(1,sum_weight)
# 
# ### looking for the rare case that a MeTu cell is equally strong connected to two differnet subgroups
# MeTu_group %>% group_by(bodyid) %>% filter(n()>1)
# 
# ## there is one neuron 5812986247 which is equally strong connected to sub3 and sub4 (3 synapses each), lets
# ## exlcude this neuron for now from our further analysis
# MeTu_group = MeTu_group[!(MeTu_group$bodyid == "5812986247"),]
# 
# ### How many MeTus are in each subgroup?
# table(MeTu_group$subgroup)



### Let's get the connectivity matrix for "MeTu > TuBu_Right” connections 
matrix = neuprint_get_adjacency_matrix(inputids = MeTu.info$bodyid, outputids = TuBu.info$bodyid)
colnames(matrix) = TuBu.info$name
rownames(matrix) = MeTu.info$bodyid

# get the neurons and their corresponding group according to precalculated MeTu-TuBu connectivity 
group <- read.csv(file = "/home/alikiz/Documents/WernerLab/hemibrain/MeTu_grouping.csv")

## get the matrix for the MeTu rows: 466 rows
MeTuconnect = matrix[rownames(matrix) %in% group$bodyid,]
## check if columns contain nonzero entries all over the column --> all columns contain nonzero features 
table(apply(MeTuconnect, 2, function(x) !all(x==0)))
heatmap(MeTuconnect, scale = "col", Rowv = NA, Colv = NA)
## create the classes with an integer encoding as following (created in alphabetic order):
## sub1_a: 1;  sub1_b: 2; sub2: 3; sub3: 4;  sub4: 5; BUT we put sub1_a and sub1_b together to sub1: 2
subgr <- as.numeric(factor(group$subgroup))
names(subgr) = group$bodyid
# TODO: fix the cbind
matMeTu = cbind(MeTuconnect, subgr)
matMeTu = as.data.frame(matMeTu)
## convert the type of subgr so as to use the PCA library later on, merge sub1_a and sub1_b to sub1: 2
matMeTu$subgr <- as.character(matMeTu$subgr)
matMeTu$subgr[matMeTu$subgr == 1] <- 2

## PCA analysis: scaled and not scaled --> both hold low percentage of explained variance in the first 5 PCs (25% and 40% respectively)
somepca <- PCA(MeTuconnect, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "#E7B800", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group")
fviz_eig(somepca, addlabels = TRUE)


somepca <- PCA(MeTuconnect, scale.unit = FALSE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group")
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group",axes = c(1,3))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(2,3))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(1,4))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$subgr , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(1,5))
fviz_eig(somepca, addlabels = TRUE)
fviz_pca_var(res.pca, col.var = "black")  # for correlation circle of the variables, important: variables need different names
# more: contribution of variables to PCi, variables according to their contributions to the principal components (dimension description)


## Iterating over the 4 subcategories, we will extract the data for each category and perform PCA on one category per time
## Firstly, let us set a variabl
print("analyze aggregated group sube for the column containing the grouping so as to exclude this column from the analysis
colgroup = 75
## Category sub1:1")
nonzcol2 <- PCAanalysis_metu(matMeTu[which(matMeTu$subgr == 2),-colgroup], "/home/alikiz/Documents/WernerLab/plots/groupsub1.pdf")
table(nonzcol2)

## Category sub2:
print("analyze group sub2")
nonzcol3 <- PCAanalysis_metu(matMeTu[which(matMeTu$subgr == 3),-colgroup], "/home/alikiz/Documents/WernerLab/plots/groupsub2.pdf")
table(nonzcol3)

## Category sub3:
print("analyze group sub3")
nonzcol4 <- PCAanalysis_metu(matMeTu[which(matMeTu$subgr == 4),-colgroup], "/home/alikiz/Documents/WernerLab/plots/groupsub3.pdf")
table(nonzcol4)

## Category sub4:
print("analyze group sub4")
nonzcol5 <- PCAanalysis_metu(matMeTu[which(matMeTu$subgr == 5),-colgroup], "/home/alikiz/Documents/WernerLab/plots/groupsub4.pdf")
table(nonzcol5)

## Category sub4:
print("analyze group 3 and 4")
nzcol5 <- PCAanalysis_metu(matMeTu[which(matMeTu$subgr == 5 & matMeTu$subgr == 4),-colgroup], "/home/alikiz/Documents/WernerLab/plots/groupsub4and3.pdf")
table(nonzcol5)

# # new matrix: comb, combination of both 
# comb = cbind(connect_mat, mat_both)
# comb_label = cbind(mydata, mat_both)
# comb_pca = PCA(new, scale.unit = TRUE, graph = FALSE)
# fviz_pca_ind(comb_pca, geom.ind = "point", col.ind= new_groups$label, palette = c("#00AFBB", "#E7B800"), axes.linetype = "blank", legend.title="Group")
# fviz_eig(comb_pca, addlabels = TRUE)

# TuBu.info = neuprint_search(".*TuBu.*")
# m =neuprint_get_adjacency_matrix(inputids = MeTu.info$bodyid, outputids = TuBu.info$bodyid)
# colnames(m) = TuBu.info$name
# rownames(m) = MeTu.info$bodyid
# meeetuu = m[rownames(m) %in% group$bodyid,]
