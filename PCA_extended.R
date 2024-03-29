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
  somepca <- PCA(nonz, scale.unit = FALSE, graph = FALSE)
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity not scaled matrix") )
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity not scaled matrix", axes=c(1,3)) )
  print(fviz_pca_ind(somepca, geom = "point", title= "PCA on connectivity not scaled matrix", axes=c(2,3)) )
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
source("/home/alikiz/Documents/WernetLab/code/loadhemibrain.R")

## getting some info for our cells of interest (MC61 & MC64 aka “MeTu") including the bodyIDs
MC61.info = neuprint_search("MC61") 
MC64.info = neuprint_search("MC64")
MeTu.info = rbind(MC61.info,MC64.info)

### Lets get some info for the TuBu neurons in the hemibrain data set
TuBu.info = neuprint_search(".*TuBu.*")
TuBu.info = TuBu.info[TuBu.info$name != "TuBu_L",] # removing all TuBus from the left hemisphere, 74 TuBu neurons in the right hemisphere


### Let's get the connectivity matrix for "MeTu > TuBu_Right” connections 
matrix = neuprint_get_adjacency_matrix(inputids = MeTu.info$bodyid, outputids = TuBu.info$bodyid)
colnames(matrix) = TuBu.info$name
rownames(matrix) = MeTu.info$bodyid

# get the neurons and their corresponding group according to precalculated MeTu-TuBu connectivity 
group <- read.csv(file = "MeTu_grouping.csv")

## get the matrix for the MeTu rows: 466 rows
MeTuconnect = matrix[rownames(matrix) %in% group$bodyid,]
## check if columns contain nonzero entries all over the column --> all columns contain nonzero features 
table(apply(MeTuconnect, 2, function(x) !all(x==0)))
heatmap(MeTuconnect, scale = "col", Rowv = NA, Colv = NA)
## create the connectome matrix as a data frame and add a column for the grouping labels
matMeTu = as.data.frame(MeTuconnect)
matMeTu$label = 0
## create the classes by selecting the rownames directly from the MeTuconnect matrix with the following integer encoding:
## sub1_a: 1;  sub1_b: 2; sub2: 3; sub3: 4;  sub4: 5; BUT we put sub1_a and sub1_b together to sub1: 2
matMeTu[ rownames(matMeTu) %in% group$bodyid[which(group$subgroup == "sub1_a")] , 75] = 1
matMeTu[ rownames(matMeTu) %in% group$bodyid[which(group$subgroup == "sub1_b")] , 75] = 2
matMeTu[ rownames(matMeTu) %in% group$bodyid[which(group$subgroup == "sub2")] , 75] = 3
matMeTu[ rownames(matMeTu) %in% group$bodyid[which(group$subgroup == "sub3")] , 75] = 4
matMeTu[ rownames(matMeTu) %in% group$bodyid[which(group$subgroup == "sub4")] , 75] = 5

## convert the type of subgr so as to use the PCA library later on, 
matMeTu$label <- as.character(matMeTu$label)

## PCA analysis: scaled and not scaled --> both hold low percentage of explained variance in the first 5 PCs (25% and 40% respectively)
somepca <- PCA(MeTuconnect, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "#E7B800", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group")
fviz_eig(somepca, addlabels = TRUE)

## merge sub1_a and sub1_b to sub1: 2
matMeTu$label[matMeTu$label == 1] <- 2
#table(matMeTu$label)        # to possibly check the 
somepca <- PCA(MeTuconnect, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group")
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group",axes = c(1,3))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(2,3))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(1,4))
fviz_pca_ind(somepca, geom.ind = "point", col.ind= matMeTu$label , palette = c("#00AFBB", "orange", "green", "purple"), axes.linetype = "blank", legend.title="Group", axes = c(1,5))
fviz_eig(somepca, addlabels = TRUE)
fviz_pca_var(res.pca, col.var = "black")  # for correlation circle of the variables, important: variables need different names
# more: contribution of variables to PCi, variables according to their contributions to the principal components (dimension description)


## Iterating over the 4 subcategories, we will extract the data for each category and perform PCA on one category per time
## Firstly, let us set a variabl
print("analyze aggregated groups for the column containing the grouping so as to exclude this column from the analysis")
colgroup = 75
## Category sub1
nonzcol2 <- PCAanalysis_metu(matMeTu[which(matMeTu$label == 2),-colgroup], "/home/alikiz/Documents/WernetLab/plots/groupsub1.pdf")
table(nonzcol2)

## Category sub2:
print("analyze group sub2")
nonzcol3 <- PCAanalysis_metu(matMeTu[which(matMeTu$label == 3),-colgroup], "/home/alikiz/Documents/WernetLab/plots/groupsub2.pdf")
table(nonzcol3)

## Category sub3:
print("analyze group sub3")
nonzcol4 <- PCAanalysis_metu(matMeTu[which(matMeTu$label == 4),-colgroup], "/home/alikiz/Documents/WernetLab/plots/groupsub3.pdf")
table(nonzcol4)

## Category sub4:
print("analyze group sub4")
nonzcol5 <- PCAanalysis_metu(matMeTu[which(matMeTu$label == 5),-colgroup], "/home/alikiz/Documents/WernetLab/plots/groupsub4.pdf")
table(nonzcol5)

## Category sub4:
print("analyze group 3 and 4")
nzcol5 <- PCAanalysis_metu(matMeTu[which(matMeTu$label == 5 & matMeTu$label == 4),-colgroup], "/home/alikiz/Documents/WernetLab/plots/groupsub4and3.pdf")
table(nonzcol5)

# new matrix: comb, combination of both
comb = cbind(connect_mat, mat_both)
comb_label = cbind(mydata, mat_both)
comb_pca = PCA(new, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(comb_pca, geom.ind = "point", col.ind= new_groups$label, palette = c("#00AFBB", "#E7B800"), axes.linetype = "blank", legend.title="Group")
fviz_eig(comb_pca, addlabels = TRUE)

TuBu.info = neuprint_search(".*TuBu.*")
m =neuprint_get_adjacency_matrix(inputids = MeTu.info$bodyid, outputids = TuBu.info$bodyid)
colnames(m) = TuBu.info$name
rownames(m) = MeTu.info$bodyid
meeetuu = m[rownames(m) %in% group$bodyid,]