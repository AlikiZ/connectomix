library(neuprintr)
library(nat)
library(ggpubr)
library(ggfortify)
library(factoextra)
library(FactoMineR)
library(corrplot)


### small functions ###
confmat <- function(actual, predicted){
  confusionmat = as.matrix(table(Actual= true_classes, Predicted=fit$cluster))
  print(confusionmat)
  accurracy = sum(diag(cm)) / sum(cm)
  metr = list("confusionmat" = cm, "accuracy" = accurracy)
  return(metr)
}


#### data loading and preprocessing ####

## connection to neuprint API for dataset = "hemibrain:v1.1" by calling another script so as not to reveal token
source("/home/aliki/WernetLab/code/loadhemibrain.R")

## getting some info for our cells of interest (MC61 & MC64 aka “MeTu") including the bodyIDs
MC61.info = neuprint_search("MC61") 
MC64.info = neuprint_search("MC64")
MC.info = rbind(MC61.info,MC64.info)

## getting some info for “TuBu" cells (main downstream target of MC61 & MC64) including the bodyIDs 
TuBu.info = neuprint_search(".*TuBu.*")
#TuBu.info = TuBu.info[complete.cases(TuBu.info),]

## connectivity matrix for "MeTu > TuBu” connections 
#neuprint_connection_table() 
#inputids = MC.info$bodyid, outputids = pre or postsynaptic cells
matrix = neuprint_get_adjacency_matrix(inputids = MC.info$bodyid, outputids = TuBu.info$bodyid)
colnames(matrix) = TuBu.info$name
rownames(matrix) = MC.info$bodyid

## bodyIDs from Aljoscha's cluster analysis for the two DRAMeTu subgroups: sub1 has 14 elements, sub2 has 39 elements
DRAMeTu_sub1_ID = c(1078097609, 1201227161, 1171548392, 1076729226, 1202250779, 1077756950, 5813016373, 1198127720, 1078102468,
                    1108796315, 1108446169, 5812987715, 1077761540, 1170179426)

DRAMeTu_sub2_ID = c(1139826933, 1077756731, 1078097817, 1198809632, 1076729233, 1078102490, 1107423158, 1078102267, 1139485872,
                    1077761532, 1078102350, 1139148877, 1140172133, 1170861854, 1201559521, 1077424293, 1201572715, 1076729262,
                    1201210029, 5901210109, 1199155111, 1077424204, 1198809646, 1171211215, 5813104106, 1078102400, 1077424267,
                    1139831143, 1139831106, 1200873417, 1232939935, 1140517287, 1108450770, 5812986562, 1078102419, 1169834028,
                    1202250762, 1201905111, 1046039977)

## get the matrix for the DRAMeTu rows: 14 + 39 = 53 rows
DRAMeTu_connect = matrix[rownames(matrix) %in% c(DRAMeTu_sub1_ID,DRAMeTu_sub2_ID),]
## check which columns contain nonzero entries all over the column, if a column contains only zero entries then remove it assuming it does not contribute to the solution
nonzerocol = apply(DRAMeTu_connect, 2, function(x) !all(x==0))
connect_mat = DRAMeTu_connect[,nonzerocol]
# dim(relevant_mat[,this]) should be 53 * 25
true_classes = rownames(connect_mat) %in% DRAMeTu_sub1_ID + 1
names(true_classes) = rownames(connect_mat)
heatmap(connect_mat)
# normalize by dividing by the column sum
scaled = scale(connect_mat, center=FALSE, scale=colSums(connect_mat))


#### algorithms ####


# K-Means Cluster Analysis, the initial pick of the center influences the output so we change the initial canter by setting nstart = various_values
fit <- kmeans(x = connect_mat, centers = 2, nstart = 1)
fit2 <- kmeans(x = connect_mat, centers = 2, nstart = 10)
# explore results of k-means e.g. size of groups
fit$size
fit2$size

# calculate evaluation measures assuming there are true labels (correct classes of neurons). having true labels is not typical in unsupervised models like clustering
# the labels come from Aljoscha's cluster analysis. if rowname of neuron is in the DRAMeTu_sub1 then next line evaluates to TRUE (equivalent to 1) + 1 = 2; (possibly TO DO summary statistics by group)
true_classes = rownames(connect_mat) %in% DRAMeTu_sub1_ID + 1
names(true_classes) = rownames(connect_mat)
# append cluster assignment at the end of matrix mydata
mydata <- data.frame(connect_mat, true_classes, fit$cluster, fit2$cluster)
mydata$label <- as.character(mydata$true_classes)
# One can compare the true labels with the predicted labels in terms of some metrics: confusion matrix (cm), accuracy (acc)
kmeans1 = confmat(true_classes, fit$cluster)

# Compute dissimilarity matrix, equivalently the distance matrix ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
res.dist <- dist(connect_mat, method = "euclidean")
# Compute hierarchical clustering with various methods ("ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC) )
res.hc <- hclust(res.dist, method = "ward.D2")
plot(res.hc, cex = 0.5)

# set a cut off by count of the desired number of clusters https://www.r-bloggers.com/2016/01/hierarchical-clustering-in-r-2/
clusterCut <- cutree(res.hc, 2)
cm2 = as.matrix(table(Actual= true_classes, Predicted=clusterCut))
cm2 

#### PCA analysis ####
# PCA - not scaled and scaled --> the second plot of not scaled data reveals a better result
pca_res <- prcomp(connect_mat, scale. = TRUE)
autoplot(pca_res, data = mydata, colour = "true_classes")
somepca <- PCA(connect_mat, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= mydata$label, palette = c("#00AFBB", "#E7B800"), axes.linetype = "blank", legend.title="Group")

pca_res <- prcomp(connect_mat, scale. = FALSE)
autoplot(pca_res, data = mydata, colour = "true_classes")
somepca <- PCA(connect_mat, scale.unit = FALSE, graph = FALSE)
fviz_pca_ind(somepca, geom.ind = "point", col.ind= mydata$label, palette = c("orange", "green3"), axes.linetype = "blank", legend.title="Group")
fviz_pca_ind(somepca, geom.ind = "point", col.ind= mydata$label, palette = c("orange", "green3"), axes.linetype = "blank", legend.title="Group", axes = c(3,4))
fviz_eig(somepca, addlabels = TRUE)
# fviz_pca_var(res.pca, col.var = "black")  # for correlation circle of the variables, important: variables need deiffrent names
# more: contribution of variables to PCi, variables according to their contributions to the principal components (dimension description)

# new matrix: comb, combination of both 
comb = cbind(connect_mat, mat_both)
comb_label = cbind(mydata, mat_both)
comb_pca = PCA(new, scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(comb_pca, geom.ind = "point", col.ind= new_groups$label, palette = c("#00AFBB", "#E7B800"), axes.linetype = "blank", legend.title="Group")
fviz_eig(comb_pca, addlabels = TRUE)


## get the summary for more neurons to check if I capture more variability so that it is easier to see the clusters
## build a matrix of the neuron information (morphology)
aotu <- neuprint_ROI_mesh("AOTU(R)")
aotu.mesh <- as.hxsurf(aotu)

a = data.frame()
for (i in 1:length(DRAMeTu_sub1_ID)) {
  someneuron <- neuprint_read_neurons(DRAMeTu_sub1_ID[i])
  someneuron_pruned <- prune_in_volume(someneuron, surf=aotu.mesh )
  info <- as.vector(summary(someneuron_pruned))
  a <- rbind(a, info)
}  

b = data.frame()
for (i in 1:length(DRAMeTu_sub2_ID)) {
  someneuron <- neuprint_read_neurons(DRAMeTu_sub2_ID[i])
  someneuron_pruned <- prune_in_volume(someneuron, surf=aotu.mesh )
  info <- as.vector(summary(someneuron_pruned))
  b <- rbind(b, info)
}  

both_summary <- rbind(a, b)
mat_both <- as.matrix(both_summary)
heatmap(mat_both[,-1], scale = "col", Rowv = NA, Colv = NA)
heatmap(mat_both[,-c(1,6)], scale = "none", Rowv = NA)

# PCA - not scaled and scaled --> when scaled PC1 has 100% information
morphpca <- PCA(mat_both[,-1], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(morphpca, geom.ind = "point", col.ind= mydata$label, palette = c("#00AFBB", "#E7B800"), axes.linetype = "blank", legend.title="Group")
var <- get_pca_var(morphpca)
head(var$contrib)
fviz_eig(morphpca, addlabels = TRUE)
fviz_pca_var(morphpca, repel = TRUE)
corrplot(var$contrib, is.corr=FALSE)    

morphpca <- PCA(mat_both[,-1], scale.unit = FALSE, graph = FALSE)
fviz_pca_ind(morphpca, geom.ind = "point", col.ind= mydata$label, palette = c("orange", "green3"), axes.linetype = "blank", legend.title="Group")

pca_res <- prcomp(mat_both, scale. = FALSE)
autoplot(pca_res, data = mydata, colour = "true_classes")

## TO DO: normal distribution is not necessary (check if the distribution is normal per column), plot more PCAs because those do not show much
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#eigenvalues-variances

# discussed on Friday
# DRAMeTu_no_1 <- neuprint_read_neurons("1078097609")
# aotu <- neuprint_ROI_mesh("AOTU(R)")
# aotu.mesh <- as.hxsurf(aotu)
# DRAMETu_pruned <- prune_in_volume(DRAMeTu_no_1,surf = aotu.mesh)
# b<-summary(DRAMETu_pruned) 
