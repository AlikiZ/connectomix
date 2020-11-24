library(neuprintr)
library(nat)
library(ggpubr)


## connection to neuprint API for dataset = "hemibrain:v1.1" by calling another script so as not to reveal token
source("/home/alikiz/Documents/connectomics/code/loadhemibrain.R")

## getting some info for our cells of interest (MC61 & MC64 aka “MeTu") including the bodyIDs
MC61.info = neuprint_search("MC61") 
MC64.info = neuprint_search("MC64")
MC.info = rbind(MC61.info,MC64.info)

## getting some info for “TuBu" cells (main downstream target of MC61 & MC64) including the bodyIDs 
TuBu.info = neuprint_search(".*TuBu.*")
#TuBu.info = TuBu.info[complete.cases(TuBu.info),]

## connectivity matrix for "MeTu > TuBu” connections 
matrix = neuprint_get_adjacency_matrix(inputids = MC.info$bodyid, outputids = TuBu.info$bodyid)
colnames(matrix) = TuBu.info$name
rownames(matrix) = MC.info$bodyid

## Looking at some basic morphological parameters for one example cell (‘1108796315’)
test = neuprint_read_neurons(1108796315)
summary(test)

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
## check which columns contain nonzero entries all over the column, if a column contains only zero entries then remove it assuming it does not contirbute to the solution
nonzerocol = apply(DRAMeTu_connect, 2, function(x) !all(x==0))
connect_mat = DRAMeTu_connect[,nonzerocol]
# dim(relevant_mat[,this]) should be 53 * 25
heatmap(connect_mat)

# connect_DRAMeTu_sub1 = matrix[rownames(matrix) %in% DRAMeTu_sub1_ID,]
# connect_DRAMeTu_sub2 = matrix[rownames(matrix) %in% DRAMeTu_sub2_ID,]
# heatmap(connect_DRAMeTu_sub1)

# cluster in two groups , possibly TO DO summary statistics by group
true_classes = rownames(connect_mat) %in% DRAMeTu_sub1_ID + 1

# K-Means Cluster Analysis, the initial pick of the center influences the output so we change the initial canter by setting nstart = various_values
fit <- kmeans(x = connect_mat, centers = 2, nstart = 1)
fit2 <- kmeans(x = connect_mat, centers = 2, nstart = 10)

# explore results of k-means e.g. size of groups
fit$size
fit2$size

# append cluster assignment at the end of matrix mydata
mydata <- data.frame(connect_mat, true_classes, fit$cluster, fit2$cluster)
# calculate the match between real and predicted group
sum((mydata$true_classes - mydata$fit.cluster)^2)
