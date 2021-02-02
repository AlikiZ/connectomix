# connectomix

The goal is to cluster the neurons of hemibrain after data intergration and processing

The script tocluster.R is the initial data processing and clustering attempt containg PCA analysis, kmeans clustering and hierarchical clustering both on connectivity and morphology matrices.

The script PCA_extended.R uses the MeTu_grouping.csv file from Emil to perform PCA on an extended MeTu - TuBu matrix. The goal is to try PCA on a sequence of MeTu - TuBu matrices: all of the MeTus listed in the MeTu_grouping.csv, with focus on each of the 4 subgroups of the MeTus.

The folder connectivity_sub1 contains a descriptive txt file and plots created with the script tocluster.R when named with the file extension smallcol.png and plots created with the script PCA_extended.R when named with the file extension bigcol.png. The files in connectivity_ext are created with the script PCA_extended.R.

