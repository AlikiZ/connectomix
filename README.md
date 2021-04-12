# connectomix

The goal is to cluster some neurons of the organism Drosophila utilizing the publicly available hemibrain dataset after data intergration and processing. The visualization of our result after the categorization, where one colour corresponds to one cateogry, should look like this (an example of categorization on other neurons):

![TuBu_neurons_in_AOTU_dorsal_top_view](https://user-images.githubusercontent.com/38488804/114377206-ae78b000-9b86-11eb-9e09-cb7e94c917ad.png)



This repository hosts a snippet of the currently ongoing work, as we are constantly expanding our analysis. Hopefully soon, we can update the repository with more results and reveal more of our ideas.
The script tocluster.R contains the data processing and clustering process using PCA analysis, kmeans clustering and hierarchical clustering both on connectivity and morphology matrices.

The script PCA_extended.R uses the MeTu_grouping.csv file provided from Emil Kind to perform PCA on an extended MeTu - TuBu matrix. 

The plots in the folder connectivity_ext are created with the script PCA_extended.R.
