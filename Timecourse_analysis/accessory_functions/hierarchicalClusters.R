#Function to produce hierarchical clustering dendrogram from clusters (or other groupings) of cells in scRNA
#cds should be monocle object
#group_by can be be any categorical value in pData on which to average values (provide as string)
#feature should be the name of the feature to generate distances from ("expression", "cell pattern weights")
#If the feature is not expression, provide the feature matrix with rownames being a perfect match to cells in cds
library(monocle3)
library(assertthat)
library(dplyr)

#Pass a matrix, a categorical over which to average, and the mapping between matrix rownames and the categories
#Mappings should be a named vector where the value/level is the group
#TODO: Right now this only works with mappings being a numeric factor because it has to be named
#TODO: transpose matrix so that cells are on columns, as would be for expression matrix
avgValueByCluster <- function(mat, avg_over, mappings){
  assertthat::assert_that(sum(rownames(mat) == names(mappings)) == dim(mat)[1],
                          msg = "Cannot average over groups, cell names are not identical")
  #Since the orders are the same, bind them so that matrix can be subset based on group value
  mat <- cbind(mat, mappings)
  groups <- unique(as.vector(mappings))
  print("Averaging values over groups provided...")
  
  print(head(mat))
  avg_by_group <- NULL
  for(i in 1:length(groups)){
    print(groups[i])
    
    avg <- colMeans(mat %>% filter(mappings == groups[i]) %>% dplyr::select(-mappings))
    avg <- as.data.frame(t(avg), row.names = paste0(groups[i]))
    
    avg_by_group <- rbind(avg_by_group, avg)
  }
  return(avg_by_group)
}

hierarchicalClusters <- function(cds_obj, group_by, feature, featureMat = NULL, min_count = 0){
  print(paste0("Performing hierarchical clustering on ",dim(cds_obj)[2], " cells"))
  print(paste0("Cells split into ",length(unique(pData(cds_obj)[,group_by]))," groups by ",group_by))      
  
  if(!is.null(featureMat)){
    assertthat::assert_that(sum(rownames(featureMat) == colnames(cds_obj)) == dim(cds_obj)[2],
                            msg = "Rownames of the feature matrix are not the same as the cells in the cell data set")
  }
  if(feature == "expression"){
    #Featuremat should be expression values, passed with rows as cells and columns as genes (ie transposed)
    cds_obj <- detect_genes(cds_obj)
    featureMat <- as.data.frame(t(as.matrix(exprs(cds_obj[fData(cds_obj)$num_cells_expressed > min_count]))))
    print(paste0(dim(cds_obj)[1], " genes reduced to ",
                 sum(fData(cds_obj)$num_cells_expressed > min_count)," genes with > ",
                 min_count, " counts"))
  }
  #TODO: generalize this to other annotations
  if(group_by == "cell_type"){
    #This line needs to be manually changed to accomodate character factors, make sure that they are named
    pData(cds_obj)$group_id <- as.factor(paste0(pData(cds_obj)[,"cluster"], " - ", pData(cds_obj)[,"cell_type"]))
    names(pData(cds_obj)$group_id) <- colnames(cds_obj)
    group_by <- "group_id"
  }
  else{
    print(paste0("grouping by", group_by))
  }
     #returns a data.frame with each row being a group, and each column being the average feature value for the group
  #not necessarily in the right order
  
  by_feature <- avgValueByCluster(featureMat,avg_over = group_by, mappings = pData(cds_obj)[,group_by])
  
  euc_distances<- dist(by_feature)
  dendrogram <- hclust(euc_distances)
  return(dendrogram)
}

###Example of how to run
if (0){
  cell_weights <- read.csv("/home/jared/ENS/Timecourse_ENS/results/NMF/MENS/MENS_A_pattern_cell_weights.csv", row.names = 1)
  
  tree <- hierarchicalClusters(cds_obj = MENS,
                              group_by = "clusters",
                              feature = "NMF",
                              featureMat = cell_weights,
                              #min_count = 10
                              )
}
