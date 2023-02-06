#Function for transferring labels (cell type annotations) from one cell data set to another, by cell name
#Provided that the cells are the same or are a subset of one another
#Provide both cds and filename to write out plots
#Return newly annotated cell data set


#ToDO: Check that the annotation is actually present
#Something is weird with plotting device... need another dev.off()?
library(monocle3)
library(dplyr)
library(tibble)
library(assertthat)
library(ggplot2)

transferAnnotations <- function(to_annotate_cds, annotated_cds, filename=NULL, by_cluster=FALSE, annotation_name){
  message(paste0(dim(annotated_cds)[2], " cells in the annotated data set"))
  message(paste0(dim(to_annotate_cds)[2], " cells in the unannotated data set"))
  message(paste0(length(dplyr::intersect(colnames(annotated_cds), colnames(to_annotate_cds))), 
               " cells shared between the two sets"))
  
  #Check if that feature already exists in the dataset to annotate, if so append a unique ID so it doesn't overwrite
  if(annotation_name %in% colnames(pData(to_annotate_cds))){
    warning("Column already exists in cell data set... appending \"transferred\" to annotation_name")
    new_annotation_name <- paste0(annotation_name, "transferred")
  }else{
    new_annotation_name <- annotation_name
  }

  #rownames_to_columns used for easier call to merge()    
  df_pData <- as.data.frame(pData(to_annotate_cds)) %>% rownames_to_column()
  df_annotations <- as.data.frame(pData(annotated_cds)[,c(annotation_name)])
  df_annotations$rowname <- colnames(annotated_cds)
  colnames(df_annotations) <- c(new_annotation_name,"rowname")
  
  #Do the merge on the column titled "rowname", keep all of the cells that we are newly annotating
  #Cells without a match will have "NA" for cell_type after annotation
  #left_join maintains the order of x
  tmp_pData <-dplyr::left_join(x = df_pData, y = df_annotations, by = "rowname")
  
  
  pData(to_annotate_cds)[,new_annotation_name] <- tmp_pData[,annotation_name]
  #Transfer of annotation complete
  
  if(by_cluster){
    print("Generating plot...")
    p<-ggplot(as.data.frame(pData(to_annotate_cds))) +
      geom_bar(aes(x=cluster,fill=cell_type),position="fill") + scale_color_brewer() +
      monocle3:::monocle_theme_opts() + ggtitle("sample representation by cluster")
    if(!is.null(filename)){
      pdf(filename, width = 10)
      print(p)
      dev.off()
    
    } else{
      print(p)
    }  
    
    #What is the most common annotation for each cluster
    grp <- c("cluster",annotation_name)
    
    
    consensus_annotation <- as.data.frame(pData(to_annotate_cds)) %>% group_by_at(grp) %>%  
      summarise(count = n())
    consensus_annotation <- na.omit(consensus_annotation) #removes rows with NA
    consensus_annotation <- consensus_annotation %>% group_by(cluster) %>% top_n(1, count)
    label <- paste0("cluster_consensus_",annotation_name)
    
    #drop counts, which was how many cells in the newly annotated set are the majority cell type for each cluster
    consensus_annotation <- consensus_annotation[c("cluster","cell_type")]
    colnames(consensus_annotation) <- c("cluster",label)
    
    #Prepare df to merge 
    df_pData <- as.data.frame(pData(to_annotate_cds))
    df_annotations <- as.data.frame(consensus_annotation)
    
    
    #Do the merge on the column titled "cluster", 
    #left_join maintains the order of x
    tmp_pData <-dplyr::left_join(x = df_pData, y = df_annotations, by = "cluster")
    
    
    pData(to_annotate_cds)[,label] <- tmp_pData[,label]
    #Transfer of annotation complete
  }
  
  return(to_annotate_cds)
    
}

#Example of how to run
# lmmp <- transferAnnotations(to_annotate_cds = lmmp, 
#                             annotated_cds = MENS,
#                             by_cluster = F,
#                             filename = "/home/jared/ENS/Timecourse_ENS/plots/transferred_cell_types.pdf",
#                             annotation_name = "MENS_subcluster")



