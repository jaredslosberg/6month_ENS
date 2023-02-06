

#Script shared from Charles Zheng 10/11/2020
#Classify cells into cell cycle stages with better resolution than Cyclone
#Gene sets to use from 
library(devtools)
library(dplyr)
library(stringr)
#devtools::install_github('danielschw188/Revelio')


assignCC <- function(log.s.m, batch, cycleGene.l = NULL,  corThres = 0.2, tolerance = 0.1) {
  if (is.null(cycleGene.l)) cycleGene.l <- lapply(as.list(Revelio::revelioTestData_cyclicGenes), function(x) na.omit(as.character(x)))
  # if (is.null(batch)) batch <- rep(1, ncol(log.s.m))
  cycleGene.l <- lapply(cycleGene.l, intersect, rownames(log.s.m))
  allgenes.v <- purrr::reduce(cycleGene.l, union)
  if (sum(allgenes.v %in% rownames(log.s.m)) < 30)  stop("No enough cell cycle genes.")## threshold?
  data.m <- log.s.m[allgenes.v, ]
  
  cc.v <- rep(NA_character_, ncol(data.m))
  for (b in unique(batch)) {
    idx <- which(batch == b)
    dat <- data.m[, idx]
    score.m <- do.call(cbind, lapply(seq_along(cycleGene.l), function(i) {
      gene <- cycleGene.l[[i]]
      mean.v <- colMeans(dat[gene, ]) 
      cor.v <- cor(as.matrix(t(dat[gene, ])), mean.v)
      cor.v[is.na(cor.v)] <- 0
      if (sum(cor.v > 0.2) < 3) warning(str_c("Batch ", b, " phase ", names(cycleGene.l)[i], " too few genes (<3)"))
      message(stringr::str_c("Batch ", b,  " phase ", names(cycleGene.l)[i], " gene:", sum(cor.v > 0.2), "\n"))
      return(colMeans(dat[gene[cor.v > 0.2], ]))
    })) %>% scale() %>% t() %>% scale() %>% t() 
    cc.v[idx] <- apply(score.m, MARGIN = 1, function(s) {
      o.v <- order(s, decreasing = TRUE)
      eta <- o.v[1]
      eta_ <- o.v[2]
      if (((abs(eta - eta_) > 1) & (abs(eta - eta_) < (length(cycleGene.l) - 1))) | ((s[eta] - s[eta_]) < tolerance)) return(NA_character_)
      return(names(cycleGene.l)[eta])
    })
    
  }
  cc.v <- factor(cc.v, levels = names(cycleGene.l))
  return(cc.v)
}

#Comments on initial run:
#Returns about half NA, because for many cells the two most "likely" (?) cycle stage are not consecutive, 
#which decreases confidence in the most likely assigned stage being correct.
#Seems to align well with cyclone for IDing populations that cyclone assigned as G2/M
#but this method looks very heterogenous on UMAP, phases all over the place

###Example of how to run
if(0){
  MENS <- readRDS("/home/jared/ENS/Timecourse_ENS/TC_MENS.rds")
  #For some reason there are duplicate gene short names, which don't work for rows.
  #Not a good solution, but drop the second occurance of the gene name if there is one
  dup.ind <- which(duplicated(fData(MENS[fData(MENS)$num_cells_expressed >= 0])$gene_short_name))
  exprs_mat <- log10(exprs(MENS[-dup.ind,])+1)
  rownames(exprs_mat) <- tolower(fData(MENS[-dup.ind,])$gene_short_name)
  
  #Retrieve gene names, subset to ones that are in my data set
  cycleGene.raw <- lapply(as.list(Revelio::revelioTestData_cyclicGenes), function(x) na.omit(as.character(x)))
  cycleGene.filt <- lapply(cycleGene.raw, function(gene_set){
    attr(gene_set ,"na.action") <- NULL
    match.ind <- tolower(gene_set) %in% rownames(exprs_mat)
    print("Genes not matched being removed:")
    print(gene_set[!match.ind])
    gene_set <- tolower(gene_set[match.ind])
    return(gene_set)})
  
  #If we have batches to run on, specify here
  sample.label <- rep(1,dim(MENS)[2])
  
  cc.label <- assignCC(log.s.m = as.matrix(exprs_mat), batch = sample.label, cycleGene.l = cycleGene.filt)                  
}