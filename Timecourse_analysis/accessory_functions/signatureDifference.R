###Function to compare expression signatures

#Reference signature will be gene vector with avg expression

#Query signatures will be a gene by sample matrix (r x c) 
#Sample could be age, binned pseudotime...

#Query and reference ids should be same format

#center and scaling occurs within genes

require(lsa)

signatureDifference <- function(reference_signature, query_signatures, metric = "cosine", center_scale_reference = F, center_scale_query = F){
  

  ref_names <- rownames(reference_signature)
  query_names <- rownames(query_signatures)
  
  shared_ids <- intersect(ref_names, query_names)
  message(paste0(length(shared_ids)," features shared between reference and query"))
  
  #filtered on shared genes
  ref_sig <- reference_signature[shared_ids,]
  query_sig <- query_signatures[shared_ids,]
  
  #after filtering, center and scale
  if(center_scale_reference){
    ref_sig <- scale(ref_sig) %>% .[,1]
  }
  
  if(center_scale_query){
    query_sig <- scale(query_sig)
  }
  
  
  #TODO: distance of the average, or average of the distances??
  ##for each gene
  distance <- apply(query_sig, MARGIN = 2, FUN = function(x){
    
    if(metric == "euclidian"){ dis <- stats::dist(rbind(x, ref_sig), method = "euclidian")}
    if(metric == "cosine"){ dis <- lsa::cosine(x, ref_sig)}
    
    return(dis)
  })
    
  
  return(distance) 
    
}


# ##test
# samp <- c(1:10)
# names(samp) <- paste0("sample", samp)
# samps <- purrr::map_df(samp, function(x){ rnorm(100, mean = x)})
# 
# samp_ref <- rnorm(100,5)
# names(samp_ref) <- 1:100
# 
# signatureDifference(samp_ref, samps, metric = "euclidian")