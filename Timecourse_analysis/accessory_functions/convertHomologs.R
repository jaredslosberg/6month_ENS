#' Converts gene names of query to match type/species of reference names (human
#' or mouse).
#'
#'Modified from Azimuth::ConvertGeneNames
#'
library(stringr)

convertHomologs <- function(gene.names, reference.names, homolog.table) {
  query.names <- gene.names
  linked <- homolog.table
  
  # remove version numbers
  ensembl.id <- '(?:ENSG|ENSMUS)'
  capture.nonversion <- paste0('(', ensembl.id, '.*)\\.(?:.*)')
  pattern <- paste0('(?:', capture.nonversion,')|(.*)')
  query.names <- Filter(
    f = function(x) !is.na(x = x),
    x = as.vector(x = t(x = str_match(string = query.names, pattern = pattern)[, 2:3]))
  )
  # determine idtype and species of query
  # 5000 because sometimes ensembl IDs are mixed with regular gene names
  query.names.sub = sample(
    x = query.names,
    size = min(length(x = query.names), 5000)
  )
  idtype <- names(x = which.max(x = apply(
    X = linked,
    MARGIN = 2,
    FUN = function(col) {
      length(x = intersect(x = col, y = query.names.sub))
    }
  )))
  species <- ifelse(test = length(x = grep(pattern = '\\.mouse', x = idtype)) > 0, 'mouse', 'human')
  idtype <- gsub(pattern = '\\.mouse|\\.human', replacement = '', x = idtype)
  message('detected inputs from ', toupper(x = species), ' with id type ', idtype)
  totalid <- paste0(idtype, '.', species)
  
  # determine idtype and species of ref
  reference.names.sub <- sample(x = reference.names, size = min(length(x = reference.names), 5000))
  idtype.ref <- names(x = which.max(x = apply(
    X = linked,
    MARGIN = 2,
    FUN = function(col) {
      length(x = intersect(x = col, y = reference.names))
    }
  )))
  species.ref <- ifelse(test = length(x = grep(pattern = '\\.mouse', x = idtype.ref)) > 0, 'mouse', 'human')
  idtype.ref <- gsub(pattern = '\\.mouse|\\.human', replacement = '', x = idtype.ref)
  message('reference rownames detected ', toupper(x = species.ref),' with id type ', idtype.ref)
  totalid.ref <- paste0(idtype.ref, '.', species.ref)
  
  if (totalid == totalid.ref) {
    return(geneNames)
  } else {
    display.names <- setNames(
      c('ENSEMBL gene', 'ENSEMBL transcript', 'gene name', 'transcript name'),
      nm = c('Gene.stable.ID', 'Transcript.stable.ID','Gene.name','Transcript.name')
    )

    # set up table indexed by query ids (totalid)
    linked.unique <- linked[!duplicated(x = linked[[totalid]]), ]
    new.indices <- which(query.names %in% linked.unique[[totalid]])
    message("Found ", length(x = new.indices), " out of ",
            length(x = query.names), " total inputs in conversion table")
    query.names <- query.names[new.indices]
    rownames(x = linked.unique) <- linked.unique[[totalid]]
    linked.unique <- linked.unique[query.names, ]
    # get converted totalid.ref
    new.names <- linked.unique[, totalid.ref]
    # remove duplicates
    notdup <- !duplicated(x = new.names)
    new.indices <- new.indices[notdup]
    new.names <- new.names[notdup]
    
    new_name_vec <- rep(NA_character_, length(gene.names))
    new_name_vec[new.indices] <- new.names
    matched_homologs <- data.frame(input = gene.names, converted = new_name_vec)
    return(matched_homologs)
  }
}