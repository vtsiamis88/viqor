library(igraph)
library(preprocessCore)
library("protr")
library(ggplot2)
library(reshape2)
library(shiny)
library(heatmaply)
library(dplyr)
library(Biostrings)

#Function that applies general filtering and aggregation operations:
#1)Run/Condition arrangment of quantitative columns
#2)Based on the type of input dataset removes the rows without peptide/protein info.
#3)Replaces 0, inf and -inf with NA.
#4)Removes peptides that their missing values exceed a user defined threshold.
#5)Aggregates duplicate peptides along with their adundances and their protein membership (if included).
general.Filtering <- function(peptide.Data.Frame, type = "include", total.quantitative.columns = 0, total.conditions = 0, maxNA = 0, log = F, column.arrangement = F){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #type: either "include" when proteins are assigned to every peptide or "missing" otherwise.
  #total.quantitative.columns: the total number of columns in the data frame that contain peptide abundances.
  #total.conditions: the total number of conditions.
  #maxNA: the maximum number of allowed missing values per peptide.
  #log: indicates whether the data are in log scale.
  #column.arrangement: F if the columns are arranged by condition indices or T if they are arranged by run indices.
  
  if(missing(peptide.Data.Frame)){
    
    stop("Error: Input data frame is missing!")
    
  }

  #1)If the dataset's columns are arranged based on the run indices then arrange it based on condition indices.
  if(column.arrangement){
    
    arrangement <- c((ncol(peptide.Data.Frame) - total.quantitative.columns + 1):ncol(peptide.Data.Frame))
    new.arrangement <- unlist(lapply(c(1:total.conditions), function(x) arrangement[seq(x,length(arrangement),total.conditions)]))
    peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)] <- peptide.Data.Frame[,new.arrangement]
    colnames(peptide.Data.Frame)[((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)] <- colnames(peptide.Data.Frame)[new.arrangement]
    
  }
  
  #2)Remove of rows with missing peptide/protein info (sequence/accession number).
  if(type == "include"){
    
    indices.Empty.Pep <- which(peptide.Data.Frame[,1] == "")
    indices.NA.Pep <- which(is.na(peptide.Data.Frame[,1]))
    indices.Empty.Prot <- which(peptide.Data.Frame[,2] == "")
    indices.NA.Prot <- which(is.na(peptide.Data.Frame[,2]))
    
    indices <- union(union(indices.Empty.Prot,indices.NA.Prot), union(indices.Empty.Pep,indices.NA.Pep))
    
    if(length(indices) > 0){
      
      peptide.Data.Frame <- peptide.Data.Frame[-indices,]
      rownames(peptide.Data.Frame) <- c(1:nrow(peptide.Data.Frame))
      
    }
    
  } else if(type == "missing"){
    
    indices.Empty.Pep <- which(peptide.Data.Frame[,1] == "")
    indices.NA.Pep <- which(is.na(peptide.Data.Frame[,1]))
    
    indices <- union(indices.Empty.Pep, indices.NA.Pep)
    
    if(length(indices) > 0){
      
      peptide.Data.Frame <- peptide.Data.Frame[-indices,]
      rownames(peptide.Data.Frame) <- c(1:nrow(peptide.Data.Frame))
      
    }
    
  } else {
    
    stop("Error: Not a valid type of dataset!")
    
  }
  
  #3)Replace 0 -Inf and Inf with NA.
  peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)][peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)] == 0] <- NA
  peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)][peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)] == -Inf] <- NA
  peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)][peptide.Data.Frame[,((ncol(peptide.Data.Frame) - total.quantitative.columns) + 1):ncol(peptide.Data.Frame)] == Inf] <- NA
  
  #4)Aggregate duplicate peptides along with their abundances and protein membership (if included).
  if(length(unique(peptide.Data.Frame[,1])) < nrow(peptide.Data.Frame)){
    
    if(log){
      
      classes <- unlist(sapply(1:ncol(peptide.Data.Frame), function(x) class(peptide.Data.Frame[,x])))
      to.unlog.columns <- intersect(which(classes == "numeric" | classes == "integer"), (ncol(peptide.Data.Frame) - total.quantitative.columns + 1):ncol(peptide.Data.Frame))
      peptide.Data.Frame[, to.unlog.columns] <- 2^peptide.Data.Frame[, to.unlog.columns]
      
    }
    
    if(type == "missing"){
      
      col.names <- colnames(peptide.Data.Frame)[1]
      
      classes <- unlist(sapply(2:ncol(peptide.Data.Frame), function(x) class(peptide.Data.Frame[,x])))
      numeric <- which(classes == "numeric" | classes == "integer") + 1
      character <- which(classes == "character") + 1 
      data.numeric <- aggregate(x = peptide.Data.Frame[, numeric]/1, by = list(peptide.Data.Frame[, 1]), FUN = sum, na.rm = T)
      
      if(length(character) > 0){
        
        data.character <- aggregate(x = peptide.Data.Frame[, character], by = list(peptide.Data.Frame[, 1]), FUN = paste, collapse = "/")
        colnames(data.character) <- c(col.names, colnames(peptide.Data.Frame)[character])
        peptide.Data.Frame <- cbind(data.character, data.numeric[, 2:ncol(data.numeric)])
        
      } else {
        
        peptide.Data.Frame <- data.numeric
        
      }
      
      colnames(peptide.Data.Frame)[1] <- col.names
      
      peptide.Data.Frame[peptide.Data.Frame == 0] <- NA
      
    } else if(type == "include"){ #Not used anymore. It should work properly when 1st and 2nd columns correspond to peps and proteins. 
      
      col.names <- colnames(peptide.Data.Frame)[1:2]
      peptide.Data.Frame <- aggregate(x = peptide.Data.Frame[,3:ncol(peptide.Data.Frame)]/1, by = list(peptide.Data.Frame[,1], peptide.Data.Frame[,2]), FUN = sum, na.rm = T)
      
      protein.Aggregation <- aggregate(x = peptide.Data.Frame[,2], by = list(peptide.Data.Frame[,1]), FUN = paste, collapse = "/")
      abundances.Aggregation <- aggregate(x = peptide.Data.Frame[,3:ncol(peptide.Data.Frame)]/1, by = list(peptide.Data.Frame[,1]), FUN = sum, na.rm = T)
      
      peptide.Data.Frame <- data.frame(protein.Aggregation, abundances.Aggregation[,2:ncol(abundances.Aggregation)])
      
      colnames(peptide.Data.Frame)[1:2] <- col.names
      
      peptide.Data.Frame[peptide.Data.Frame == 0] <- NA
      
    }
    
    if(log){
      
      classes <- unlist(sapply(1:ncol(peptide.Data.Frame), function(x) class(peptide.Data.Frame[,x])))
      to.log.columns <- intersect(which(classes == "numeric" | classes == "integer"), (ncol(peptide.Data.Frame) - total.quantitative.columns + 1):ncol(peptide.Data.Frame))
      peptide.Data.Frame[, to.log.columns] <- log2(peptide.Data.Frame[, to.log.columns])
      
    }
    
  }
  
  #5)Removing peptides that have more NAs than maxNA.
  indices <- which(apply(peptide.Data.Frame, 1, function(x) sum(is.na(x))) > maxNA)
  
  if(length(indices) > 0){
    
    peptide.Data.Frame <- peptide.Data.Frame[-indices,]
    row.names(peptide.Data.Frame) <- c(1:nrow(peptide.Data.Frame))
    
  }
  
  #6) Remove peptides with X symbol
  unknown <- which(grepl("X", peptide.Data.Frame[,1]))
  
  if(length(unknown) > 0){
    
    peptide.Data.Frame <- peptide.Data.Frame[-unknown,]
    rownames(peptide.Data.Frame) <- 1:nrow(peptide.Data.Frame)
    
  }
  
  return(peptide.Data.Frame)
}

#Function that constructs a matrix for the conditions and their samples based of the experimental setup.
samples.And.Groups <- function(peptide.Data.Frame, total.conditions, total.quantitative.columns, reference.condition){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #total.conditions: the total number of conditions in the dataset (including reference condition)
  #total.quantitative.columns: the total number of columns in the data frame that contain peptide abundances.
  #reference.condition: the index of the reference condition.
  #labels: a matrix in which the columns and rows correspond to the name of conditions and quantitative columns identifiers.
  
  if(missing(total.conditions)){
    
    total.conditions <- 1
    
  }
  
  if(missing(total.quantitative.columns)){
    
    stop("Error: total quantitative columns attribute is missing!")
    
  }
  
  if(missing(reference.condition)){
    
    reference.condition <- 0
    
  }
  
  quantitative.columns.names <- colnames(peptide.Data.Frame[,(ncol(peptide.Data.Frame) - total.quantitative.columns + 1):ncol(peptide.Data.Frame)])
  
  labels <- matrix(ncol = total.conditions, nrow = total.quantitative.columns/total.conditions)
  colnames(labels) <- c(1:total.conditions)
  
  condition.counter <- 1
  
  for(i in 1:total.conditions){
    
    if(i != reference.condition){
      
      labels[,i] <- quantitative.columns.names[1:(total.quantitative.columns/total.conditions)]
      quantitative.columns.names <- quantitative.columns.names[-c(1:(total.quantitative.columns/total.conditions))]
      colnames(labels)[i] <- paste("Condition_", condition.counter, sep = "")
      condition.counter <- condition.counter + 1
      
    } else {
      
      labels[,i] <- quantitative.columns.names[1:(total.quantitative.columns/total.conditions)]
      quantitative.columns.names <- quantitative.columns.names[-c(1:(total.quantitative.columns/total.conditions))]
      colnames(labels)[i] <-  "Reference"
      
    }
  }
  
  return(labels)
  
}

#Function that finds the types of modification in the dataset. Recognizes 3 different modification expression types.
PTM.recognition <- function(peptide.Data.Frame){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #peptide.sequences: the peptide sequences.
  #model: 0 -> no modifications found, 1 -> "()", 2 -> "[]", 3 -> [a-z].
  #recognized.modifications: list that contains the types of modifications found and the model they were expressed.
  
  if(missing(peptide.Data.Frame)){
    
    stop("Error! Peptide dataframe is missing!")
    
  }
  
  peptide.sequences <- peptide.Data.Frame[,1]
  
  modification.types <- NULL
  model <- 0
  
  peptide.string <- paste(peptide.sequences, collapse = "")
  
  if(grepl("\\(", peptide.string)){
    
    expression <- "\\(.*?\\)"
    substitution <- "\\(|\\)"
    model <- 1
    
  } else if(grepl("\\[", peptide.string)){
    
    expression <- "\\[.*?\\]"
    substitution <- "\\[|\\]"
    model <- 2
    
  } else if(grepl("[a-z]", peptide.string)) {
    
    expression <- "[a-z]*[a-z]"
    substitution <- ""
    model <- 3
    
  } else {
    
    model <- 0
    
  }
  
  if(model != 0){
    
    modification.types <- unique(unlist(lapply(regmatches(peptide.string, gregexpr(expression, peptide.string)), function(x) gsub(substitution, "", x))))
    
  }
  
  recognized.modifications <- list(model = model, modifications = modification.types)
  
  return(recognized.modifications)
  
}

#Function that extracts the selected modification types and their position for each modified peptide, as well as the index of the peptide in the dataset.
PTM.extraction <- function(peptide.Data.Frame, keep.modifications, all.modifications, model){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #keep.modifications: the types of modifications that will procceed in PTMs quantification.
  #all.modifications: all the types of modifications found in the peptide data set.
  #model: 1 -> "()", 2 -> "[]", 3 -> [a-z].
  #sequence: a vector that contains all the modified peptide sequences.
  #counterpart: a vector that contains the counterpart of all modified peptides (counterpart: the unmodified version of that peptide)
  #ptm.type: a list that stores the type(s) of PRTM for every peptide.
  #ptm.position: a list that contains the position of each PTM of the sequence of each corresponding modified peptide.
  #ptm.index: the index of the modified peptide in the dataset.
  #ptm.list: a list that contains the above data structures.
  
  if(missing(peptide.Data.Frame) | missing(keep.modifications) | missing(all.modifications) | missing(model)){
    
    stop("Error: one or more of the input arguments are missing!")
    
  }
  
  peptide.sequences <- peptide.Data.Frame[,1]
  remove.modifications <- setdiff(all.modifications, keep.modifications)
  
  sequence <- vector(mode = "character", length = 0)
  counterpart <- vector(mode = "character", length = 0)
  ptm.type <- list()
  ptm.position <- list()
  ptm.index <- vector(mode = "numeric", length = 0)
  
  if(model == 1){
    
    expression <- "\\(.*?\\)"
    position.char <- "\\("
    substitution <- "\\(|\\)"
    mod.substitution.par <- "\\(\\)"
    
  } else if(model == 2){
    
    expression <- "\\[.*?\\]"
    position.char <- "\\["
    substitution <- "\\[|\\]"
    mod.substitution.par <- "\\[\\]"
    
  } else if(model == 3){
    
    expression <- "[a-z]*[a-z]"
    position.char <- "([a-z])"
    substitution <- ""
    mod.substitution.par <- ""
    
  } else {
    
    stop("Error: wrong input model!")
    
  }
  
  
  if(length(remove.modifications) > 0){
    
    mod.substitution <- paste(remove.modifications, collapse = "|")
    
  } else {
    
    mod.substitution <- ""
    
  }
  
  peptide.sequences <- unlist(lapply(peptide.sequences, function(x) gsub(mod.substitution, "", x)))
  peptide.sequences <- unlist(lapply(peptide.sequences, function(x) gsub(mod.substitution.par, "", x)))
  
  for(i in 1:length(peptide.sequences)){
    
    ptm <- as.character(unlist(regmatches(peptide.sequences[i], gregexpr(expression, peptide.sequences[i]))))
    
    if(length(ptm) > 0){
      
      temp.sequence <- peptide.sequences[i]
      temp.position <- vector(mode = "numeric", length = length(ptm))
      
      for(j in 1:length(ptm)){
        
        temp.position[j] <- unlist(gregexpr(position.char, temp.sequence))[1]
        
        temp.sequence <- sub(paste("\\", ptm[j], sep = "", collapse = ""), "", temp.sequence)
        
      }
      
      sequence[length(sequence) + 1] <- peptide.sequences[i]
      counterpart[length(counterpart) + 1] <- temp.sequence
      ptm.type[[length(ptm.type) + 1]] <- gsub(substitution, "", ptm)
      ptm.position[[length(ptm.position) + 1]] <- temp.position
      ptm.index[length(ptm.index) + 1] <- i
      
    }
  }

  #Modified peptide sequences will remain intact. Leucine substituted sequences will be used only for visualization in VIQoR plot.
  ptm.list <- list(sequences = sequence, sequences.sub = gsub("[I]|[J]", "L", sequence), counterpart = counterpart, counterpart.sub = gsub("[I]|[J]", "L", counterpart),index = ptm.index, type = ptm.type, position = ptm.position)
  
  return(ptm.list)
  
}

#Function that searches modified peptide sequences against a FASTA file specified by user and returns their mapping.
mapping.Modified.Peptides <- function(counterparts, fasta){
  
  #counterparts: a set of modified peptide counterparts.
  #fasta: a fasta file that contains protein sequences. Peptides will be searched against that file.
  #mapping: the final mapping of the counterpart sequences.
  
  if(missing(counterparts) | missing(fasta)){
    
    stop("Error: required input attributes are missing in mapping.Modified.Peptides().")
    
  }
  
  shiny::withProgress(message = "Modified peptide mapping: ", min = 0, max = 3, value = 0, {
  
    counterparts <- sort(unique(counterparts))
  
    temp_counterpart_sequences_file <- tempfile(pattern = "counterpart_sequences", fileext = ".txt")
    write.table(counterparts, temp_counterpart_sequences_file, row.names = F, col.names = F, quote = F)
  
    fasta <-Biostrings::readAAStringSet(fasta, format="fasta", use.names=TRUE)
    fasta <- chartr(old = "I", new = "L", fasta)
    temp_fasta_file <- tempfile(pattern = "fasta_file", fileext = ".txt")
    Biostrings::writeXStringSet(fasta, temp_fasta_file, format = "fasta")
    
    shiny::incProgress(1, detail = "files prepared")
  
    temp_counterparts_mapping_file <- tempfile(pattern = "counterparts_mapping", fileext = ".csv")
    
    cat("\n")
    cat("Modified peptide mapping started:\n")
    cat("\n")
    
    #Separator for lib in class path (Windows)
    lib.separator <- ";"
    
    if(substr(Sys.info()['sysname'], 1, 1) != "W"){
      
      lib.separator <- ":" #(Linux)
      
    }
    
    system(paste("java -cp utilities-4.12.9/utilities-4.12.9.jar", lib.separator, "utilities-4.12.9/* com.compomics.util.experiment.identification.protein_inference.executable.PeptideMapping -p ", temp_fasta_file, " ", temp_counterpart_sequences_file, " ", temp_counterparts_mapping_file, sep = "", collapse = ""))

    cat("Modified peptide mapping finished!\n")
    
    shiny::incProgress(1, detail = "mapping is done")
    
    mapping <- read.csv(temp_counterparts_mapping_file, header=F, sep = ",", stringsAsFactors = FALSE, comment.char = "#")
    colnames(mapping) <- c("Counterpart", "Protein","Position")
    
    fasta <- protr::readFASTA(file = temp_fasta_file, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- matrix(cbind(sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), unlist(fasta)), ncol = 2)
    
    protein.Sequences <- unlist(lapply(mapping[,2], function(x) fasta[which(fasta[,1] == x),2]))
    
    mapping <- data.frame(mapping, protein.Sequences, stringsAsFactors = FALSE)
  
    duplicates <- which(duplicated(mapping[,c(1,2)]))
  
    if(length(duplicates) > 0){
    
      mapping <- mapping[-duplicates,]
    
    }
  
    unknown <- which(grepl("X", mapping$Counterpart))
  
    if(length(unknown) > 0){
    
      mapping <- mapping[-unknown,]
    
    }
  
    rownames(mapping) <- c(1:length(mapping$Counterpart))
  
    file.remove(temp_counterparts_mapping_file, temp_fasta_file, paste(temp_fasta_file, ".cui", sep = "", collapse = ""), temp_counterpart_sequences_file)
  
    shiny::incProgress(1, detail = "post-processing is done")
    
  })
  
  mapping$Position <- mapping$Position + 1
  
  return(mapping)
  
}

#Function that performs per sample log zero center normalization.
log.Zero.Center.Normalization <- function(peptide.Data.Frame, labels, log = FALSE, method = "median"){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #labels: a matrix in which the columns and rows correspond to the name of conditions and quantitative columns identifiers.
  #log: FALSE is peptide.Data.Frame is not log transformed, elsewise TRUE.
  #method: the method which the dataset will  be normalized. (average, median, quantile)
  
  if(missing(peptide.Data.Frame) | missing(labels)){
    
    stop("Error: missing arguments!")
    
  }
  
  
  labels <- as.vector(labels)
  
  if(log == FALSE){
    
    peptide.Data.Frame[labels] <- log2(peptide.Data.Frame[labels])
    
  }
  
  if(method == "average"){
    
    peptide.Data.Frame[labels] <- scale(peptide.Data.Frame[labels], scale = FALSE)
    
  } else if(method == "median"){
    
    median.per.column <- apply(peptide.Data.Frame[labels], 2, median, na.rm = TRUE)
    
    for(i in 1:length(labels)){
      
      peptide.Data.Frame[labels[i]] <- peptide.Data.Frame[labels[i]] - median.per.column[labels[i]]
      
    }
    
  } else if(method == "quantile"){
    
    peptide.Data.Frame[labels] <- preprocessCore::normalize.quantiles(as.matrix(peptide.Data.Frame[labels]))
    
    peptide.Data.Frame[labels] <- scale(peptide.Data.Frame[labels], scale = FALSE)

  } else {
    
    stop("Error: zero center normalization method is not valid!")
    
  }
  
  return(peptide.Data.Frame)
  
}

#Function that returns the protein groups. Combines two types of graph initiation.
protein.Grouping <- function(peptide.Data.Frame, fasta, parsimony = "soft"){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #fasta: a fasta file that contains protein sequences. Peptides will be searched against that file.
  #parsimony: either "soft", "strict" or "none".
  #protien.Groups: a list of protein groups and their corresponding peptides derived by the parsimonious protein inference algorithm. 
  
  if(missing(peptide.Data.Frame)){
  
    stop("Error: required input argument is missing!")  
    
  }
  
  if(missing(fasta)){ #Since the app requires CSV file + FASTA file for all cases, this condition is never TRUE.
    
    grouping.Graph <- mapping.By.Accession(peptide.Data.Frame = peptide.Data.Frame)
    
  } else {
    
    grouping.Graph <- mapping.By.Fasta(peptide.Data.Frame = peptide.Data.Frame, fasta = fasta)
    
  }
  
  if(parsimony == "soft"){
    
    #protein.Groups <- soft.Parsimony(grouping.Graph)
    protein.Groups <- fast.soft.Parsimony(grouping.Graph)
    
  } else if(parsimony == "strict"){
    
    #protein.Groups <- strict.Parsimony(grouping.Graph)
    protein.Groups <- fast.strict.Parsimony(grouping.Graph)
    
  } else if(parsimony == "no"){
    
    #protein.Groups <- no.Parsimony(grouping.Graph)
    protein.Groups <- fast.no.Parsimony(grouping.Graph)
    
  } else {
    
    stop("Error: Not a valid parsimony method!")
    
  }
  
  return(protein.Groups)
}

#Function that searches peptide sequences against a FASTA file specified by user and generates an undirected graph.
mapping.By.Fasta <- function(peptide.Data.Frame, fasta){
  
  #peptide.Data.Frame: a data frame that contains peptide quantitative data.
  #fasta: a fasta file that contains protein sequences. Peptides will be searched against that file.
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  
  if(missing(peptide.Data.Frame) | missing(fasta)){
    
    stop("Error: Missing arguments!")
    
  }
  
  shiny::withProgress(message = "Peptide mapping: ", min = 0, max = 4, value = 0, {
    
    peptide.Sequences <- sort(unique(peptide.Data.Frame[,1]))
    temp_sequences_file <- tempfile(pattern = "peptide_sequences", fileext = ".txt")
    write.table(peptide.Sequences, temp_sequences_file, row.names = F, col.names = F, quote = F)
    
    fasta <- Biostrings::readAAStringSet(fasta, format = "fasta", use.names=TRUE)
    fasta <- chartr(old = "I", new = "L", fasta)
    temp_fasta_file <- tempfile(pattern = "fasta_file", fileext = ".txt")
    Biostrings::writeXStringSet(fasta, temp_fasta_file, format = "fasta")
    
    shiny::incProgress(1, detail = "files prepared")
  
    temp_mapping_file <- tempfile(pattern = "mapping", fileext = ".csv")
    
    cat("\n")
    cat("Peptide mapping started:\n")
    cat("\n")
    
    #Separator for lib in class path (Windows)
    lib.separator <- ";"
    
    if(substr(Sys.info()['sysname'], 1, 1) != "W"){
     
      lib.separator <- ":" #(Linux)
      
    }
    
    system(paste("java -cp utilities-4.12.9/utilities-4.12.9.jar", lib.separator, "utilities-4.12.9/* com.compomics.util.experiment.identification.protein_inference.executable.PeptideMapping -p ", temp_fasta_file, " ", temp_sequences_file, " ", temp_mapping_file, sep = "", collapse = ""))

    cat("Peptide mapping finished!\n")
    
    mapping <- read.csv(temp_mapping_file, header = F, sep = ",", stringsAsFactors = FALSE, comment.char = "#")
    colnames(mapping) <- c("Peptide", "Protein", "Position")

    fasta <- protr::readFASTA(file = temp_fasta_file, legacy.mode = TRUE, seqonly = FALSE)
    fasta <- matrix(cbind(sub(".*[|]([^.]+)[|].*", "\\1", names(fasta)), unlist(fasta)), ncol = 2)
  
    protein.Sequences <- unlist(lapply(mapping[,2], function(x) fasta[which(fasta[,1] == x),2]))
  
    mapping <- data.frame(mapping, protein.Sequences, stringsAsFactors = FALSE)

    duplicates <- which(duplicated(mapping[,c(1,2)]))
  
    if(length(duplicates) > 0){
    
      mapping <- mapping[-duplicates,]
  
    }
    
    unknown <- which(grepl("X", mapping[,1]))
    
    
    if(length(unknown) > 0){
      
      mapping <- mapping[-unknown,]
      
    }

    rownames(mapping) <- c(1:length(mapping$Peptide))
    mapping$Position <- mapping$Position + 1

    file.remove(temp_mapping_file, temp_fasta_file, temp_sequences_file, paste(temp_fasta_file, ".cui", sep = "", collapse = ""))
    
    shiny::incProgress(1, detail = "post-processing is done")
  
    unique.indices <- !duplicated(mapping$Protein)
    attributes <- data.frame(name = c(mapping$Protein[unique.indices], unique(mapping$Peptide)), pname = c(mapping$Protein[unique.indices], unique(mapping$Peptide)), sequence = c(mapping$protein.Sequences[unique.indices], rep(NA,length(unique(mapping$Peptide)))), stringsAsFactors = F)
    relations <- data.frame(from = mapping$Protein, to = mapping$Peptide, stringsAsFactors = F)
  
    graph <- igraph::graph_from_data_frame(relations, directed = F, vertices = attributes)
    
    shiny::incProgress(1, detail = "undirected graph is generated")
    
    graph_n_mapping <- list(graph = graph, mapping = mapping)
  
  })
  
  return(graph_n_mapping)
  
}

#Function that creates an undirected graph based on given protein accession numbers. #NOT USED
mapping.By.Accession <- function(peptide.Data.Frame){
  
  #Deprecated for now.
  
}

#Implementation for a "strict" parsimonious algorithm for protein inference.
#Operating on the graph
strict.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of the parsimonious algorithm.
  
  graph <- graph_n_mapping$graph

  proteins <- igraph::V(graph)$pname[!is.na(igraph::V(graph)$sequence)]
  peptides <- igraph::V(graph)$pname[is.na(igraph::V(graph)$sequence)]
  cc <- igraph::clusters(graph)
  
  groups <- list(proteins = NULL, peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)
  
  shiny::withProgress(message = "Protein Inference: ", min = 0, max = cc$no,value = 0, {
  
    cat("\n")
    cat("Strict parsimony started:\n")
    
    for(i in 1:cc$no){
    
      if(i %% 100 == 0){
        
        cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
        
      }
    
      peptides.In.Group <- intersect(igraph::V(graph)$pname[cc$membership == i], peptides)
      proteins.In.Group <- intersect(igraph::V(graph)$pname[cc$membership == i], proteins)
    
      if(length(proteins.In.Group) == 1){
      
        groups$proteins[length(groups$proteins) + 1] <- proteins.In.Group
        groups$peptides[length(groups$proteins)] <- list(peptides.In.Group)
      
      
      } else if(length(proteins.In.Group) > 1){
      
        reported <- vector()
        
        while(length(peptides.In.Group) != 0){
        
          remaining.proteins <- setdiff(proteins.In.Group, reported)
        
          number.of.peptides <- c(NA)
        
          for(k in 1:length(remaining.proteins)){
          
            number.of.peptides[k] <- length(intersect(peptides.In.Group, igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(remaining.proteins[k], igraph::V(graph)$pname)])))]))
          
          }
        
          peptides.per.protein <- data.frame(name = remaining.proteins, edges = number.of.peptides, stringsAsFactors = FALSE)
        
          peptides.per.protein <- peptides.per.protein[order(peptides.per.protein$edges, decreasing = TRUE),]
        
          p <- peptides.per.protein$name[1]
          reported[length(reported) + 1] <- p
        
          peptides.to.group <- intersect(peptides.In.Group,igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(p, igraph::V(graph)$pname)])))])
          peptides.to.group <- sort(peptides.to.group, decreasing = FALSE)
      
          if(nrow(peptides.per.protein) > 1){
          
            for(m in 2:nrow(peptides.per.protein)){
            
              if(peptides.per.protein$edges[m] == peptides.per.protein$edges[1]){
              
                adjust.peps <- intersect(peptides.In.Group, igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(peptides.per.protein$name[m], igraph::V(graph)$pname)])))])
                adjust.peps <- sort(adjust.peps, decreasing = FALSE)
              
                if(all(peptides.to.group == adjust.peps)){
                
                  p <- paste(c(p, peptides.per.protein$name[m]), collapse = "|")
                  reported[length(reported) + 1] <- peptides.per.protein$name[m]
                  
                }
              }
            }
          }
        
          # print(paste0("Compoment: ", i))
          # print(paste0("Peps in group: ", length(peptides.In.Group)))
          # print(paste0("Prots in group: ", length(remaining.proteins)))
          # print("           ")
          
          peptides.In.Group <- setdiff(peptides.In.Group, peptides.to.group)
          groups$proteins[length(groups$proteins) + 1] <- p
          groups$peptides[length(groups$proteins)] <- list(peptides.to.group)
        
        }
      }
    
      shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", cc$no, " connected components.", collapse = ""))
    }
    
    cat("Strict parsimony finished.\n")
  
  })
  
  return(groups)
}

#Common R structures
fast.strict.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of the parsimonious algorithm.
  
  graph <- graph_n_mapping$graph
  
  #Get connected component subgraphs
  cc.list <- decompose.graph(graph)
  
  #Connected component membership for plotting
  cc <- igraph::clusters(graph)
  
  groups <- list(proteins = vector(mode = "character", length = 0), peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)
  
  shiny::withProgress(message = "Protein Inference: ", min = 0, max = length(cc.list), value = 0, {
  
  cat("\n")
  cat("Strict parsimony started:\n")
  
  for(i in 1:length(cc.list)){
    
    if(i %% 100 == 0){
      
      cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
      
    }
    
    #Graph to DF
    cc.DF <- igraph::as_data_frame(cc.list[[i]], what = c("edges"))
    
    peptides.In.Group <- unique(cc.DF$to)
    proteins.In.Group <- unique(cc.DF$from)
    
    if(length(proteins.In.Group) == 1){
      
      groups$proteins <- c(groups$proteins, proteins.In.Group)
      groups$peptides[length(groups$proteins)] <- list(sort(peptides.In.Group, decreasing = FALSE))
      
      
    } else if(length(proteins.In.Group) > 1){
      
      reported.Proteins <- vector(mode = "numeric", length = 0)
      
      while(length(peptides.In.Group) != 0){
        
        remaining.Proteins <- setdiff(proteins.In.Group, reported.Proteins)
        
        peptides.Per.Protein <- lapply(remaining.Proteins, function(x) cc.DF$to[cc.DF$from == x])
        names(peptides.Per.Protein) <- c(remaining.Proteins)
        
        rank <- data.frame(accession = names(peptides.Per.Protein), number.Of.Peptides = lengths(peptides.Per.Protein))
        rank <- rank[order(rank$number.Of.Peptides, decreasing = TRUE), ]
        
        protein.To.Group <- rank$accession[1]
        reported.Proteins <- c(reported.Proteins, protein.To.Group)
        
        peptides.To.Group <- sort(unlist(peptides.Per.Protein[[rank$accession[1]]]), decreasing = FALSE)
        
        if(nrow(rank) > 1){
          
          #Check if the other leading proteins have the same exactly peptides like the reported protein.
          candidates.To.Merge <- rank$accession[c(FALSE, rank$number.Of.Peptides[2:nrow(rank)] ==  rank$number.Of.Peptides[1])]
          which.To.Merge <- unlist(lapply(candidates.To.Merge, function(x) all(peptides.Per.Protein[[x]] %in% peptides.To.Group)))
          
          if(sum(which.To.Merge) >= 1){
            
            protein.To.Group <- paste(sort(c(protein.To.Group, candidates.To.Merge[which.To.Merge]), decreasing = FALSE), collapse = "|")
            reported.Proteins <- c(reported.Proteins, candidates.To.Merge[which.To.Merge])
            
          }
          
        }
        
        #Reduce peptide pool
        peptides.In.Group <- setdiff(peptides.In.Group, peptides.To.Group)
        
        #Remove the peptides from the cc.DF to reduce search time
        cc.DF <- cc.DF[-unlist(sapply(peptides.To.Group, function(x) which(cc.DF$to == x))), ]
        
        #Report protein group and corresponding peptides
        groups$proteins <- c(groups$proteins, protein.To.Group)
        groups$peptides[length(groups$proteins)] <- list(peptides.To.Group)
        
      }
      
    }
    
    shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", length(cc.list), " connected components.", collapse = ""))
    
  }
  
  cat("Strict parsimony finished.\n")
  
  })
  
  return(groups)
  
}

#Implementation for a "soft" parsimonious algorithm for protein inference.
#Operating on the graph
soft.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of the parsimonious algorithm.
  
  graph <- graph_n_mapping$graph
  
  proteins <- igraph::V(graph)$pname[!is.na(igraph::V(graph)$sequence)]
  peptides <- igraph::V(graph)$pname[is.na(igraph::V(graph)$sequence)]
  cc <- igraph::clusters(graph)

  groups <- list(proteins = NULL, peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)
  
  shiny::withProgress(message = "Protein Inference: ", min = 0, max = cc$no,value = 0, {
    
    cat("\n")
    cat("Soft parsimony started:\n")
  
    for(i in 1:cc$no){
    
      if(i %% 100 == 0){
        
        cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
        
      }
    
      peptides.In.Group <- intersect(igraph::V(graph)$pname[cc$membership == i], peptides)
      proteins.In.Group <- intersect(igraph::V(graph)$pname[cc$membership == i], proteins)
    
      if(length(proteins.In.Group) == 1){
      
        groups$proteins[length(groups$proteins) + 1] <- proteins.In.Group
        groups$peptides[length(groups$proteins)] <- list(peptides.In.Group)
      
      } else if(length(proteins.In.Group) > 1){
      
        reported.proteins <- vector()
        reported.peptides <- vector()
      
        while(length(setdiff(proteins.In.Group, reported.proteins)) > 0){
        
          candidate.proteins <- setdiff(proteins.In.Group, reported.proteins)
        
          number.of.peptides <- c(NA)
        
          for(k in 1:length(candidate.proteins)){
          
            number.of.peptides[k] <- length(intersect(peptides.In.Group, igraph::V(graph)$pname[as.integer(neighbors(graph, as.integer(igraph::V(graph)[match(candidate.proteins[k], igraph::V(graph)$pname)])))]))
              
          }
        
          peptides.per.protein <- data.frame(name = candidate.proteins, edges = number.of.peptides, stringsAsFactors = FALSE)
        
          peptides.per.protein <- peptides.per.protein[order(peptides.per.protein$edges, decreasing = TRUE),]
        
          p <- peptides.per.protein$name[1]
        
          current.peptides <- igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(p, igraph::V(graph)$pname)])))]
          current.peptides <- sort(current.peptides, decreasing = FALSE)
        
          if(nrow(peptides.per.protein) > 1){
          
            for(m in 2:nrow(peptides.per.protein)){
            
              adjust.peps <- igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(peptides.per.protein$name[m], igraph::V(graph)$pname)])))]
              adjust.peps <- sort(adjust.peps, decreasing = FALSE)
            
              if(length(current.peptides) == length(adjust.peps)){
              
                if(all(current.peptides == adjust.peps)){
                
                  p <- paste(c(p, peptides.per.protein$name[m]), collapse = "|")
                
                }
              }
            
              if(length(setdiff(adjust.peps, current.peptides)) == 0){
              
                reported.proteins[length(reported.proteins) + 1] <- peptides.per.protein$name[m]
              
              }
            }
          }
        
          if(length(setdiff(current.peptides, reported.peptides)) > 0){
          
            groups$proteins[length(groups$proteins) + 1] <- p
            groups$peptides[length(groups$proteins)] <- list(current.peptides)
          
            reported.peptides <- union(reported.peptides, current.peptides)
          
            peptides.In.Group <- setdiff(peptides.In.Group, current.peptides)
          
          }
        
          reported.proteins[length(reported.proteins) + 1] <- p 
        
        }
      }
    
      shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", cc$no, " connected components.", collapse = ""))
    }
    
    cat("Soft parsimony finished.\n")
    
  })
  
  return(groups)
}

#Common R structures
fast.soft.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of the parsimonious algorithm.
  
  graph <- graph_n_mapping$graph
  
  #Get connected component subgraphs
  cc.list <- decompose.graph(graph)
  
  #Connected component membership for plotting
  cc <- igraph::clusters(graph)
  
  groups <- list(proteins = vector(mode = "character", length = 0), peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)
  
  shiny::withProgress(message = "Protein Inference: ", min = 0, max = length(cc.list), value = 0, {
  
  cat("\n")
  cat("Soft parsimony started:\n")
  
  rank.list <- list()
  
  for(i in 1:length(cc.list)){
    
    if(i %% 100 == 0){
      
      cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
      
    }
    
    #Graph to DF
    cc.DF <- igraph::as_data_frame(cc.list[[i]], what = c("edges"))
    
    peptides.In.Group <- unique(cc.DF$to)
    proteins.In.Group <- unique(cc.DF$from)
    
    if(length(proteins.In.Group) == 1){
      
      groups$proteins <- c(groups$proteins, proteins.In.Group)
      groups$peptides[length(groups$proteins)] <- list(sort(peptides.In.Group, decreasing = FALSE))
      
      
    } else if(length(proteins.In.Group) > 1){
      
      reported.Proteins <- vector(mode = "numeric", length = 0)
      reported.Peptides <- vector(mode = "numeric", length = 0)
      
      while(length(setdiff(proteins.In.Group, reported.Proteins)) > 0){
        
        candidate.proteins <- setdiff(proteins.In.Group, reported.Proteins)
        
        peptides.Per.Protein <- lapply(candidate.proteins, function(x) cc.DF$to[cc.DF$from == x])
        names(peptides.Per.Protein) <- c(candidate.proteins)
        
        edges <- lapply(peptides.Per.Protein, function(x) length(intersect(peptides.In.Group, x)))
        names(edges) <- c(candidate.proteins)
        
        rank <- data.frame(accession = names(edges), number.Of.Peptides = unlist(edges))
        rank <- rank[order(rank$number.Of.Peptides, decreasing = TRUE), ]
        
        protein.To.Group <- rank$accession[1]
        cc.DF <- cc.DF[!grepl(protein.To.Group, cc.DF$from), ]
        reported.Proteins[length(reported.Proteins) + 1] <- protein.To.Group
        
        current.peptides <- sort(unlist(peptides.Per.Protein[[rank$accession[1]]]), decreasing = FALSE)
        
        if(nrow(rank) > 1){
          
          for (m in 2:nrow(rank)){
            
            adjusent.peps <- peptides.Per.Protein[[rank$accession[m]]]
            
            if(length(setdiff(adjusent.peps, current.peptides)) == 0){
              
              print(adjusent.peps)
              print(current.peptides)
              
              protein.To.Group <- paste(c(protein.To.Group, rank$accession[m]), collapse = "|")
              reported.Proteins[length(reported.Proteins) + 1] <- rank$accession[m]
              cc.DF <- cc.DF[!grepl(rank$accession[m], cc.DF$from), ]
              
            }
            
            
          }
          
        }
        
        if(length(setdiff(current.peptides, reported.Peptides)) > 0){
          
          groups$proteins[length(groups$proteins) + 1] <- protein.To.Group
          groups$peptides[length(groups$proteins)] <- list(current.peptides)
          
          reported.Peptides <- union(reported.Peptides, current.peptides)
          
          peptides.In.Group <- setdiff(peptides.In.Group, current.peptides)

        }
        
      }
      
    }
    
    shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", length(cc.list), " connected components.", collapse = ""))
    
  }
  
  cat("Strict parsimony finished.\n")
  
  })
  
  return(groups)
  
}

#Function that creates protein groups without using any protein inference algorithm.
#Operating on the graph
no.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of no parsimony.
  
  graph <- graph_n_mapping$graph
  
  proteins <- igraph::V(graph)$pname[!is.na(igraph::V(graph)$sequence)]
  
  cc <- igraph::clusters(graph)
  
  groups <- list(proteins = NULL, peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)

  shiny::withProgress(message = "Protein Inference: ", min = 0, max = cc$no,value = 0, {
    
    cat("\n")
    cat("No parsimony started:\n")
  
    for(i in 1:cc$no){
    
      if(i %% 100 == 0){
        
        cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
        
      }
  
      proteins.In.Group <- intersect(igraph::V(graph)$pname[cc$membership == i], proteins)
    
      for(j in 1:length(proteins.In.Group)){
      
        groups$proteins[length(groups$proteins) + 1] <- proteins.In.Group[j]
      
        adjust.peps <- igraph::V(graph)$pname[as.integer(igraph::neighbors(graph, as.integer(igraph::V(graph)[match(proteins.In.Group[j], igraph::V(graph)$pname)])))]
        adjust.peps <- sort(adjust.peps, decreasing = FALSE)
      
        groups$peptides[length(groups$proteins)] <- list(adjust.peps)
      
      }
    
      shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", cc$no, " connected components.", collapse = ""))
    }
    
    cat("No parsimony finished.\n")
    
  })
  
  return(groups)
}

#Common R structures
fast.no.Parsimony <- function(graph_n_mapping){
  
  #graph: an undirected graph in which every peptide and protein are represented by nodes and their relation by edges.
  #groups: a minimal list of protein groups as a result of the parsimonious algorithm.
  
  graph <- graph_n_mapping$graph
  
  #Get connected component subgraphs
  cc.list <- decompose.graph(graph)
  
  #Connected component membership for plotting
  cc <- igraph::clusters(graph)
  
  groups <- list(proteins = vector(mode = "character", length = 0), peptides = list(), cc = cc, g = graph, mapping = graph_n_mapping$mapping)
  
  shiny::withProgress(message = "Protein Inference: ", min = 0, max = cc$no,value = 0, {
  
  cat("\n")
  cat("No parsimony started:\n")
  
  for(i in 1:cc$no){
    
    if(i %% 100 == 0){
      
      cat(paste("| ", i, " connected components processed.\n", collapse = "", sep = ""))
      
    }
    
    #Graph to DF
    cc.DF <- igraph::as_data_frame(cc.list[[i]], what = c("edges"))
    
    proteins.In.Group <- unique(cc.DF$from)
    
    peptides.Per.Protein <- lapply(proteins.In.Group, function(x) sort(cc.DF$to[cc.DF$from == x], decreasing = F))
    names(peptides.Per.Protein) <- c(proteins.In.Group)
    
    #Report protein group and corresponding peptides
    groups$proteins <- c(groups$proteins, proteins.In.Group)
    groups$peptides <- c(groups$peptides, peptides.Per.Protein)
    
    shiny::incProgress(amount = 1, detail = paste("processed ", i, " out of ", cc$no, " connected components.", collapse = ""))
  }
  
  cat("No parsimony finished.\n")
  
  })
  
  return(groups)
}

#Implementation of fast-FARMS function as it is used in Diffacto.
fast.Farms <- function(probes, weight, mu, max_iter, force_iter, min_noise, fill_nan){
  
  if(missing(weight)){
    
    weight <- 0.5
    
  }	
  
  if(missing(mu)){
    
    mu <- 0
    
  }	
  
  if(missing(max_iter)){
    
    max_iter <- 1000
    
  }	
  
  if(missing(force_iter)){
    
    force_iter <- FALSE
    
  }
  
  if(missing(min_noise)){
    
    min_noise <- 0.0001
    
  }
  
  if(missing(fill_nan)){
    
    fill_nan <- 0
    
  }
  
  readouts <- as.matrix(probes)
  
  readouts[is.na(readouts)] <- fill_nan
  
  #Normalize and transform X
  X <- t(readouts)
  X <- t(t(X) - colMeans(X, na.rm = T))
  xsd <- apply(X, 2, function(x) sd(x, na.rm = T) * sqrt((length(x) - 1) / length(x)))
  xsd[xsd < min_noise] <- 1
  X <- t(t(X)/xsd)
  X[!is.finite(X)] <- 0
  
  n_samples <- nrow(X)
  n_features <- ncol(X)
  C <- crossprod(X, X)/n_samples
  
  #Positive definite
  C <- (C+t(C))/2
  C[which(C < 0)] <- 0
  
  #Robustness
  SVD <- svd(C)
  U <- SVD$u
  s <- SVD$d
  V <- t(SVD$v)
  s[s<min_noise] <- min_noise
  C <- U %*% diag(s) %*% V
  diag(C) <- abs(diag(C))
  
  #Initiation
  lamda <- sqrt(0.75*diag(C))
  psi <- diag(C) - lamda^2
  old_psi <- psi
  alpha <- weight * n_features
  E <- 1
  
  for(i in 1:max_iter){
    #E Step
    phi <- (1/psi)*lamda
    a <- as.vector(1+crossprod(lamda,phi))
    eta <- phi/a
    zeta <- C %*% eta
    E <- 1 - as.vector(eta) %*% lamda + as.vector(eta) %*% zeta
    
    #M Step
    lamda = zeta/(c(E) + as.vector(psi)*alpha)
    psi <- diag(C) - as.vector(zeta)[1] * lamda + psi * alpha * lamda * (mu - lamda)
    psi[psi < min_noise^2] <- min_noise^2
    
    if(!force_iter){
      
      if(max(abs(psi-old_psi))/max(abs(old_psi)) < min_noise/10){
        
        break
        
      }
    }
    
    old_psi <- psi
    
  }
  
  loading <- as.vector(sqrt(E))*lamda
  phi <- loading/psi
  weights <- loading/max(loading)
  weights <- round(weights,5)
  noise <- 1/as.vector(1+crossprod(loading,phi))
  
  loading.noise <- list("loadings" = weights, "noise" = noise)
  
  return(loading.noise)
  
}

#Function that returns the protein groups, probes (peptide abundances per protein group) and the corresponding assigned weights from fast-FARMS. 
probes.And.Weights <- function(peptide.abundances, protein.groups, labels, rescale, weight, mu, force){
  
  if(missing(peptide.abundances) | missing(protein.groups) | missing(labels)){
    
    stop("Error: One or more required attributes are missing!")
    
  }
  
  if(missing(rescale)){
    
    rescale <- TRUE
    
  }
  
  if(missing(force)){
    
    force <- FALSE
    
  }
  
  group.names <- vector(mode = "character", length = 0)
  
  if(is.element("Reference", colnames(labels))){
    
    group.names <- setdiff(colnames(labels), "Reference")
    
  } else if(!is.element("Reference", colnames(labels))){
    
    group.names <- colnames(labels)
    
  }
  
  probes_n_weights <- list(protein.group = vector(mode = "character", length = length(which(lengths(protein.groups$peptides) > 1))),
                           probes = vector(mode = "list", length = length(which(lengths(protein.groups$peptides) > 1))),
                           peptide.indices = vector(mode = "list", length = length(which(lengths(protein.groups$peptides) > 1))),
                           weights = vector(mode = "list", length = length(which(lengths(protein.groups$peptides) > 1))),
                           noise = vector(mode = "numeric", length = length(which(lengths(protein.groups$peptides) > 1))) )
  
  protein.counter <- 1
  
  shiny::withProgress(message = "Weights calculation: ", min = 0, max = length(protein.groups$proteins), value = 0, {
  
    for(i in 1:length(protein.groups$proteins)){
    
      peps <- unlist(protein.groups$peptides[i])
      
      peps.count <- length(peps)
    
      protein.df <- data.frame(matrix(ncol = length(as.vector(labels)), nrow = peps.count))
      colnames(protein.df) <- as.vector(labels)
    
      p.indices <- vector(mode = "numeric", length = peps.count)
      
      for(j in 1:length(peps)){
      
        p.index <- which(peptide.abundances[,1] == peps[j])
        p.indices[j] <- p.index
        protein.df[j, as.vector(labels)] <- peptide.abundances[p.index, as.vector(labels)]
      
      }
      
      probes.abs <- protein.df
    
      if(length(as.vector(labels)) > length(as.vector(labels[, group.names])) & !rescale){
      
        reference.abundance <- protein.df[, as.vector(labels[, "Reference"])]
        
        if(length(as.vector(labels[, "Reference"])) > 1){
          
          reference.abundance <- rowMeans(reference.abundance, na.rm = T)
          
        } else {
          
          reference.abundance <- reference.abundance
          
        }

        reference.abundance[is.na(reference.abundance)] <- mean(unlist(protein.df), na.rm = T)

      } else if(rescale) {
      
        reference.abundance <- rowMeans(protein.df, na.rm = T)
      
      } else {
      
        reference.abundance <- rep(0, peps.count)
      
      }
    
      protein.df <- protein.df - reference.abundance
      
      if(peps.count == 1){
      
        loading <- 1
        noise <- 1
      
      } else {
        
        farms.results <- fast.Farms(probes =  protein.df, weight = weight, mu = mu, max_iter = 1000, force_iter = force )
        loading <- as.vector(farms.results$loadings)
        noise <- farms.results$noise
      
      }
    
      if(peps.count > 1){
      
        probes_n_weights$protein.group[protein.counter] <- protein.groups$proteins[i]
        protein.df <- cbind(peps, protein.df)
        colnames(protein.df)[1] <- "peptide"
        probes_n_weights$probes[[protein.counter]] <- protein.df
        probes_n_weights$peptide.indices[[protein.counter]] <- p.indices
        probes_n_weights$weights[[protein.counter]] <- loading
        probes_n_weights$noise[protein.counter] <- noise
        
        protein.counter <- protein.counter + 1
      
      }
    
      shiny::incProgress(amount = 1, detail = paste("processed ", i," out of ", length(protein.groups$proteins), " protein groups." , collapse = ""))
      
    }
  })
  
  return(probes_n_weights)
}

#Function that returns the protein quantification based on specific weight threshold.
quantification.Weight.Threshold <- function(probes_n_weights, threshold){
  
  #probes_n_weights: a list object that contains the protein group identifiers, the probes that correspond to these groups and the peptide weights.
  #threshold: a weight threshold cutoff.
  #quantification.df.rep.thr: protein expression for the specific threshold
  
  if(missing(threshold)){
    
    threshold <- 0.7
    
  }
  
  #quantification.df.rep.thr <- t(as.matrix(sapply(c(1:length(probes_n_weights$protein.group)), function(x) expr.Average(probes_n_weights$probes[[x]][2:ncol(probes_n_weights$probes[[x]])]))))
  quantification.df.rep.thr <- t(as.matrix(sapply(c(1:length(probes_n_weights$protein.group)), function(x) expr.Threshold.on.Weights(unlist(probes_n_weights$weights[x]), probes_n_weights$probes[[x]][2:ncol(probes_n_weights$probes[[x]])], threshold = threshold))))
  quantification.df.rep.thr <- data.frame(Protein = probes_n_weights$protein.group, quantification.df.rep.thr)
  colnames(quantification.df.rep.thr)[2:ncol(quantification.df.rep.thr)] <- colnames(probes_n_weights$probes[[1]])[2:ncol(probes_n_weights$probes[[1]])] 
  
  return(quantification.df.rep.thr)
  
}

#Function that performs weighted average summarization with NA values exclusion.
expr.Threshold.on.Weights <- function(weights, peptide_abundances, threshold){
  
  #weights: fast-FARMS peptide weights for a given protein group.
  #peptide_abundances: the abundances of peptides corresponding to the weights.
  #threshold: a weight threshold cutoff.
  #expr: protein expression.
  
  if(missing(threshold)){
    
    threshold = 0.7
    
  }
  
  weights[weights < threshold] <- NA
  
  abd_w <- peptide_abundances * weights
  
  expr <- vector(mode = "numeric", length = ncol(peptide_abundances))
  
  for(i in 1:ncol(abd_w)){
    
    weights_t <- weights
    weights_t[which(is.na(peptide_abundances[i]))] <- NA
    
    if(all(is.na(weights_t))){
      
      expr[i] <- NA
      
    } else {
      
      abd_w_t <- abd_w[i]/sum(weights_t,na.rm = T)
      
      if(all(is.na(abd_w_t))){
        
        expr[i] <- NA
        
      } else {
        
        expr[i] <- sum(abd_w_t, na.rm = T)
        
      }
      
    }
    
  }
  
  return(expr)
  
}

#Function that performs weighted average summarization. #NOT USED
expr.Threshold.on.Weights2 <- function(weights, peptide_abundances, threshold){
  
  #weights: fast-FARMS peptide weights for a given protein group.
  #peptide_abundances: the abundances of peptides corresponding to the weights.
  #threshold: a weight threshold cutoff.
  #expr: protein expression.
  
  if(missing(threshold)){
    
    threshold = 0.7
    
  }
  
  weights[weights < threshold] <- NA
  abd_w <- (peptide_abundances * weights)/sum(weights, na.rm = T)
  
  expr <- vector(mode = "numeric", length = ncol(peptide_abundances))
  
  for(i in 1:ncol(peptide_abundances)){
    
    if(all(is.na(abd_w[i]))){
      
      expr[i] <- NA
      
    } else {
      
      expr[i] <- sum(abd_w[i], na.rm = T)
      
    }
  }
  
  return(expr)
  
}

#Function that performs average summarization. #NOT USED
expr.Average <- function(peptide_abundances){
  
  
  return(colMeans(peptide_abundances, na.rm = T))
  
}

#Function that calls quantification.Weight.Threshold() function for multiple weights in parallel trend.
calculateParallel <- function(probes_n_weights, cores){
  
  #probes_n_weights: list that contains probes and weights for all protein groups.
  #cores: total number of cores to be used.
  #expr: list that contains the protein expressions for all weight cutoffs.
  
  if(missing(probes_n_weights)){
    
    stop("Error: required input attribute is missing!")
    
  } 
  
  if(missing(cores)){
    
    cores <- parallel::detectCores() - 1
    
  }
  
  pnw <- probes_n_weights
  thresholds <- seq(from = 0, to = 1, by = 0.1)
  expr <- vector(mode = "list", length = length(thresholds))
  
  detected.Cores <- parallel::detectCores()
  
  if(detected.Cores <= cores){
    
    if(substr(Sys.info()['sysname'], 1, 1) == "W"){
      
      cl <- parallel::makePSOCKcluster(detected.Cores - 1)
    
    } else {
     
      cl <- parallel::makeCluster(detected.Cores - 1)
       
    }
      
  } else {
    
    if(substr(Sys.info()['sysname'], 1, 1) == "W"){
      
      cl <- parallel::makePSOCKcluster(cores)
    
    } else {

      cl <- parallel::makeCluster(cores)

    }
  }
  
  parallel::setDefaultCluster(cl)
  parallel::clusterExport(cl, c("pnw","thresholds", 'quantification.Weight.Threshold', 'expr.Threshold.on.Weights'), envir = environment())
  
  expr <- parallel::parLapply(cl, thresholds, function(x) quantification.Weight.Threshold(pnw, x))
  
  parallel::stopCluster(cl)
  
  return(expr)
  
}

#Function that summarizes protein expressions by either total summation of peptides or weighted summation of peptides.
quantification.By.Summation <- function(peptide.abundances, labels, log = F, normalization.method = "median", probes_n_weights){
  
  if(missing(peptide.abundances) | missing(probes_n_weights) | missing(labels)){
    
    stop("Error: One or more required arguments are missing!")
    
  }
  
  #Create probes of absolute intensities
  absolute.probes <- lapply(probes_n_weights$peptide.indices, function(x){return(peptide.abundances[x, as.vector(labels)])})
  
  #Summarize proteins by total summation
  expr.total.sum <- do.call(rbind, lapply(absolute.probes, function(x){colSums(x, na.rm = T)}))
  expr.total.sum[expr.total.sum == 0] <- NA
  expr.total.sum <- data.frame(probes_n_weights$protein.group, expr.total.sum)
  colnames(expr.total.sum)[1] <- "Protein"
  rownames(expr.total.sum) <- 1:nrow(expr.total.sum)
  
  #Summarize proteins by weighted summation
  expr.weighted.sum <- do.call(rbind, lapply(1:length(absolute.probes), function(x){ return(colSums(absolute.probes[[x]] * probes_n_weights$weights[[x]], na.rm = T ))}))
  expr.weighted.sum[expr.weighted.sum == 0] <- NA
  expr.weighted.sum <- data.frame(probes_n_weights$protein.group, expr.weighted.sum)
  colnames(expr.weighted.sum)[1] <- "Protein"
  rownames(expr.weighted.sum) <- 1:nrow(expr.weighted.sum)
  
  #Log2 transform
  expr.total.sum[,2:ncol(expr.total.sum)] <- log2(expr.total.sum[,2:ncol(expr.total.sum)])
  expr.weighted.sum[,2:ncol(expr.weighted.sum)] <- log2(expr.weighted.sum[,2:ncol(expr.weighted.sum)])
  
  #Normalization
  if(normalization.method == "average"){

    colMeans.total <- colMeans(expr.total.sum[,2:ncol(expr.total.sum)], na.rm = T)
    expr.total.sum[,2:ncol(expr.total.sum)] <- data.frame(apply(expr.total.sum[,2:ncol(expr.total.sum)], 2, function(x) (x/mean(x, na.rm = T)) * mean(colMeans.total, na.rm = T)))
    colMeans.weighted <- colMeans(expr.weighted.sum[,2:ncol(expr.weighted.sum)], na.rm = T)
    expr.weighted.sum[,2:ncol(expr.weighted.sum)] <- data.frame(apply(expr.weighted.sum[,2:ncol(expr.weighted.sum)], 2, function(x) (x/mean(x, na.rm = T)) * mean(colMeans.weighted, na.rm = T)))

  } else if(normalization.method == "median"){

    colMedian.total <- apply(expr.total.sum[,2:ncol(expr.total.sum)], 2, median, na.rm = T)
    expr.total.sum[,2:ncol(expr.total.sum)] <- data.frame(apply(expr.total.sum[,2:ncol(expr.total.sum)], 2, function(x) (x/median(x, na.rm = T)) * mean(colMedian.total, na.rm = T)))
    colMeadian.weighted <- apply(expr.weighted.sum[,2:ncol(expr.weighted.sum)], 2, median, na.rm = T)
    expr.weighted.sum[,2:ncol(expr.weighted.sum)] <- data.frame(apply(expr.weighted.sum[,2:ncol(expr.weighted.sum)], 2, function(x) (x/median(x, na.rm = TRUE)) * mean(colMeadian.weighted)))

  } else if(normalization.method == "quantile"){

    expr.total.sum[,2:ncol(expr.total.sum)] <- preprocessCore::normalize.quantiles(as.matrix(expr.total.sum[,2:ncol(expr.total.sum)]))
    expr.weighted.sum[,2:ncol(expr.weighted.sum)] <- preprocessCore::normalize.quantiles(as.matrix(expr.weighted.sum[,2:ncol(expr.weighted.sum)]))

  } else {

    stop("Error: normalization method is not valid!")

  }
  
  return(list(total = expr.total.sum, weighted = expr.weighted.sum))
  
}

#Creates a static line plot in ggplot. NOT USED ANYMORE
static.Line.Plot <- function(probes, weights, group, expr, threshold, colors, xlabel){
  
  #probes: dataframe that contains the peptide quantitative data.
  #weights: a vector for peptide weights.
  #group: the protein group name.
  #expr: protein expression data.
  #threshold: FARMS weight threshold that is applied for the input data.
  #colors: a palette with distinctive colors.
  #g: ggplot graph.
  
  if(missing(probes) | missing(weights) | missing(group) | missing(expr) | missing(threshold) | missing(xlabel)){
    
    stop("Error: missing arguments in static.Line.Plot()!")
    
  }
  
  if(length(colors) < (length(weights) + 1)){
    
    colors <- c(colors[1], rep(colors[2:length(colors)], ceiling((length(weights) + 1)/length(colors))))
    warning("Warning: color palette argument in static.Line.Plot() is extended.")
    
  }
  
  name <- c(as.character(probes[,1]), "Z")
  visualization.df <- cbind(Name = name, rbind(probes[2:ncol(probes)], expr[2:length(expr)]))
  
  palette <- colors[2:(length(weights) + 1)]
  palette <- c(palette, colors[1])
  palette[(which(weights < threshold))] <- "#484848"
  
  long.df <- reshape2::melt(visualization.df, id = "Name")
  long.df$row.id <- rep(1:nrow(visualization.df), ncol(visualization.df) - 1)
  long.df$size <- rep(c(rep(0, nrow(visualization.df) - 1), 1),ncol(visualization.df) - 1)
  long.df$opacity <- rep(c(weights,1), ncol(visualization.df) - 1)
  long.df$threshold <- rep(c(as.numeric(weights < threshold), 0), ncol(visualization.df) - 1)
  
  g <- ggplot2::ggplot(data = long.df, ggplot2::aes(x = variable, y = value, group = factor(row.id), colour = Name, size = factor(size), alpha = opacity)) +
                       ggplot2::geom_point() +
                       ggplot2::geom_line(aes(linetype = factor(threshold))) +
                       ggplot2::scale_color_manual(breaks = name, values = palette, labels = c(name[1:(length(name) - 1)], "Protein")) +    
                       ggplot2::scale_size_manual(values = c(0.5, 1), labels = c("Peptide abundance", "Protein abundance")) +
                       #ggplot2::scale_x_discrete(labels = c(1:(ncol(visualization.df) - 1)), expand = c(0, 1)) +
                       ggplot2::scale_x_discrete(labels = colnames(visualization.df[2:ncol(visualization.df)])) +
                       ggplot2::scale_y_continuous(limits = c(min(na.omit(long.df$value)) - 0.1, max(na.omit(long.df$value)) + 0.1), expand = c(0, 0)) +
                       ggplot2::scale_linetype_manual(values = c(1, 2), labels = c("Coherent", "Eliminated")) +
                       ggplot2::labs(linetype = "Peptide status", colour = "Trace colors", size = "Size", alpha = "FARMS weight") +
                       ggplot2::ggtitle(paste("Abundance line plot of protein group", group, "and its corresponding peptides, for FARMS weight threshold", threshold, collapse = "")) +
                       ggplot2::xlab(xlabel) + ggplot2::ylab("Abundance (log2)") + 
                       ggplot2::theme_minimal(base_size = 6) +
                       ggplot2::theme(plot.title = element_text(size = 10))
  
  return(g)
  
}

#Creates an interactive line plot of peptide and protein abundances in plotly (per runs\conditions). 
interactive.Line.Plot <- function(probes, weights, expr, threshold, colors, xlabel, peptide_mapping){
  
  #probes: dataframe that contains the peptide quantitative data.
  #weights: a vector for peptide weights.
  #group: the protein group name.
  #expr: protein expression data.
  #threshold: FARMS weight threshold that is applied for the input data.
  #colors: a palette with distinctive colors.
  #p: plotly graph.
  
  if(missing(probes) | missing(weights) | missing(expr) | missing(threshold) | missing(xlabel)){
    
    stop("Error: missing arguments in static.Line.Plot()!")
    
  }
  
  
  if(length(colors) < (length(weights) + 1)){
    
    colors <- c(colors[1], rep(colors[2:length(colors)], ceiling((length(weights) + 1)/length(colors))))
    warning("Warning: color palette argument in static.Line.Plot() extended.")
    
  }
  
  visualization.df <- rbind(expr[2:length(expr)], probes[2:ncol(probes)])
  
  opacity <- round(weights, 1)
  opacity[opacity == 0] <- 0.1
  peptides <- as.character(probes[,1])
  
  peptide_positions <- peptide_mapping$Position[which(peptide_mapping$Protein == unlist(strsplit(as.character(expr[1,1]),"\\|"))[1])]
  palette <- c(colors[1], colors[(order(peptide_positions)+1)])
  
  #palette <- colors[1:(nrow(visualization.df))]
  palette[(which(weights < threshold)+1)] <- "#484848"
  linetype <- rep("solid",length(peptides))
  linetype[which(weights < threshold)] <- "dot"
  
  p <- plotly::plot_ly()
  
  if(!all(is.na(as.numeric(visualization.df[1, ])))){
    
    ribbons_min <- unlist(lapply(c(1:ncol(visualization.df)), function(x) if(length(which(is.finite(visualization.df[(which(weights >= threshold) + 1), x]))) > 1){
      
      as.numeric(visualization.df[1, x]) - sd(visualization.df[(which(weights >= threshold) + 1), x], na.rm = T)
    
    } else {as.numeric(visualization.df[1, x])}))
    
    ribbons_max <- unlist(lapply(c(1:ncol(visualization.df)), function(x) if(length(which(is.finite(visualization.df[(which(weights >= threshold) + 1), x]))) > 1){
      
      as.numeric(visualization.df[1, x]) + sd(visualization.df[(which(weights >= threshold) + 1), x], na.rm = T)
      
    } else {as.numeric(visualization.df[1, x])}))
    
    p <- plotly::add_ribbons(p, 
                             x = c(1:ncol(visualization.df)), 
                             y = as.numeric(visualization.df[1,]),
                             ymin = ribbons_min,
                             ymax = ribbons_max,
                             opacity = 0.05,
                             line = list(color = "#ff0066"),
                             fillcolor = "#ff0066",
                             showlegend = FALSE,
                             hoverinfo = "none")
    
  }
  
  for(i in 2:nrow(visualization.df)){
    
    p <-plotly::add_trace(p, 
                          x = c(1:ncol(visualization.df)), 
                          y = as.numeric(visualization.df[i,]),
                          opacity = opacity[i-1],
                          name = peptides[i-1], 
                          type = "scatter", 
                          mode = "lines+markers",
                          line = list(color = palette[i], width = 1, dash = linetype[i-1]),
                          marker = list(color = palette[i], size = 3),
                          visible = T,
                          hoverinfo = "text",
                          text = unlist(lapply(c(1:ncol(visualization.df)), function(x) paste(" Peptide: ", peptides[i-1], "<br>", 
                                                                                              "Weight: ", weights[i-1], "<br>", 
                                                                                              "Abundance: ", visualization.df[i, x], "<br>", 
                                                                                              gsub(".$", "", xlabel) ,": ", colnames(visualization.df)[x]))))
  }
  
  p <- plotly::add_trace(p, 
                         x = c(1:ncol(visualization.df)),
                         y = as.numeric(visualization.df[1,]),
                         opacity = 1,
                         name = "Protein",
                         type = "scatter",
                         mode = "lines+markers",
                         line = list(color = palette[1], width = 2),
                         marker = list(color = palette[1], size = 4),
                         visible = T,
                         hoverinfo = "text",
                         text = unlist(lapply(c(1:ncol(visualization.df)), function(x) paste(" Protein: ", expr[1, 1], "<br>",
                                                                                             "Adundance: ", visualization.df[1, x], "<br>",
                                                                                             gsub(".$", "", xlabel), ": ", colnames(visualization.df)[x]))))
  
  
  plot.Margin <- list(l = 50, r = 10, b = 50, t = 50, pad = 0)
  
  p <- plotly::layout(p, 
                      title = paste("Line plot of protein ", expr[1, 1]," and peptide abundances.", sep = "", collapse = ""),
                      font = list(size = 10),
                      margin = plot.Margin,
                      xaxis = list(title = xlabel, tickvals=c(1:(ncol(probes) - 1)), ticktext = colnames(probes[2:ncol(probes)]), size = 9),
                      yaxis = list(title = "Abundance (log2)", size = 9),
                      hoverlabel = list(namelength = -1))
  
  p <- plotly::config(p, 
                      displayModeBar = T,
                      showLink = F, 
                      displaylogo = F)
  
  p$elementId <- NULL
  
  return(p)
  
}

#Creates an interactive VIQoR plot for protein and peptide abundance fold change between two conditions.
VIQoR.Plot <- function(protein_group_index, condition1, condition2, protein_expression, probes_n_weights, threshold, peptide_mapping, pep_colors, modified_peptide_mapping, modifications, data, mod_colors, labels, enzyme = "None"){
  
  #condition1: control
  #condition2: treatment
  
  if(missing(protein_group_index) | missing(condition1) | missing(condition2) | missing(protein_expression) | missing(probes_n_weights) | missing(threshold) | missing(peptide_mapping)){
    
    stop("Error: missing input argument in VIQoR.Plot().")
    
  }
  
  if(missing(labels)){
    
    p_expr1 <- protein_expression[[(threshold * 10) + 1]][protein_group_index, condition1]
    p_expr2 <- protein_expression[[(threshold * 10) + 1]][protein_group_index, condition2]
    
    peptide_expr1 <- probes_n_weights$probes[[protein_group_index]][, condition1]
    peptide_expr2 <- probes_n_weights$probes[[protein_group_index]][, condition2]
    
    type.legend <- "runs"
    
  } else {
    
    p_expr1 <- protein_expression[protein_group_index, condition1]
    p_expr2 <- protein_expression[protein_group_index, condition2]
    
    probes <- sapply(c(1:ncol(labels)), function(x) rowMeans(probes_n_weights$probes[[protein_group_index]][, labels[, x]], na.rm = T))
    colnames(probes) <- colnames(labels)
    
    peptide_expr1 <- probes[, condition1]
    peptide_expr2 <- probes[, condition2]
    
    type.legend <- "conditions"
    
  }

  
  if(is.na(p_expr2 - p_expr1)){
    
    #stop("Error: the protein expression of one or both comparing conditions doesn't exist.")
    return(NULL)
    
  }
    
  protein_name <- unlist(strsplit(probes_n_weights$protein.group[protein_group_index], "\\|"))[1]
    
  protein_sequence <- peptide_mapping$protein.Sequences[which(peptide_mapping$Protein == protein_name)[1]]
  
  peptide_names <- as.character(probes_n_weights$probes[[protein_group_index]][, 1])
  peptide_weights <- probes_n_weights$weights[[protein_group_index]]
    
  #Peptide mapping =\= to peptides in protein group. Should keep only these that belong to protein group from the mapped peptides.
  all_mapped_peptides <- peptide_mapping[which(peptide_mapping$Protein == protein_name), ]
  
  peptide_difference <- setdiff(all_mapped_peptides$Peptide, peptide_names)
  
  if(length(peptide_difference) > 0){
    
    peptides_to_remove <- unlist(lapply(1:length(peptide_difference), function(x) which(all_mapped_peptides$Peptide == peptide_difference[x])))
    peptide_positions <- all_mapped_peptides$Position[-peptides_to_remove]
    
  } else {
    
    peptide_positions <- all_mapped_peptides$Position
    
  }
  
  if(length(pep_colors) < (length(peptide_weights) + 1)){
    
    pep_colors <- c(pep_colors[1], rep(pep_colors[2:length(pep_colors)], ceiling((length(peptide_weights) + 1)/length(pep_colors))))
    warning("Warning: color palette argument in VIQoR.Plot() extended.")
    
  }
  
  peptide_colors <- pep_colors[2:(length(peptide_names) + 1)]
  peptide_colors <- c(pep_colors[1], peptide_colors[order(peptide_positions)])
    
  colors <- c("#33ccff", rep("#66ff66",length(peptide_names)))
  
  #Find missed cleavages
  if(enzyme != "None"){
    
    enzyme <- dplyr::case_when(
      
      enzyme == "Trypsin" ~ "trypsin",
      enzyme == "Trypsin Strict" ~ "trypsin.strict",
      enzyme == "Chymotrypsin (High)" ~ "chymotrypsin.h",
      enzyme == "Chymotrypsin (Low)" ~ "chymotrypsin.l",
      enzyme == "Pepsin (pH1.3)" ~ "pepsin.2",
      enzyme == "Pepsin (pH>2)" ~ "pepsin.1.3",
      enzyme == "Lys-C" ~ "lysC",
      enzyme == "Arg-C" ~ "argC"
      
    )
    
    digest <- lapply(peptide_names, function(x) fastDigest(sequence = x, enzyme = enzyme))
    digest.length <- lengths(digest)
    colors[(which(digest.length > 0) + 1)] <- "#e8a917"
    
  }

  colors[(which(peptide_weights < threshold) + 1)] <- "#484848"

  log2FC.df <- data.frame(Sequences = c(protein_sequence, peptide_names), 
                          Length = unlist(sapply(c(protein_sequence, peptide_names), function(x) nchar(x))), 
                          Position = c(1, peptide_positions), 
                          Colors = colors,
                          Ratio = c(p_expr2 - p_expr1, peptide_expr2 - peptide_expr1), 
                          Opacity = c(1, peptide_weights),
                          PeptideColors = peptide_colors,
                          stringsAsFactors = F)
  
  coverage.vector <- rep(0,log2FC.df$Length[1])
  
  p <- plotly::plot_ly()
    
  for(i in 1:(nrow(log2FC.df))){
      
    if(is.na(log2FC.df$Ratio[i])) next
      
    p <- plotly::add_trace(p,
                           type = "scatter",
                           x = c(log2FC.df$Position[i], log2FC.df$Position[i] + log2FC.df$Length[i]),
                           y = c(log2FC.df$Ratio[i], log2FC.df$Ratio[i]),
                           mode = "lines",
                           line = list(color = log2FC.df$Colors[i], width = 10),
                           showlegend = F,
                           hoverinfo = "text",
                           opacity = max(log2FC.df$Opacity[i], 0.1),
                           text = paste("Sequence: ", if(log2FC.df$Length[i] > 150){gsub(".{120}$", "", paste(substring(log2FC.df$Sequences[i], seq(1, log2FC.df$Length[i], 150), if(log2FC.df$Length[i] %% 150 == 0){seq(150, log2FC.df$Length[i], 150)} else {c(seq(150, log2FC.df$Length[i], 150), log2FC.df$Length[i])}), paste("<br>", paste(rep("&nbsp;", 20), sep = "", collapse = ""), sep = "", collapse = "") ,sep = "", collapse = ""))} else {paste(log2FC.df$Sequences[i], "<br>")},
                                        "Ratio: ", log2FC.df$Ratio[i], "<br>",
                                        "Position: ", log2FC.df$Position[i], "<br>",
                                        "Length: ", log2FC.df$Length[i], "<br>",
                                        "Weight: ", log2FC.df$Opacity[i])
      )
    
    if(i != 1){
      
      p <- plotly::add_trace(p,
                             type = "scatter",
                             x = c(log2FC.df$Position[i] + (log2FC.df$Length[i] - 2), log2FC.df$Position[i] + log2FC.df$Length[i]),
                             y = c(log2FC.df$Ratio[i], log2FC.df$Ratio[i]),
                             mode = "lines",
                             line = list(color = log2FC.df$PeptideColors[i], width = 10),
                             opacity = max(log2FC.df$Opacity[i], 0.1),
                             showlegend = F,
                             hoverinfo = "none")
    
      coverage.vector[(log2FC.df$Position[i] + 1):(log2FC.df$Position[i] + log2FC.df$Length[i])] <- 1
      
    }
    
  }
  
  coverage <- round((sum(coverage.vector)/length(coverage.vector)) * 100, digits = 2 )
    
  if(!missing(modified_peptide_mapping) & !missing(modifications) & !missing(data) & !missing(mod_colors)){
      
    mod_peptide_names <- unique(modified_peptide_mapping$Counterpart[which(grepl(protein_name, modified_peptide_mapping$Protein))])
      
    if(length(mod_peptide_names) != 0){
        
      mod_peptide_positions <- modified_peptide_mapping$Position[which(modified_peptide_mapping$Protein == protein_name)]
      
      mod_peptide_indices <- reshape2::melt(lapply(mod_peptide_names, function(x) which(grepl(x, modifications$counterpart.sub))))
        
      mod_peptide_expr1 <- data[modifications$index[mod_peptide_indices[, 1]], condition1]
      mod_peptide_expr2 <- data[modifications$index[mod_peptide_indices[, 1]], condition2]

      log2FC.df.mod <- data.frame(Counterpart = mod_peptide_names[mod_peptide_indices[, 2]],
                                  Sequences = modifications$sequences.sub[mod_peptide_indices[, 1]],
                                  Length = unlist(sapply(mod_peptide_names[mod_peptide_indices[, 2]], function(x) nchar(x))),
                                  Ratio = mod_peptide_expr2 - mod_peptide_expr1,
                                  Position = mod_peptide_positions[mod_peptide_indices[, 2]],
                                  Colors = rep("#ffe6ff", nrow(mod_peptide_indices)),
                                  stringsAsFactors = F)
        
      mod_type_indices <- reshape2::melt(modifications$type[mod_peptide_indices[, 1]])
        
      log2FC.df.type <- data.frame(Modification = mod_type_indices[, 1],
                                   Peptide_index = mod_type_indices[, 2],
                                   Modification_position = melt(modifications$position[mod_peptide_indices[, 1]])[, 1],
                                   Color = unlist(lapply(mod_type_indices[, 1], function(x) mod_colors$Color[which(mod_colors$Modification == x)])),
                                   stringsAsFactors = F)
        
      for(i in 1:(nrow(log2FC.df.mod))){
          
        if(is.na(log2FC.df.mod$Ratio[i])) next
          
          p <- plotly::add_trace(p,
                                 type = "scatter",
                                 x = c(log2FC.df.mod$Position[i], log2FC.df.mod$Position[i] + log2FC.df.mod$Length[i]),
                                 y = c(log2FC.df.mod$Ratio[i], log2FC.df.mod$Ratio[i]),
                                 mode = "lines",
                                 line = list(color = log2FC.df.mod$Colors[i], width = 10),
                                 showlegend = F,
                                 hoverinfo = "text",
                                 text = paste(" Sequence: ", log2FC.df.mod$Sequences[i], "<br>",
                                              "Counterpart: ", log2FC.df.mod$Counterpart[i], "<br>",
                                              "Ratio: ", log2FC.df.mod$Ratio[i], "<br>",
                                              "Position: ", log2FC.df.mod$Position[i], "<br>",
                                              "Length: ", log2FC.df.mod$Length[i], "<br>")
          )
          
          temp_type_df <- log2FC.df.type[which(log2FC.df.type$Peptide_index == i),]
          
          for(j in 1:(nrow(temp_type_df))){
            
            p <- plotly::add_trace(p,
                                   type = "scatter",
                                   x = c((log2FC.df.mod$Position[i] + temp_type_df$Modification_position[j] - 1), (log2FC.df.mod$Position[i] + temp_type_df$Modification_position[j])),
                                   y = c(log2FC.df.mod$Ratio[i], log2FC.df.mod$Ratio[i]),
                                   mode = "lines",
                                   line = list(color = temp_type_df$Color[j], width = 10),
                                   showlegend = F,
                                   hoverinfo = "text",
                                   text = paste(" Modification: ", temp_type_df$Modification[j], "<br>",
                                                "Position: ", temp_type_df$Modification_position[j])
            )
          }
      }
    }
  }
    
  m <- list(l = 20, r = 20, b = 50, t = 50, pad = 0)
  
  p <- plotly::layout(p, 
                      title = paste("VIQoR plot of protein group ", probes_n_weights$protein.group[protein_group_index], " for ", condition2, " over ", condition1, collapse = ""),
                      font = list(size = 10),
                      margin = m,
                      autosize = T,
                      xaxis = list(title = paste(paste("Quantified peptides found in both ", type.legend, ":", collapse = "", sep = ""), length(which(!is.na(log2FC.df$Ratio))) - 1, " out of ", length(peptide_names), "       Amino acid position       ", "Coverage: ", coverage, " %", sep = "", collapse = ""), size = 9),
                      yaxis = list(title = "Fold change (log2)", size = 9))

  p <- plotly::config(p, 
                      displayModeBar = T,
                      showLink = F, 
                      displaylogo = F)
  
  p$elementId <- NULL
  
  return(p)
}

# Function to perform enzymatic digestion on a single amino acid sequence.
# - Supports 8 different cleavage rules, for trypsin (considering proline or just cleavage after lysine and arginine), chymotrypsin (high and low specificity),
#   pepsin (pH 2 and pH 1.3), lysC and argC.
#   Enzyme rules according to: https://www.nature.com/articles/srep22286/tables/1
#                               https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html#exceptions
fastDigest <- function(sequence, enzyme = "trypsin"){
  
  cleavage_rules <- dplyr::case_when(
    
    enzyme == "trypsin" ~ "(?!(RP|KP))(?=(K|R))(?!(K|R)$)",
    enzyme == "trypsin.strict" ~ "(?=(K|R))(?!(K|R)$)",
    enzyme == "chymotrypsin.h" ~ "(?!(FP|YP|PY|WP))(?=(F|Y|W))(?!(F|Y|W)$)",
    enzyme == "chymotrypsin.l" ~ "(?!(FP|YP|PY|WP|LP|MP))(?=(F|Y|W|L|P))(?!(F|Y|W|L|P)$)",
    enzyme == "pepsin.2" ~ "(?=(F|L|W|Y|A|E|Q))(?!(F|L|W|Y|A|E|Q)$)",
    enzyme == "pepsin.1.3" ~ "(?=(F|L))(?!(F|L)$)",
    enzyme == "lysC" ~ "(?=(K))(?!(K)$)",
    enzyme == "argC" ~ "(?!(RP))(?=(R))(?!(R)$)"
    
  )
  
  if(!is.na(cleavage_rules)){
    
    cleave <- unlist(gregexpr(pattern = cleavage_rules, text = sequence, perl = T))

  } else {
    
    stop("Error: selected enzyme is invalid!")
    
  }
  
  if(cleave[1] == -1){
    
    return(NULL)
    
  } else {
    
    start <- c(1, cleave + 1)                          
    stop <- c(cleave, nchar(sequence))
    
    return(substring(sequence, start, stop))
    
  }
  
}

#Function that creates abundance correlation heatmaps for either protein expression (per run\condition) or peptide abundances (per run\condition).
correlation.Heatmap <- function(expr, clust_method, distance_method, n_cluster, gradient, threshold, label){
  
  if(missing(expr) | missing(threshold) | missing(label)){
    
    stop("Error: missing input argument in protein.Expression.Heatmap().")
    
  }
  
  if(missing(clust_method)){
    
    clust_method <- "Complete linkage"
    
  }
  
  if(missing(distance_method)){
    
    distance_method <- "Euclidean"
    
  }
  
  if(missing(n_cluster)){
    
    n_cluster <- "Any"
    
  }
  
  if(missing(gradient)){
    
    gradient <- "2.png"
    
  }
  
  clust_method <- dplyr::case_when(
    clust_method == "Ward's" ~ "ward.D2",
    clust_method == "Single linkage" ~ "single",
    clust_method == "Complete linkage" ~ "complete",
    clust_method == "UPGMA" ~ "average",
    clust_method == "WPGMA" ~ "mcquitty",
    TRUE ~ "complete"
  )
  
  distance_method <- dplyr::case_when(
    distance_method == "Euclidean" ~ "euclidean",
    distance_method == "Maximum" ~ "maximum",
    distance_method == "Manhattan" ~ "manhattan",
    distance_method == "Canberra" ~ "canberra",
    distance_method == "Binary" ~ "binary",
    distance_method == "Minkowski" ~ "minkowski",
    TRUE ~ "euclidean"
  )
  
  gradient <- dplyr::case_when(
    gradient == "1.png" ~ c("#00416A", "#799F0C", "#FFE000"),
    gradient == "2.png" ~ c("#C6FFDD", "#FBD786", "#f7797d"),
    gradient == "3.png" ~ c("#1E9600", "#FFF200", "#FF0000"),
    gradient == "4.png" ~ c("#77A1D3", "#79CBCA", "#E684AE"),
    TRUE ~ c("#C6FFDD", "#FBD786", "#f7797d")
  )
  
  if(n_cluster == "Any"){
    
    n_cluster <- NA
    
    if((ncol(expr) - 1) <= 2){
      
      n_cluster <- 2
      
    }
    
  } else {
    
    n_cluster <- strtoi(n_cluster)
    
  }
  
  if(label != "Peptides"){
    
    title <- paste("Protein expression correlation heatmap for weight threshold ", threshold, collapse = "")
    expr <- expr[, 2:ncol(expr)]
    
  } else {
    
    title <- paste("Peptide abundance correlation heatmap for weight threshold ", threshold, collapse = "")
    
  }
  
  remove <- which(unlist(lapply(1:ncol(expr), function(y) length(which(!is.na(expr[,y]))))) <  2)
  
  if(length(remove) > 0){
    
    expr <- data.frame(expr[, -remove])
    
    if(is.finite(n_cluster) & (n_cluster - length(remove)) > 0){
        
      n_cluster <- n_cluster - length(remove)
   
    } else {
      
      n_cluster <- 2
      
    }
    
  }
  
  if(ncol(expr) >= 2){
  
    zero.sd <- which(apply(expr, 2, sd, na.rm = TRUE) == 0)
      
    if(length(zero.sd) == 1){
        
      expr[1,zero.sd] <- expr[1,zero.sd] + 0.00005
        
    }
      
    cor.data <- cor(expr, use = "pairwise.complete.obs", method = "pearson")
    distance <- stats::dist(cor.data, method = distance_method)
    error <- try(stats::as.dendrogram(stats::hclust(distance, method = clust_method)), silent = TRUE)

    if(class(error) != "try-error"){
    
      dendrogram <- stats::as.dendrogram(stats::hclust(distance, method = clust_method))
      dendrogram <- dendextend::seriate_dendrogram(dendrogram, distance)
      cor.data <- cor.data[ ,stats::order.dendrogram(dendrogram)]
      
      p <- heatmaply::heatmaply(cor.data, 
                                hclust_method = clust_method, 
                                dist_method  = distance_method,
                                k_col = n_cluster,
                                k_row = n_cluster,
                                limits = c(-1, 1), 
                                margins = c(50, 50, 50, 10),
                                fontsize_row = 7,
                                fontsize_col = 7,
                                column_text_angle = 90,
                                dendrogram = "row",
                                xlab = label,
                                ylab = label,
                                main = title,
                                scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = gradient[1], mid = gradient[2], high = gradient[3], midpoint = 0, limits = c(-1, 1))
                                )

      p <- plotly::config(p, 
                          displayModeBar = T,
                          showLink = F, 
                          displaylogo = F)
      
      p$x$layout$font$size <- 10
      
    } else {
        
      p <- NULL
        
    }
  
  } else {
    
    p <- NULL
    
  }  
    
  return(p)
  
}