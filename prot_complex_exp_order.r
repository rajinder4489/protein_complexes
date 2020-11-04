rm(list=ls())

#libraries
library(dbnR)
library(igraph)
library(reshape2)
library(ggplot2)
library(gridExtra)

#set path
setwd("C://Users/r.gupta/Desktop/prot_complex_new/")


###Functions start###

#Takes as input the complexes where atleat for one of the subunits the stiochiometry os 0 and changes it to 1
#Returns the complex as received with changed stoichiometry

change_stoi <- function(cmplx)
{
	changed_stoi <- gsub("\\(0\\)", "(1)", cmplx)
	return(changed_stoi)
}

#Takes as input one complex at a time and complex file
#Extracts the sub-units and stoichiometry from the complex
#Calls itself when a sub-unit is another complex and resolves the stoichiometry for the complex in complex by multiplying their stoichiometry. As CPX-1: X1(2), X2(2), CPX-2(2) and CPX-2: X3(3), X4(3). The resulting stoichiometry for CPX-1 will be X1(2), X2(2), X3(3*2), X4(3*2)
#Removes complexes with small molecules and nucleic acids
#Returns list of sub-units and their stoichiometry

subunits_stoichiometry <- function(cmplx, all_cmplx)
{
  ret <- list()
  
  cmplx[2] <- change_stoi(cmplx[2])
  
  if(!grepl("\\(0\\)|CHEBI|URS", cmplx[2]))
  {
    ret$status <- 1     #setting the status to 1
    
    sub_stoi <- unlist(strsplit(as.vector(unlist(cmplx[2])), "\\|"))
    ret$sub <- gsub("\\(.\\)", "", sub_stoi)  #subunits to return
    
    stoi <- gsub(paste(ret$sub, collapse = "|"), "", sub_stoi)
    ret$stoi <- as.numeric(gsub("\\(|\\)", "", stoi))   #stoichiometry to return
    
    for(j in 1:length(ret$sub))
    {
      subunits <- ret$sub[j]
      if(grepl("CPX", subunits))
      {
        sub_complex <-  all_cmplx[all_cmplx$Complex.ac == subunits, ]  #when a subunit is also a complex
        sub_ss <- subunits_stoichiometry(sub_complex, all_cmplx)  #calling the function from within the function
        
        if(sub_ss$status == 1)  #status returned when computing for the complex which is a subunit 
        {
          ret$sub <- c(ret$sub, sub_ss$sub)
          ret$stoi <- c(ret$stoi, sub_ss$stoi*ret$stoi[j])  #stoichiometry of the subunits from the complex which is a subunit in another complex is given by the stoichiometry of the subunits in the complex times stoichiometry of the complex in the main complex
        }
        else
        {
          ret$status <- 0   #because for one the subunits of the complex which is a subunit in another complex did not pass the filters
        }
      }
    }
  }
  else
  {
    ret$status <- 0
  }
  return(ret)
}

#Takes as input the subunits of a complex (one complex at a time) and biomart_file
#Returns the Ensembl Transcript ids for the Uniprot identifiers of the subunits and also the count of the subunits that were mapped successfully.


mapping_to_ensembl <- function(subu, rel)	#subu: subunit
{
  ret <- list()
  ret$status = 1
  
  ret$mapping <- rel[rel$UniProtKB.Swiss.Prot.ID %in% subu, ]
  
  if(length(unique(ret$mapping$UniProtKB.Swiss.Prot.ID)) != length(subu))
  {
    ret$status = 0    #when all sub-units are not mapped to transcript ids
  }
  return(ret)
}


#Takes as input Ensembl gene ids for a complex (one complex at a time) and APPRIS_file
#When one Uniprot id is mapped to multiple Ensembl ids, they are returned as a group 
#Returns the principal transcript for the ensembl genes

select_transcripts <- function(map, rel)
{
  ret <- list()
  ret$status <- 1
  uniprot_ids <- unique(map$UniProtKB.Swiss.Prot.ID)
  no_uniprot_ids <- length(unique(map$UniProtKB.Swiss.Prot.ID))
  
  #counting the subunits from APPRIS and longest transcript
  appris_count <- 0
  longest_count <- 0
  
  appris_order_of_selection <- c("principal1", "principal2", "principal3", "principal4", "alternative1", "alternative2")
  
  ptr <- map[FALSE,]   #principal transcripts
  
  for(appris in appris_order_of_selection)
  {
    if(nrow(map)>0)
    {
      selection <- map[map$APPRIS.annotation == appris, ]
      
      appris_count <- appris_count + nrow(selection)
      
      if(nrow(selection)>0)
      {
        ptr <- rbind(ptr, selection)
      
        map <- map[!grepl(paste(unique(selection$Gene.stable.ID), collapse="|"), map$Gene.stable.ID),]    #to remove the gene ids for which transcript has been selected
      }
    }
  }
  
  if(no_uniprot_ids != length(unique(ptr$UniProtKB.Swiss.Prot.ID)))
  {
    if(sum(!uniprot_ids %in% ptr$UniProtKB.Swiss.Prot.ID) > 0)
    {
      missing_ids <- uniprot_ids[!uniprot_ids %in% ptr$UniProtKB.Swiss.Prot.ID]
      
      longest_count <- longest_count + length(missing_ids)
      
      for(id in missing_ids)
      {
        selecting_on_length <- rel[rel$UniProtKB.Swiss.Prot.ID %in% id, ]
        longest_len <- max(selecting_on_length$Transcript.length..including.UTRs.and.CDS.)
        
        selected <- selecting_on_length[selecting_on_length$Transcript.length..including.UTRs.and.CDS. %in% longest_len, ]
        ptr <- rbind(ptr, selected)
      }
      
      if(no_uniprot_ids != length(unique(ptr$UniProtKB.Swiss.Prot.ID)))
      {
        ret$status <- 0
      }
    }
  }
  ret$ptr <- unique(ptr)
  ret[["count"]][["appris"]] <- appris_count
  ret[["count"]][["longest"]] <- longest_count
  
  return(ret)
}


#Takes principal transcripts for a complex (one subunit of one complex at a time) and time-series expression data
#Where multiple principal transcripts are present (due to one Uniprot id mapped to many Ensembl gene ids), the expression for the subunit is calculated by summing the expression for all transcripts
#Returns the expression of the subunits

subunits_exp <- function(ptr, expr, rel)
{
  ret <- list()
  ret$status <- 1
  pt_exp <- expr[expr$tids %in% unique(ptr$Transcript.stable.ID), ]
  
  pt_exp_uniprot <- merge(pt_exp, rel[,c("Transcript.stable.ID", "UniProtKB.Swiss.Prot.ID")], by.x= "tids", by.y="Transcript.stable.ID")
  
  ret$pt_exp <- pt_exp_uniprot
  
  if(nrow(pt_exp) > 0)
  {
    sub_exp <- aggregate(pt_exp_uniprot[, -which(colnames(pt_exp_uniprot) %in% c("tids", "UniProtKB.Swiss.Prot.ID"))], by=list(subunits=pt_exp_uniprot$UniProtKB.Swiss.Prot.ID), FUN=sum)
    ret$sub_exp <- sub_exp
  }
  else
  {
    ret$status = 0
  }
  
  if(length(unique(ptr$Transcript.stable.ID)) != nrow(pt_exp))
  {
    ret$status = 0
  }
  
  return(ret)
}


#Takes as input the expression of the sub-units and their stoichiometry (one complex at a time)
#Returns the expression of the complex

prot_complex_exp <- function(subs_exp, subunits, stoichiometry, tp)
{
  ret <- list()
  
  dff <- as.data.frame(cbind(subunits, as.integer(stoichiometry)))
  colnames(dff) <- c("subunits", "stoichiometry")
  
  exp_stoic <- merge(subs_exp, dff, by="subunits")
  
  for_norm <- exp_stoic[, -which(colnames(exp_stoic) %in% c("subunits", "stoichiometry"))]
  exp_normalized <- for_norm/c(as.integer(exp_stoic$stoichiometry))  #normalized for stoichiometry
  
  tp_exp <- data.frame(matrix(NA, nrow = nrow(exp_normalized), ncol = length(tp)+1))
  colnames(tp_exp) <- c("subunits", tp)
  
  tp_exp[,1] <- exp_stoic$subunits
  
  for(x in tp)
  {
    tp_exp[, x] <- rowMeans(exp_normalized[, grepl(x, colnames(exp_normalized))])
  }
  
  pc_exp <- apply(tp_exp[, -which(colnames(tp_exp) %in% "subunits")], 2, FUN=min)
  
  ret$pt_exp <- exp_stoic
  ret$pt_exp_norm <- exp_normalized
  ret$tp_exp <- tp_exp
  ret$pc_exp <- pc_exp
  
  return(ret)
}


#Takes as input the expression of the subunits

#1
#Using dbnR
#Returns dbn model and plot

network_dbnR <- function(subs_exp, no_tp)
{
	ret <- list()

	dff <- as.data.frame(t(subs_exp))
	colnames(dff) <- dff[1,]
	dff <- dff[-1,]
	dff <- as.data.frame(sapply(dff, as.double))

#  dbn <- learn_dbn_struc(dff, size=no_tp)
  
	out <- tryCatch(
	{
		dbn <- learn_dbn_struc(dff, size=no_tp)
		dbn_plot <- plot_dynamic_network(dbn, offset = 1)
		ret$dbn <- dbn
		ret$dbn_plot <- dbn_plot
		return(ret)
	},
	error=function(cond) {
		message("learn_dbn_struc returned errors")
		message("Here's the original error message:")
		message(cond)
		return(0)
	},
		warning=function(cond) {
		message("learn_dbn_struc returned warning")
		message("Here's the original warning message:")
		message(cond)
		return(0)
	},
	finally={
#		print("Here")
	}
	)
	return(out)
}


#2
#bnstruct

netwrok_bnstruct <- function()
{
  
}


#3
#bnlearn
netwrok_bnlearn <- function()
    
#Takes as input the adjacency matrix
#Selects the highest probability edges for each time point
#Uses highest probability edges for generate the assembly order
#Returns the order of assembly and images

finding_order <- function()
{}

###Functions end###


###Main starts###

#Reading files
complex_file <- read.delim(file="complexPortal_homo_sapiens.tsv", sep="\t", stringsAsFactors=F)

biomart_file <- read.delim(file="biomart.txt", sep="\t", stringsAsFactors = F)    #also has details on APPRIS and transcript length

exp_file <- read.delim(file="ConUNTR_isoform_fpkm_liver_hecatos.txt", sep="\t", stringsAsFactors = F)

timepoints <- c("000", "002", "008", "024", "072", "168", "240", "336") #time points
timepoints_name <- c("0 hrs", "2 hrs", "8 hrs", "24 hrs", "72 hrs", "168 hrs", "240 hrs", "336 hrs") #time points for naming in graphs
#Processing starts

counter1 = counter2 = counter3 = counter4 = counter5 = counter6 = counter7 = 0
complex_expression <- complex_expression2 <- complex_network <- c()  #to make a list of complexes for which networks are built

all_results <- list()

for(i in 1:nrow(complex_file))
{
  complex_name <- complex_file[i,1]
  
  ss <- subunits_stoichiometry(complex_file[i,c(1,5)], complex_file[,c(1,5)])  #returned subunits and stoichiometry (ss). If the complex did not pass the filters, NULL is returned
  
  if(ss$status == 1)   #When all sub-units passed the filters; map the Uniprot ids to transcript ids; counter1
  {
    #remove the complexes from the list of subunits and their corresponding stoichiometry
    
    ss_filtered <- ss
    
    #remove complex ids, if any 
    for(j in 1: length(ss$sub))
    {
      subunits <- ss$sub[j]
      if(grepl("CPX", subunits))
      {
        ss_filtered$sub <- ss_filtered$sub[-which(grepl("CPX", subunits))]
        ss_filtered$stoi <- ss_filtered$stoi[-which(grepl("CPX", subunits))]
      }
    }
    
    all_results[[complex_name]][["ss_filtered"]] <- ss_filtered
    
    mapping <- mapping_to_ensembl(ss_filtered$sub, biomart_file)
    
    all_results[[complex_name]][["mapping"]] <- mapping
    
    if(mapping$status == 1)		#counter2
    {
      principal_transcript <- select_transcripts(mapping$mapping, biomart_file)
      
      all_results[[complex_name]][["principal_transcript"]] <- principal_transcript
      
      if(principal_transcript$status == 1)		#counter3
      {
        sub_exp <- subunits_exp(principal_transcript$ptr, exp_file, biomart_file)
        
        all_results[[complex_name]][["sub_exp"]] <- sub_exp
        
        
        
        if(sub_exp$status == 1)		#counter4
        {
          #changing Uniprot ids to Names
          sub_exp_names <- merge(sub_exp$sub_exp, unique(biomart_file[, c("UniProtKB.Swiss.Prot.ID", "UniProtKB.Gene.Name.symbol")]), by.x = "subunits", by.y = "UniProtKB.Swiss.Prot.ID")
          sub_exp_names <- sub_exp_names[, -which(colnames(sub_exp_names) %in% "subunits")]
          sub_exp_names <- sub_exp_names[,c(ncol(sub_exp_names),1:(ncol(sub_exp_names)-1))]
          
          all_results[[complex_name]][["sub_exp"]][["sub_exp_names"]] <-  sub_exp_names
          
          #finding the protein complex expression
          
          complex_expression <- c(complex_expression, complex_name)
          pc_exp <- prot_complex_exp(sub_exp$sub_exp, ss_filtered$sub, ss_filtered$stoi, timepoints)
          
          all_results[[complex_name]][["pc_exp"]] <- pc_exp
          
          #Generating the dbn
          if(nrow(as.data.frame(sub_exp$sub_exp))>1 & 
             sum(rowSums(sub_exp$sub_exp==0)/(ncol(sub_exp$sub_exp)-1)*100 < 50) == nrow(sub_exp$sub_exp))  #if a heteromer, only then generate the network; counter 5
          {
            complex_expression2 <- c(complex_expression2, complex_name)#complexes for which expression calculation makes sense
            
            network <- network_dbnR(sub_exp_names, length(timepoints))
            
            if(is.list(network))
            {
              complex_network <- c(complex_network, complex_name) #complexes for which netwroks were built
              all_results[[complex_name]][["dbnR"]] <- network
              counter7 = counter7 + 1 #number of complexes for which dbns formed
            }
           else
            {
              counter6 = counter6 + 1 # number of complexes for which dbn not formed, because no connection could be predicted
            }
          }
          else
          {
            counter5 = counter5 + 1 #when a homomer or no expression in more than 50% cases
          }
        }
        else
        {
          counter4 = counter4 + 1
        }
      }
      else
      {
        counter3 = counter3 + 1
      }
    }
    else
    {
      counter2 = counter2 + 1 #counts number of complexes discarded in step 2
    }
  }
  else
  {
    counter1 = counter1 + 1 #counts number of complexes discarded in step 1
  }
}
#Processing ends

print(paste(counter1, counter2, counter3, counter4, counter5, counter6, counter7))
print(complex_network)

#making expression plots
for(cmplx in complex_expression)
#for(cmplx in ids_sel)
{
  naam <- gsub("/", " ", complex_file[complex_file$Complex.ac %in% cmplx, 2])
  
  sub_expr <- all_results[[cmplx]][["sub_exp"]][["sub_exp_names"]]
  
  sub_expr_avg <- data.frame(matrix(NA, nrow = nrow(sub_expr), ncol = length(timepoints)+1))
  colnames(sub_expr_avg) <- c("subunits", timepoints)
  sub_expr_avg[,1] <- sub_expr$UniProtKB.Gene.Name.symbol
  
  for(x in timepoints)
  {
    sub_expr_avg[, x] <- rowMeans(sub_expr[, grepl(x, colnames(sub_expr))])
  }
  
  pc_expr <- c(cmplx, all_results[[cmplx]][["pc_exp"]][["pc_exp"]])
  
  all_expr <- rbind(sub_expr_avg, pc_expr)
  colnames(all_expr) <- c("subunits", timepoints_name)
#  rownames(all_expr) <- all_expr$subunits
#  all_expr <- all_expr[, -which(colnames(all_expr) %in% "subunits")]
  all_expr$type <- c(rep("Subunit", nrow(all_expr)-1), "Protein complex")
  
  melted <- melt(as.data.frame(all_expr), id.vars = c("subunits", "type"))
  melted[,"value"] <- ifelse(melted[, "value"] == "0", "1", melted[,"value"])  
  melted[,"value"] <- as.numeric(melted[,"value"])		#if to take log  melted[,"value"] <- log(as.numeric(melted[,"value"]))
  
  melted$order <- factor(unique(melted$subunits), levels = c(sort(unique(melted$subunits[-which(melted$subunits %in% cmplx)])), cmplx))
  
  
  p <- ggplot(data=melted, aes(x=order, y=value, fill=type)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "bottom") +
#    theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") +
    facet_grid(. ~ variable , space="free_x", scales="free_x", switch="x") +
#    labs(title = "", x ="" , y="log(FPKM)", fill = "Type")
    labs(title = "", x ="" , y="FPKM", fill = "Type")
  
  p_grid <- grid.arrange(p, layout_matrix = rbind(c(1)))
  ggsave(file=paste("results/", cmplx, "_", naam, "_expr", ".jpeg", sep=""), plot = p, device = NULL, path = NULL, scale = 1, width = 10, height = 4, units = "in", dpi = 200, limitsize = TRUE) 
}

#Making network images
counterx <- 0

for(cmplx in complex_network)
{
  if (length(all_results[[cmplx]][["ss_filtered"]][["sub"]]) > 2)
  {
    counterx <- counterx + 1
    
    naam <- gsub("/", " ", complex_file[complex_file$Complex.ac %in% cmplx, 2])
    mat <- as.matrix(unlist(all_results[[cmplx]][["dbnR"]][["dbn_plot"]]))   #matrix obtained from dbn plot
    
    edges_nodes <- mat[grep("edges.to|edges.from", rownames(mat)), ]
    
    melted_en <- melt(edges_nodes)
    
    dff <- cbind(melted_en[ c(1:(nrow(melted_en)/2)),], melted_en[ c(((nrow(melted_en)/2)+1):nrow(melted_en)),])
    colnames(dff) <- c("From", "To")
    
    adj_mat <- get.adjacency(graph.edgelist(as.matrix(dff), directed=TRUE))
    
    graph_adj_mat <- graph_from_adjacency_matrix(adj_mat, mode = "directed", weighted = NULL,
                                diag = TRUE, add.colnames = NULL, add.rownames = NA)
    jpeg(filename = paste("images/", cmplx, "_", naam, "_net-dbnR", ".jpeg", sep=""))
    plot(graph_adj_mat)
    dev.off()
	}
}

aa = bb = cc = 0
for(cmplx in complex_expression)
{
  appris <- all_results[[cmplx]]$principal_transcript$count$appris
  longest <- all_results[[cmplx]]$principal_transcript$count$longest
  if(appris > 0 & longest == 0)
  {
   aa = aa+1 
  }
  else if(longest > 0 & appris == 0)
  {
    bb = bb+1
  }
  else
  {
    cc = cc+1
  }
}
print(paste(aa, bb, cc, sep=" "))

#bnstruct
#run the other program named bnstruct.r  


###Main ends###
