library(igraph)

# note that this uses the tsv which doesn't have the header at the top 
broad_database <- read.csv("broad_repurp_data.tsv",header=TRUE, sep="\t")
tot_len <- length(broad_database[,1])

getTargetList <- function(set_string) {
  val <- strsplit(set_string,"\\|")[[1]]
  return(val)
}

# returns number of overlaps between source and other
numOverlaps <- function(source_vector, other_vector) {
  count_overlap <- 0
  all_total <- c(source_vector, other_vector)
  count_overlap <- length(all_total) - length(unique(all_total))
  return(count_overlap)
}


# for making source_drug_target_gene list for cytoscape
cyto_net <- vector()

i = 1
while (i <= tot_len) {
  compound <- broad_database$pert_iname[i]
  phase <- broad_database$clinical_phase[i]
  meth_of_action <- broad_database$moa[i]
  preproc_target <- broad_database$target[i]
  disease_area <- broad_database$disease_area[i]
  indication <- broad_database$indication[i]
  
  target_list <- getTargetList(preproc_target)
  if (length(target_list) < 1) {
    i <- i + 1
    next
  }
  for (targets in target_list) {
    cyto_net <- append(cyto_net, paste(compound,"drug",targets,"gene", sep="\t"))
  }
  
  i <- i + 1
}
write(cyto_net, "BROAD_drug_to_gene.tsv")

# -----------------------------------

# BELOW: makes a drug to drug network based on the number of gene-products shared
# between two compounds
# also has the basis for an igraph, which has a different separator

k = 1
drug_to_drug <- vector()
drug_to_drug_igraph <- vector()
vals_set <- c()
THRESHOLD <- 5
while (k <= tot_len) {
  compound1 <- broad_database$pert_iname[k]
  disease_area1 <- broad_database$disease_area[k]
  cur_list <- broad_database$target[k]
  
  # genes saved as a single list string, split by getTargetList on all rows
  target_lists <- lapply(broad_database$target, getTargetList)
  j = 1
  for (target_list in target_lists) {
    if (length(target_list) < 1) {
      j <- j + 1
      next
    }
    val <- numOverlaps(getTargetList(cur_list), target_list)
    if (val > THRESHOLD) {
      vals_set <- append(vals_set, val)
      compound2 <- broad_database$pert_iname[j]
      disease_area2 <- broad_database$disease_area[j]
      if (disease_area2 == "") {
        disease_area2 <- "None"
      }
      drug_to_drug <- append(drug_to_drug, paste(compound1,disease_area1,val, compound2,disease_area2, sep="\t"))
      drug_to_drug_igraph <- append(drug_to_drug_igraph, paste(compound1, compound2, sep=" "))
    }
    j <- j + 1
  }
  k <- k + 1
}
write(drug_to_drug,"BROAD_drug_to_drug_weighted_net_thr5.tsv")


# roundabout way of getting data into igraph
write(drug_to_drug_igraph,"temp_drug_to_drug.txt")
drug_to_drug_table <- read.table("temp_drug_to_drug.txt")
BROAD_graph <- igraph::graph_from_data_frame(drug_to_drug_table)
# still thinking about what to do with igraph analysis

