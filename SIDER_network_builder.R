### for making the SIDER network

# Takes about 30 minutes on Intel i7-10750H, 2.6Ghz, 16GB (though capped to one core)

ids_to_sider <- read.csv("meddra_freq.tsv",header=FALSE, sep="\t")
ids_to_name <- read.csv("drug_names.tsv",header=FALSE, sep="\t")

SIDER_network <- vector()

i <- 1
while (i <= dim(ids_to_sider)[1]) {
  
  cur_row <- ids_to_sider[i,]
  
  current_drug <- ""
  current_sideff <- ""
  current_weight <- 0
  # only taking preferred term, exact values, and not factoring placebo
  # comparison
  min_fr <- as.double(cur_row$V6)
  max_fr <- as.double(cur_row$V7)
  if (min_fr == max_fr) {
    current_weight <- max_fr
    if (cur_row$V8 == "PT") {
      current_sideff <- cur_row$V10
      
      j <- 1
      while (j <= dim(ids_to_name)[1]) {
        drug_name_row <- ids_to_name[j,]
        if (cur_row$V1 == drug_name_row$V1) {
          current_drug <- drug_name_row$V2
          
          new_row <- paste(current_drug,current_sideff,current_weight,sep="\t")
          SIDER_network <- append(SIDER_network, new_row)
          break
        }
        j <- j + 1
      }
    }
  }
  i <- i + 1
}

length(SIDER_network)
write(SIDER_network,"SIDER_drug_to_sideff_weighted_nothres.tsv")


# -----------------------------------

# reads in SIDER network data as table and sets thresholds based on risk percentage

split_SIDER <- read.csv("SIDER_drug_to_sideff_weighted_nothres.tsv", sep="\t",header=FALSE)

thres_five <- .05
SIDER_thres_five <- split_SIDER[as.double(split_SIDER$V3) > .05,]
thres_ten <- .10
SIDER_thres_ten <- split_SIDER[as.double(split_SIDER$V3) > thres_ten,]
thres_twenty <- .20
SIDER_thres_twenty <- split_SIDER[as.double(split_SIDER$V3) > thres_twenty,]

out_array <- vector()
j <- 1
while (j <= dim(SIDER_thres_twenty)[1]) {
  out_array <- append(out_array, paste(SIDER_thres_twenty$V1[j],"drug",SIDER_thres_twenty$V2[j],"se",SIDER_thres_twenty$V3[j],sep="\t"))
  j <- j + 1
}

write(out_array,"SIDER_drug_to_sideff_weighted_20pthres.tsv")
