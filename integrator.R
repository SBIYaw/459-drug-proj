# Integrates BROAD, SIDER and ADR

# Reading ADR
adr_database <- read.csv("ADR_paper_supplement_4_rankings.csv")
z_scores <- adr_database$rank_score / adr_database$rank_stdev_out_of_2929

# looks into values
#val<- lm(1:length(z_scores) ~ z_scores)
#plot(val)

overlap_set <- broad_database$pert_iname[broad_database$pert_iname %in% split_SIDER$V1]
ov_SIDER <- split_SIDER[split_SIDER$V1 %in% overlap_set,]
colnames(ov_SIDER) <- c("drug","side_effect","se_frequency")

drug_riskscores <- read.csv("drug_riskscores_2.csv")
ov_risk <- drug_riskscores[drug_riskscores$drug %in% overlap_set,]


# NOTE: issue with ext_broad seems to be causing issue with BROAD_with_risk.csv
#
# unified_drug_risk made using the BROAD_database extended into multiple rows
## per compound

#ov_BROAD <- ext_broad[ext_broad$drug %in% overlap_set, ]
#final_overlaps <- ov_risk$drug[ov_risk$drug %in% overlap_set]
#ov_risk <- drug_riskscores[drug_riskscores$drug %in% final_overlaps,]
#ov_BROAD <- ov_BROAD[ov_BROAD$drug %in% final_overlaps, ]
#
## maps risk scores to each row by repeating the rows in ov_risk according to the
### count of each drug in the extended BROAD table
#
#count_table <- table(ov_BROAD$drug)
#rep_vector <- c()
#j <- 1
#while(j <= length(names(count_table))){
#  rep_vector <- append(rep_vector, count_table[j])
#  j <- j + 1
#}
#rep_vector <- unname(rep_vector)
#matching_riskscores <- ov_risk[rep(seq_len(nrow(ov_risk)), rep_vector), ]
#
#unified_drug_risk <- cbind(ov_BROAD,matching_riskscores)
#unified_drug_risk <- unified_drug_risk[, !duplicated(colnames(unified_drug_risk))]
#
## for network viewing
#write.csv(unified_drug_risk, file="BROAD_with_risk.csv",row.names=FALSE)

# flat_drug_risk binds risk scores to BROAD_database and keeps multiple genes in
## each row (can be readily checked with "in" in python and "grepl" in R)

ov_BROAD2 <- broad_database[broad_database$pert_iname %in% overlap_set,]

final_overlaps2 <- ov_risk$drug[ov_risk$drug %in% ov_BROAD2$pert_iname]
ov_BROAD2 <- ov_BROAD2[ov_BROAD2$pert_iname %in% final_overlaps2,]
ov_BROAD2 <- ov_BROAD2[order(ov_BROAD2$pert_iname),]
ov_risk2 <- drug_riskscores[drug_riskscores$drug %in% final_overlaps2,]

# for comparisons
colnames(ov_BROAD2) <- c("drug","clinical_phase","moa","target","disease_area","indication")
flat_drug_risk <- cbind(ov_BROAD2, ov_risk2)
flat_drug_risk <- flat_drug_risk[, !duplicated(colnames(flat_drug_risk))]
# for processing overlaps
write.csv(flat_drug_risk, file="flatBROAD_with_risk.csv",row.names=FALSE)
