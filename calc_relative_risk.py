# NOTE, needs flatBROAD_with_risk.csv to run

import pandas as pd
import argparse

# VARIABLES
DRUG_POS = 0
GENES_POS = 3
RISK_SCORE_POS = 6

# manages inputs coming in to the python script
inputParser = argparse.ArgumentParser(description="gives the rankings of compounds by risk score associated with common genes")
inputParser.add_argument("-i", type=str, nargs=1, help="target", required=True)
inputParser.add_argument("-t", type=str, nargs=1, help="target type", required=True)

# parses the inputs into variables
args = inputParser.parse_args()

target = args.i[0]
target_type = args.t[0]

#target = "ADRB1"
#target_type = "gene"

def readCSV(filename, hasHeader):
	output = []
	
	file_in = open(filename, "r")
	filelines = file_in.readlines()
	
	isFirst = hasHeader
	for line in filelines:
		if (isFirst):
			isFirst = False
			continue
		proc_line = line.replace("\n","")
		output.append(proc_line.split(","))
	return(output)
	
broad_risk = readCSV("flatBROAD_with_risk.csv",True)
# print(len(broad_risk)) # DEBUG

gene_inputs = []
if (target_type == "drug"):
	for row in broad_risk:
		if (target in row[DRUG_POS]):
			gene_row = row[GENES_POS].replace("\'","")
			gene_row = row[GENES_POS].replace("\"","")
			associated_genes = gene_row.split("|")
			for assoc_gene in associated_genes:
				gene_inputs.append(assoc_gene)
			break
	if (len(gene_inputs) < 1):
		print("Drug target {} not found. Exiting ...".format(target))
		quit()
	print("Rankings for drug {} ...".format(target))
elif (target_type == "gene"):
	gene_inputs.append(target)
else:
	print("Target type for {} unrecognized. Exiting ...".format(target))
	quit()

raw_risk_col = "Raw Risk Score"
scaled_risk_col = "Scaled Risk Score"
drug_col = "Drugs"
# for each gene, gives ranked list of compounds ordered by raw risk score
for gene_target in gene_inputs:

	# adds space
	print("")

	drugs_risk_array = []

	for row in broad_risk:
		if (gene_target in row[GENES_POS]):
			cur_drug = (row[DRUG_POS]).replace("\"","")
			if(target_type == "drug" and target != cur_drug):
				drugs_risk_array.append([cur_drug, float(row[RISK_SCORE_POS])])
			
	analysis_array = pd.DataFrame(drugs_risk_array, columns=[drug_col,raw_risk_col])
	print("{} associated compounds found for {}".format(len(analysis_array.index), gene_target))
	print("------------------------------------")
	
	if(len(analysis_array) < 1):
		continue
	
	analysis_array[scaled_risk_col] = analysis_array[raw_risk_col]/min(analysis_array[raw_risk_col])
	
	analysis_array.sort_values([scaled_risk_col], ascending=True, inplace=True)
	#print(min(analysis_array[raw_risk_col]))
	#print(analysis_array)
	#quit()
	
	analysis_array.reset_index()
	count = 1
	for index, row in analysis_array.iterrows():
		drug = row[0]
		raw = row[1]
		scaled = row[2]
		#drug = analysis_array.loc[count - 1][0]
		#raw = analysis_array.loc[count - 1][1]
		#scaled = analysis_array.loc[count - 1][2]
		print("Rank {}: {} with relative risk of {}".format(count, drug, scaled))
		count = count + 1
	
	
