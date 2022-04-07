# NOTE, needs flatBROAD_with_risk.csv to run
# NOTE, nonfunctional due to numpy array issue, still working on it

import numpy as np
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
			associated_genes = row[GENES_POS].split("|")
			for assoc_gene in associated_genes:
				gene_inputs.append(assoc_gene)
			break
	if (len(gene_inputs)):
		print("Drug target {} not found. Exiting ...".format(target))
		quit()
		
elif (target_type == "gene"):
	gene_inputs.append(target)
else:
	print("Target type for {} unrecognized. Exiting ...".format(target))
	quit()


# for each gene, gives ranked list of compounds ordered by raw risk score
for gene_target in gene_inputs:
	count=0

	# --- DEBUG: cause of issue --- #
	analysis_array = np.array() #?
	quit()
	# --- DEBUG: cause of issue --- #
	
	for row in broad_risk:
		if (gene_target in row[GENES_POS]):
			analysis_array <- np.append(analysis_array, (row[DRUG_POS], row[RISK_SCORE_POS]))
			
	analysis_array = analysis_array[np.argsort(analysis_array[:,2])]
	print("{} associated compounds found for {}".format(len(analysis_array), gene_target))
	
	count = 1
	for output_row in analysis_array:
		print("Rank {}: {} - {}".format(count, output_row[0], output_row[1]))
		count = count + 1

	
