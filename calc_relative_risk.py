# NOTE, needs flatBROAD_with_risk.csv to run

import pandas as pd
import argparse
import sys
from val import *

# VARIABLES
DRUG_POS = 0
GENES_POS = 3
RISK_SCORE_POS = 6

# manages inputs coming in to the python script
'''
inputParser = argparse.ArgumentParser(description="gives the rankings of compounds by risk score associated with common genes")
inputParser.add_argument("-i", type=str, nargs=1, help="target", required=True)
inputParser.add_argument("-t", type=str, nargs=1, help="target type", required=True)

# parses the inputs into variables
args = inputParser.parse_args()

target = args.i[0]
target_type = args.t[0]
'''
#target = "ADRB1"
#target_type = "gene"

def main(target, target_type):
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
				
				orig_risk = float(row[RISK_SCORE_POS])
				for assoc_gene in associated_genes:
					gene_inputs.append(assoc_gene)
				break
		if (len(gene_inputs) < 1):
			print("Drug target {} not found. Exiting ...".format(target))
			sys.exit()
		print("Rankings for drug {} ...".format(target))
	elif (target_type == "gene"):
		gene_inputs.append(target)
	else:
		print("Target type for {} unrecognized. Exiting ...".format(target))
		sys.exit()

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

		analysis_array[scaled_risk_col] = round(100 * analysis_array[raw_risk_col]/orig_risk, 2)
		
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

			# validation for drug targets
			if (target_type == "drug"):
				approved_diseases = validate(target, drug)
				if len(approved_diseases) == 0:
					print("Rank {}: {} with relative risk of {}% of {} | Validation: {}".format(count, drug, scaled, target, 'FP'))
				else:
					diseases = ','.join([str(disease) for disease in approved_diseases])
					print("Rank {}: {} with relative risk of {}% of {} | Validation: {}, can be repurposed for diseases: {}".format(count, drug, scaled, target, 'TP', diseases))
				
			else:
				print("Rank {}: {} with relative risk of {}% of {}".format(count, drug, scaled, target))
			count = count + 1
		

if __name__ == "__main__":
	inputParser = argparse.ArgumentParser(description="gives the rankings of compounds by risk score associated with common genes")
	inputParser.add_argument("-i", type=str, nargs=1, help="target", required=True)
	inputParser.add_argument("-t", type=str, nargs=1, help="target type", required=True)

	# parses the inputs into variables
	args = inputParser.parse_args()

	target = args.i[0]
	target_type = args.t[0]
    
	main(target, target_type)
