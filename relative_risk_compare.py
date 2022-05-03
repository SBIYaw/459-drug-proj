from math import comb
# note: requires python 3.8+

# uses hypergeometric distribution p-value calculation on overlaps
def processOverlap(genes_A, genes_B):
	#print("Overlaps: {}, {}".format(genes_A, genes_B))

	genes_A = genes_A.replace("\"", "")
	genes_B = genes_B.replace("\"", "")
	genes_A_array = genes_A.split("|")
	genes_B_array = genes_B.split("|")

	count_overlap = 0
	on_AEffect = 0
	on_BEffect = 0
	#print(genes_A_array)
	#print(genes_B_array)
	if (genes_A_array[0] == "" or genes_B_array[0] == ""):
		return(0,0,0)

	for A_gene in genes_A_array:
		for B_gene in genes_B_array:
			if (A_gene == B_gene):
				count_overlap = count_overlap + 1

	#on_AEffect = (count_overlap / (len(genes_A_array) + len(genes_B_array)  - count_overlap)) #+ (.1 * count_overlap)
	#on_BEffect = (count_overlap / len(genes_B_array)) #+ (.1 * count_overlap)

	# n is genes in A, B is genes in B, b is count_overlap, and N is either total universe (20,000+ human genes or 654 from our entire gene subset)
	#N = 654

	# hg19 total
	N = 23271
	on_AEffect = (comb(len(genes_B_array), count_overlap) * comb(N - len(genes_B_array), len(genes_A_array) - count_overlap))/ (comb(N, len(genes_A_array)))
	on_BEffect = (comb(len(genes_A_array), count_overlap) * comb(N - len(genes_A_array), len(genes_B_array) - count_overlap))/ (comb(N, len(genes_B_array)))

	return(count_overlap, on_AEffect, on_BEffect)

def readTable(filename):
	file_var = open(filename, "r")
	file_lines = file_var.readlines()[1:]

	return(file_lines)

def writeTable(filename, dataset, header):
	file_var = open(filename, "w")

	header=",".join(header)
	file_var.write(header)
	file_var.write("\n")
	file_var.flush()

	for line in dataset:
		file_var.write(line)
		file_var.write("\n")
	file_var.flush()
	file_var.close()

def processGeneNames(geneset_names):
	output = []
	for gene_list in geneset_names:
		gene_list2 = gene_list.split("|")
		for gene in gene_list2:
			output.append(gene)
	output = set(output)
	return(list(output))


def processLine(line):
	orig_line = line.replace("\n","")
	orig_line = orig_line.replace("\"","")
	line_fields = orig_line.split(",")

	genes = line_fields[3]
	risk_score = float(line_fields[6])

	return(orig_line, genes, risk_score)


broad_risk = readTable("flatBROAD_with_risk.csv")
ddr_out = []
ddr = []

countA = 0
for line_A in broad_risk:
	data_line_A, genes_A, risk_score_a = processLine(line_A)
	countB = 0
	for line_B in broad_risk:
		if (countA == countB):
			countB = countB + 1
			continue
		data_line_B, genes_B, risk_score_b = processLine(line_B)

		overlaps, on_AEffect, on_BEffect = processOverlap(genes_A, genes_B)
		#print(overlaps)
		# on_AEffect == on_BEffect

		if (overlaps > 0):
			added_data = str(on_AEffect * risk_score_b) + "," + str(overlaps) # + "," + str(on_BEffect)
			new_line = data_line_A + "," + added_data + "," + data_line_B
			new_line.replace("\"","")
			ddr_out.append(new_line)
			ddr.append(new_line.split(","))
		countB = countB + 1
		#quit()
	countA = countA + 1

header_row = ["drug_a", "clinical_phase_a","moa_a","target_a","disease_area_a","indication_a","risk_score_a", "scaled_risk_score", "total_overlaps","drug_b", "clinical_phase_b","moa_b","target_b","disease_area_b","indication_b","risk_score_b"]

writeTable("DDR_database_risk_fixed.csv", ddr_out, header_row)
print("Wrote out DDR database!")

# flatten database into simple network - done
ddr_net = []
ddr_net_s = []

for data_row in ddr:
	#print(data_row)
	target_a = data_row[3]
	target_a = target_a.replace("\"", "")
	targets_a = target_a.split("|")
	#print(targets_a)

	target_b = data_row[12]
	target_b = target_b.replace("\"", "")
	targets_b = target_b.split("|")

	#print(targets_b)

	overlaping_targets = []
	for gene_a in targets_a:
		for gene_b in targets_b:
			if (gene_a == gene_b):
				overlaping_targets.append(gene_a)
	#print(overlaping_targets)
	for overlap in overlaping_targets:
		new_row = ",".join(data_row)
		new_row = new_row + "," + str(overlap)
		ddr_net.append(new_row)

		ddr_net_s.append("{},{},{},{}".format(data_row[0],data_row[7],overlap,data_row[9]))


new_header = header_row
new_header.append("gene")
writeTable("DDR_full_net_risk_fixed.csv", ddr_net, header_row)
print("Wrote out DDR full network!")

short_header = ["drug_a","scaled_risk_score","gene","drug_b"]
writeTable("DDR_simple_net_risk_fixed.csv", ddr_net_s, short_header)
print("Wrote out DDR simple network!")

print("Done!")
quit()

# IGNORE BELOW

# filter by disease area
DISEASE_AREA_A_POS = 4
DISEASE_AREA_B_POS = 13

filtered_ddr_out = []

for row_a in ddr:
	da_a = row_a[DISEASE_AREA_A_POS]
	da_b = row_a[DISEASE_AREA_B_POS]

	da_b = da_b.split("|")
	if (da_a in da_b):
		new_row = ",".join(row_a)
		filtered_ddr_out.append(new_row)

writeTable("DDR_filtered_by_disease_area.csv", filtered_ddr_out, header_row)
