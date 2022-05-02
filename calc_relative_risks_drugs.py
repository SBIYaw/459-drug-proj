import pandas as pd
import argparse
import sys
from val import *

DRUG_B_POS = 9
LOG_REL_RISK_POS = -1
COOL_LINE = '★゜・。。・゜゜・。。・゜☆゜・'
COOL_LINE1 = '*★*――――*★*'

def main(target):

    # SUBSET DDR WRT TARGET DRUG
    ddr = get_DD_relations()
    drug_df = ddr.loc[ddr['drug_a'] == target]
    
    if drug_df.empty:
        print("Drug target {} not found. Exiting ...".format(target))
		
        sys.exit()

    # RANK DRUGS WRT scaled_risk_scores (taking natural log of scaled_risk_scores and sorting them to be in descending order)
    pd.options.mode.chained_assignment = None
    drug_df['scaled_risk_score_log'] = np.log(drug_df.loc[:,'scaled_risk_score'].to_numpy(dtype=float))
    drug_df = drug_df.sort_values(by='scaled_risk_score_log', ascending=False)

    print('\n' + COOL_LINE1*3 + ' 𝕽𝖆𝖓𝖐𝖎𝖓𝖌𝖘 𝖋𝖔𝖗 {} '.format(target) + COOL_LINE1*3 + '\n')
    # GET TP AND FP DRUGS BY VALIDATING USING REPODB AND BROAD
    ranks = 1
    for i, candidate_drug_info in drug_df.iterrows():
        candidate_drug = candidate_drug_info[DRUG_B_POS]
        rel_risk_log = candidate_drug_info[LOG_REL_RISK_POS]

        # validation
        approved_diseases = validate_bypair(target, candidate_drug)
        if len(approved_diseases) == 0:
            print("Rank {}: {} with relative risk of {:.3f} of {} | Validation: {}".format(ranks, candidate_drug, rel_risk_log, target, 'FP (つ﹏<。)'))
        else:
            diseases = ','.join([str(disease) for disease in approved_diseases])
            print("Rank {}: {} with relative risk of {:.3f} of {} | Validation: {}, can be repurposed for diseases: {}".format(ranks, candidate_drug, rel_risk_log, target, 'TP (◠‿◠✿)', diseases))
		
        ranks += 1	
    print('\n'+COOL_LINE1 * 9 + '\n')	

if __name__ == "__main__":
	inputParser = argparse.ArgumentParser(description="gives the rankings of drugs by relative risk score associated with common genes")
	inputParser.add_argument("-i", type=str, nargs=1, help="target", required=True)
	
	args = inputParser.parse_args()
	target = args.i[0]
	main(target)


