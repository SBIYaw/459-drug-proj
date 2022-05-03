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

    tp_value = 0
    fp_value = 0

    # save state (not necessary for the function)
    with open('DR_save_state_random.txt', 'a') as f:
        #f.write('readme')
        for i, candidate_drug_info in drug_df.iterrows():
            candidate_drug = candidate_drug_info[DRUG_B_POS]
            rel_risk_log = candidate_drug_info[LOG_REL_RISK_POS]

            # validation
            approved_diseases = validate_bypair(target, candidate_drug)
            if len(approved_diseases) == 0:
                print("Rank {}: {} with relative risk of {:.3f} of {} | Validation: {}".format(ranks, candidate_drug, rel_risk_log, target, 'FP (つ﹏<。)'))
                fp_value += 1/ranks
                f.write(target + ' ' + candidate_drug  + ' [' +  np.asarray(drug_df[drug_df['drug_b'] == candidate_drug]['scaled_risk_score'])[0]  + '] ' +  str(ranks)  + ' ' +  'False \n')
            else:
                diseases = ','.join([str(disease) for disease in approved_diseases])
                print("Rank {}: {} with relative risk of {:.3f} of {} | Validation: {}, can be repurposed for diseases: {}".format(ranks, candidate_drug, rel_risk_log, target, 'TP (◠‿◠✿)', diseases))
                tp_value += 1/ranks
                f.write(target + ' ' + candidate_drug  + ' [' +  np.asarray(drug_df[drug_df['drug_b'] == candidate_drug]['scaled_risk_score'])[0]  + '] ' +  str(ranks)  + ' ' +  'True \n')
            
            ranks += 1
    
        f.write('-------------------------------------------------------------------------------\n')
    f.close()
    
    #print('\n TP:', round(tp_value,3))
    #print('\n FP:', round(fp_value,3))
    precision = tp_value / (tp_value + fp_value)
    print('\n𝕻𝖗𝖊𝖈𝖎𝖘𝖎𝖔𝖓:･ﾟ✧ ', round(precision,3))

    print('\n'+COOL_LINE1 * 16 + '\n')	
    return precision

if __name__ == "__main__":
	inputParser = argparse.ArgumentParser(description="gives the rankings of drugs by relative risk score associated with common genes")
	inputParser.add_argument("-i", type=str, nargs=1, help="target", required=True)
	
	args = inputParser.parse_args()
	target = args.i[0]
	main(target)


