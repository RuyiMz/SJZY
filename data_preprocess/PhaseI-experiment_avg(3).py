# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:44:40 2019

@author: 52665
"""

import pandas as pd
pd.set_option('display.max_columns',20)

#0. Getting the data and writing it to file
record = 0
LINCS_PhaseII_CTRP_file = open("LINCS_PhaseI_CTRP/LINCS_PhaseI_CTRP_dosetemp_24h.txt","w")
LINCS_PhaseII_CTRP_file.write("record\t"+"signature\t"+"cell_line\t"+"drug\t"+"dose\t"+"time\t"+"cell_viability\n")


LINCS_PhaseII_CTRP =pd.read_table(r"LINCS_PhaseI_CTRP\LINCS_PhaseI_CTRP_dose_24h.txt",sep="	");
#1. Using the collection method to filter the signature id and remove the duplicate signature
distinct_sigids = set(LINCS_PhaseII_CTRP[:]["sig_id"])
#print(len(distinct_sigid))
for sigid in distinct_sigids:
    #sigid="LPROT003_A549_6H:E18"
    #print(sigid)
    LINCS_CTRP_subset = LINCS_PhaseII_CTRP[LINCS_PhaseII_CTRP["sig_id"] == sigid]
    LINCS_CTRP_subset = LINCS_CTRP_subset.reset_index(drop=True)
    #print(len(LINCS_CTRP_subset))
    #print(LINCS_CTRP_subset)
    
    #2. Number of screening concentrations
    #print(len(set(LINCS_CTRP_subset["cpd_conc_umol"])))
    cpd_conc_singles = set(LINCS_CTRP_subset["cpd_conc_umol"])
    for cpd_conc in cpd_conc_singles:
        #Writing the record information to the file
        LINCS_PhaseII_CTRP_file.write(str(record)+"\t")
        LINCS_PhaseII_CTRP_file.write(LINCS_CTRP_subset["sig_id"][0]+"\t")
        LINCS_PhaseII_CTRP_file.write(LINCS_CTRP_subset["ccl_name"][0]+"\t")
        LINCS_PhaseII_CTRP_file.write(LINCS_CTRP_subset["broad_cpd_id"][0]+"\t")
        LINCS_PhaseII_CTRP_file.write(str(cpd_conc)+"\t")
        LINCS_PhaseII_CTRP_file.write("6\t")
        
        select_conc = LINCS_CTRP_subset[LINCS_CTRP_subset["cpd_conc_umol"] == cpd_conc]
        #select_conc_length = len(select_conc)
        #select_conc = select_conc.reset_index(drop=True)
        #print(select_conc["cpd_pred_pv"].mean())
        LINCS_PhaseII_CTRP_file.write(str(select_conc["cpd_pred_pv"].mean())+"\n")
        #print(select_conc)
    
        record = record + 1
        
        
        #break

LINCS_PhaseII_CTRP_file.close()