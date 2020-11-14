# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 10:07:53 2019

@author: 52665
"""

import pandas as pd
import math
pd.set_option('display.max_columns',20)

LINCS_CTRP_dose_file = open("LINCS_PhaseI_CTRP/LINCS_PhaseI_CTRP_dose_24h.txt","w")
LINCS_CTRP_dose_logfile = open("LINCS_PhaseI_CTRP/LINCS_PhaseI_CTRP_dose_log_24h.txt","w")

LINCS_CTRP_dose_file.write("new_matched_record_id\t"+"old_matched_record_id\t"+"sig_id\t"+"pert_id\t"+"pert_iname\t"+"pert_type\t"+"cell_id\t"+"pert_idose\t"+"pert_itime\t"+"experiment_id\t"+"conc_pts_fit\t"+"apparent_ec50_umol\t"+"pred_pv_high_conc\t"+"area_under_curve\t"+"master_cpd_id\t"+"cpd_name\t"+"broad_cpd_id\t"+"top_test_conc_umol\t"+"master_ccl_id\t"+"ccl_name\t"+"cpd_pred_pv\t"+"cpd_avg_pv\t"+"cpd_conc_umol"+"\n")
LINCS_CTRP_dose_logfile.write("LINCS_CTRP_oldid"+"\n")

#1. Reading the LINCS_PhaseI data
LINCS_PhaseII_CTRP =pd.read_table(r"LINCS_PhaseI_CTRP\LINCS_PhaseI_CTRP_24h.txt",sep="	");
#print(len(LINCS_PhaseII_CTRP))
#aa = set(LINCS_PhaseII_CTRP[:]["sig_id"])
#print(len(aa))
#print(LINCS_PhaseII_CTRP[657:662])
#print(LINCS_PhaseII_CTRP.dtypes)

#2. Reading the CTRP data
CTRP_dose = pd.read_table(r"CTRP\CTRPv2.0_2015_ctd2_ExpandedDataset\v20.data.per_cpd_post_qc.txt",sep="	");
#print(CTRP_dose[0:5])
#print(CTRP_dose.dtypes)

#3. Filtering the corresponding concentration data
LINCS_CTRP_newRecord_id = 0
for LINCS_CTRP_index in range(len(LINCS_PhaseII_CTRP)):
    LINCS_dose_value = float(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_idose"])
    LINCS_CTRP_concNum = LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"conc_pts_fit"]
    LINCS_CTRP_experiment = LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"experiment_id"]
    LINCS_CTRP_cpd = LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"master_cpd_id"]
    
    #print(str(LINCS_CTRP_index)+"\n")
    CTRP_dose_filter1 = CTRP_dose[CTRP_dose["experiment_id"] == LINCS_CTRP_experiment] 
    CTRP_dose_filter2 = CTRP_dose_filter1[CTRP_dose["master_cpd_id"] == LINCS_CTRP_cpd]
    
    #Renumbering the filtered CTRP_dose_filter data
    CTRP_dose_filter3 = CTRP_dose_filter2.reset_index(drop=True)
    #print(CTRP_dose_filter2)
    #print(len(CTRP_dose_filter3))
    for CTRP_dose_index in range(len(CTRP_dose_filter3)):
        CTRP_dose_value = float(CTRP_dose_filter3.at[CTRP_dose_index,"cpd_conc_umol"])
        #Calculating the concentration difference
        dose_difference = abs(math.log10(CTRP_dose_value) - math.log10(LINCS_dose_value))
        if(dose_difference <= 0.2):
            LINCS_CTRP_dose_file.write(str(LINCS_CTRP_newRecord_id)+"\t")
            
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"matched_record_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"sig_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_iname"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_type"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"cell_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_idose"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pert_itime"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"experiment_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"conc_pts_fit"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"apparent_ec50_umol"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"pred_pv_high_conc"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"area_under_curve"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"master_cpd_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"cpd_name"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"broad_cpd_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"top_test_conc_umol"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"master_ccl_id"])+"\t")
            LINCS_CTRP_dose_file.write(str(LINCS_PhaseII_CTRP.at[LINCS_CTRP_index,"ccl_name"])+"\t")
            
            LINCS_CTRP_dose_file.write(str(CTRP_dose_filter3.at[CTRP_dose_index,"cpd_pred_pv"])+"\t")
            LINCS_CTRP_dose_file.write(str(CTRP_dose_filter3.at[CTRP_dose_index,"cpd_avg_pv"])+"\t")
            LINCS_CTRP_dose_file.write(str(CTRP_dose_filter3.at[CTRP_dose_index,"cpd_conc_umol"])+"\n")
            
            LINCS_CTRP_newRecord_id = LINCS_CTRP_newRecord_id + 1
    #Saving the data running status to log file
    LINCS_CTRP_dose_logfile.write(str(LINCS_CTRP_index)+"\n")
    print(str(LINCS_CTRP_index)+"\n")
    #if(LINCS_CTRP_index == 2):
    #    break

LINCS_CTRP_dose_file.close()
LINCS_CTRP_dose_logfile.close()