# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 11:58:50 2019

@author: 52665

1 Matching the LINCS-L1000-PhaseII and CTRP data
"""

import pandas as pd
import math
pd.set_option('display.max_columns',15)

def CTRP_match():
    #1. Reading the file
    curves_post_qc = pd.read_table(r"CTRP\CTRPv2.0_2015_ctd2_ExpandedDataset\v20.data.curves_post_qc.txt",sep='	')
    curves_post_qc_filter = curves_post_qc[['experiment_id','conc_pts_fit','apparent_ec50_umol','pred_pv_high_conc','area_under_curve','master_cpd_id']]
    #print(curves_post_qc_filter)
    
    per_experiment = pd.read_table(r"CTRP\CTRPv2.0_2015_ctd2_ExpandedDataset\v20.meta.per_experiment.txt",sep='	')
    per_experiment_filter = per_experiment[['experiment_id','master_ccl_id']]
    per_experiment_filter = per_experiment_filter.drop_duplicates()
    #print(per_experiment_filter[360:380])
    
    per_compound = pd.read_table(r"CTRP\CTRPv2.0_2015_ctd2_ExpandedDataset\v20.meta.per_compound.txt",sep='	')
    per_compound_filter = per_compound[['master_cpd_id','cpd_name','broad_cpd_id','top_test_conc_umol']]
    #print(per_compound_filter)
    
    per_cell_line = pd.read_table(r"CTRP\CTRPv2.0_2015_ctd2_ExpandedDataset\v20.meta.per_cell_line.txt",sep='	')
    per_cell_line_filter = per_cell_line[['master_ccl_id','ccl_name']]
    #print(per_cell_line_filter)
    
    #2.Merging the data
    #2.1 Matching the drug information according to master_cpd_id
    viability_cpd = pd.merge(curves_post_qc_filter,per_compound_filter,on=['master_cpd_id','master_cpd_id'])
    #print(viability_cpd[0:5])
    
    #2.2 Matching the cell line based on the experiment_id
    viability_cpd_cclid = pd.merge(viability_cpd,per_experiment_filter,on=['experiment_id','experiment_id'])
    #print(viability_cpd_cclid[0:5])
    
    #2.3 Pairing the cell line information according to master_ccl_id
    viability_cpd_ccl = pd.merge(viability_cpd_cclid,per_cell_line_filter,on=['master_ccl_id','master_ccl_id'])
    print(len(viability_cpd_ccl))

    return viability_cpd_ccl

#0.Creating the data files and log files
LINCS_PhaseII_CTRP_file = open("LINCS_PhaseII_CTRP/LINCS_PhaseII_CTRP_24h.txt","w")
LINCS_PhaseII_CTRP_logfile = open("LINCS_PhaseII_CTRP/LINCS_PhaseII_CTRP_log_24h.txt","w")
LINCS_PhaseII_CTRP_file.write("matched_record_id\t"+"sig_id\t"+"pert_id\t"+"pert_iname\t"+"pert_type\t"+"cell_id\t"+"pert_idose\t"+"pert_itime\t"+"experiment_id\t"+"conc_pts_fit\t"+"apparent_ec50_umol\t"+"pred_pv_high_conc\t"+"area_under_curve\t"+"master_cpd_id\t"+"cpd_name\t"+"broad_cpd_id\t"+"top_test_conc_umol\t"+"master_ccl_id\t"+"ccl_name"+"\n")
LINCS_PhaseII_CTRP_logfile.write("CTRP_curves_post_id"+"\n")

viability_cpd_ccl_info = CTRP_match()
#print(viability_cpd_ccl_info.at[1,'apparent_ec50_umol'])

#1.Reading the LINCS-L1000 file
sig_info = pd.read_csv("LINCS\LINCS Phase II Data--GEO GSE70138\GSE70138_Broad_LINCS_sig_info.txt", sep="\t")
#1.1 Data processing
for sig_info_index in range(len(sig_info)):
    dose_value = sig_info.at[sig_info_index,"pert_idose"].strip(" um")
    sig_info.at[sig_info_index,"pert_idose"] = dose_value
#print(sig_info)
sig_info_exclude_dose = sig_info[sig_info["pert_idose"] != "-666"]
sig_info_exclude_dose = sig_info_exclude_dose[sig_info["pert_idose"] != "0.0"]
sig_info_exclude_dose = sig_info_exclude_dose[sig_info["pert_itime"] == "24 h"]
sig_info_exclude_dose = sig_info_exclude_dose.reset_index(drop=True)
#print(sig_info_exclude_dose)

#2.Pairing with the CTRP data to find the sig_id
record_id = 0
all_matched_sig_id = []
for CTRP_index in range(len(viability_cpd_ccl_info)):
    CTRP_cell_line = viability_cpd_ccl_info.at[CTRP_index,'ccl_name'].strip().encode("utf-8")
    CTRP_drug = viability_cpd_ccl_info.at[CTRP_index,'broad_cpd_id'].strip().encode("utf-8")
    #CTRP_dose = float(viability_cpd_ccl_info.at[CTRP_index,'top_test_conc_umol'])
    for LINCS_index in range(len(sig_info_exclude_dose)):
        LINCS_cell_line = sig_info_exclude_dose.at[LINCS_index,'cell_id'].strip().encode("utf-8")
        LINCS_drug = sig_info_exclude_dose.at[LINCS_index,'pert_id'].strip().encode("utf-8")
        #LINCS_dose = float(sig_info_exclude_dose.at[LINCS_index,'pert_idose'])
        #Calculating the concentration difference
        #dose_difference = abs(math.log10(CTRP_dose) - math.log10(LINCS_dose))
        if ((CTRP_cell_line == LINCS_cell_line) and (CTRP_drug == LINCS_drug)):
            LINCS_PhaseII_CTRP_file.write(str(record_id)+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'sig_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'pert_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'pert_iname'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'pert_type'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'cell_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'pert_idose'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(sig_info_exclude_dose.at[LINCS_index,'pert_itime'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'experiment_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'conc_pts_fit'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'apparent_ec50_umol'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'pred_pv_high_conc'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'area_under_curve'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'master_cpd_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'cpd_name'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'broad_cpd_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'top_test_conc_umol'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'master_ccl_id'])+"\t")
            LINCS_PhaseII_CTRP_file.write(str(viability_cpd_ccl_info.at[CTRP_index,'ccl_name'])+"\n")
            record_id = record_id + 1
            all_matched_sig_id.append(sig_info_exclude_dose.at[LINCS_index,'sig_id'])
    LINCS_PhaseII_CTRP_logfile.write(str(CTRP_index)+"\n")
    print(CTRP_index,":",len(all_matched_sig_id))

print(all_matched_sig_id)

#3.Closing the file stream
LINCS_PhaseII_CTRP_file.close()
LINCS_PhaseII_CTRP_logfile.close()








