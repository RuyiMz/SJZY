# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 22:01:43 2019

@author: 52665
"""

import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

#sig_ids = ["LPROT001_MCF7_6H:A13","LPROT001_A375_6H:P09","LPROT002_A375_6H:B09"];

#1.Looking for the landmark genes
gene_info = pd.read_csv("LINCS\LINCS Phase II Data--GEO GSE70138\GSE92742_Broad_LINCS_gene_info.txt",sep="\t",dtype=str)
landmark_gene_row_ids = gene_info["pr_gene_id"][gene_info["pr_is_lm"] == "1"]

#2.Reading the dosetemp file
sig_ids = []
dose_temp = pd.read_csv("LINCS_PhaseII_CTRP\LINCS_PhaseII_CTRP_dosetemp_24h.txt",sep="\t",dtype=str);
#print(dose_temp['signature'])
sig_ids = set(dose_temp['signature'])
sig_ids = list(sig_ids)
#print(sig_ids)
#for sig_index in range(len(dose_temp)):
#    sig_ids.append(dose_temp.at[sig_index,'signature'].strip());
#sig_ids = sig_ids.sort()

landmark_only_gctoo = parse("LINCS\LINCS Phase II Data--GEO GSE70138\GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx",rid=landmark_gene_row_ids,cid=sig_ids)

#print(landmark_only_gctoo.data_df.shape)
#print(landmark_only_gctoo.data_df.T)
df = pd.DataFrame(landmark_only_gctoo.data_df.T)
df.to_csv("LINCS_PhaseII_CTRP/LINCS_PhaseII_gene_expression_24h.csv");


