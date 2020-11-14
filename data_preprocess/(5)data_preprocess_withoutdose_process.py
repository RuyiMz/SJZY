# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:17:06 2019

@author: 52665
"""

import pandas as pd
pd.set_option('display.max_columns',15)

#创建文件进行数据记录写入
newRecord_file = open(r"LINCS_CTRP_withoutdose\LINCS_CTRP_withoutdose_24h.txt","w")
newRecord_file.write("record\t"+"signature\t"+"cell_line\t"+"drug\t"+"dose\t"+"time\t"+"cell_viability\n")

dosetemp = pd.read_table(r"LINCS_CTRP_withoutdose\LINCS_CTRP_dosetemp_24h.txt",sep='	')
#print(len(dosetemp['signature']))
dosetemp_unique = dosetemp['signature'].unique()
record_id = 0
for i in range(len(dosetemp_unique)):
    one_signature = dosetemp_unique[i]
    record = dosetemp[dosetemp['signature'] == one_signature]
    record_viability_mean = record['cell_viability'].mean()
    record_dose_mean = record['dose'].mean()
    #print(record['cell_viability'].mean())
    #print(record['dose'].mean())
    newRecord_file.write(str(record_id)+"\t")
    #print(record['signature'].tolist()[0])
    newRecord_file.write(record['signature'].tolist()[0]+"\t")
    newRecord_file.write(record['cell_line'].tolist()[0]+"\t")
    newRecord_file.write(record['drug'].tolist()[0]+"\t")
    newRecord_file.write(str(record_dose_mean)+"\t")
    newRecord_file.write('24'+"\t")
    newRecord_file.write(str(record_viability_mean)+"\n")
    
    record_id = record_id + 1

newRecord_file.close()