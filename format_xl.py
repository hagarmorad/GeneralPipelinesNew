#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 05:44:31 2022

@author: hagar
"""
import pandas as pd 
from xlsxwriter.utility import xl_col_to_name

def save_format_xl(df,num_samples):

        
    df["group_change"] = 0
    for index, row in df.iterrows():
        if row["R/S"] == "S":
            continue
        aa = row["aa_group"].split(",")
        if not aa[0].split("(")[0].strip() == aa[0].split("(")[1].strip():
            df.at[index,'group_change']=1
    
            
    writer = pd.ExcelWriter('reports/mutations.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', index=False)
    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']
    (max_row, max_col) = df.shape
    
    
    red_format = workbook.add_format({'font_color': 'red'})
    grey_format = workbook.add_format({'bg_color':   '#C1C1C1'})
    yellow_format = workbook.add_format({'bg_color':   '#FFFB00'})
    
    worksheet.conditional_format(0,max_col-2, max_row, max_col-2, {'type':     'formula',
                                       'criteria': "=$" + xl_col_to_name(max_col-1) +"2=1",
                                       'format':   red_format})
    
    worksheet.conditional_format(0,1, max_row, max_col-1, {'type':     'cell',
                                    'criteria': 'equal to',
                                    'value':    '"N"',
                                    'format':   grey_format})
                                                                   
    worksheet.conditional_format(0,1, max_row, max_col-1, {'type':     'cell',
                                    'criteria': 'equal to',
                                    'value':    '"-"',
                                    'format':   grey_format})
    worksheet.conditional_format(1,4, max_row, 4+num_samples-1, {'type':     'formula',
                                       'criteria': "=NOT($D2=E2)",
                                       'format':   yellow_format})
    
    worksheet.conditional_format(0,1, max_row, max_col-1, {'type':     'cell',
                                    'criteria': 'equal to',
                                    'value':    '"X"',
                                    'format':   grey_format})
    
    xl_aa_start = 4+num_samples-1+2
    xl_aa_end = 4+num_samples-1+2+num_samples
    
    worksheet.conditional_format(1,xl_aa_start+1, max_row, xl_aa_end, {'type':     'formula',
                                       'criteria': "=NOT($"+xl_col_to_name(xl_aa_start)+"2="+xl_col_to_name(xl_aa_start+1)+"2)",
                                       'format':   yellow_format})


    writer.save()


