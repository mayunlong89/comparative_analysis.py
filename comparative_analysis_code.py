#!/usr/bin/env python
#Author: Yunlong Ma
#Usage: Code for comparative analysis which is a part of an analysis of integrative genomics analysis of CAD


import numpy


#Part I Load data

#Read data 
f1 = open("Zeller_sig_CAD.txt","r+")
f2 = open("Dixon_sig_CAD.txt","r+")
f3 = open("MAGMA_sig_CAD.txt","r+")
f4 = open("MAGMA_sig_NULL.txt","r+")


#Part II Extract genes according to different P values

#Set up 3 thresholds of P values as 3 comparative points
#P values: 0.05, 0.01, and 0.001 from  Sherlock analysis

Geneset_1 =[]
Geneset_2 =[]
Geneset_3 =[]
 
#f1 Sherlock analysis of integrating GWAS with Zeller et al. eQTL data.
for line in f1:
    if line[0:3] !="Gene":
        dd = line.strip().split()
        Geneset_1.append(dd[0])      
        if dd[1]!="P" and float(dd[1])<0.01:
            Geneset_2.append(dd[0])          
        if dd[1]!="P" and float(dd[1])<0.001:
            Geneset_3.append(dd[0])

Geneset2_1 =[]
Geneset2_2 =[]
Geneset2_3 =[]
        
#f2 Sherlock analysis of integrating GWAS with Dixon et al. eQTL data.
for line2 in f2:
    if line2[0:3] !="Gene":
        dd2 = line2.strip().split()
        Geneset2_1.append(dd[0])      
        if dd2[1]!="P" and float(dd2[1])<0.01:
            Geneset2_2.append(dd2[0])          
        if dd2[1]!="P" and float(dd2[1])<0.001:
            Geneset2_3.append(dd2[0])


#f3 MAGMA analysis of GWAS on CAD (MAGMA-based P value <0.05)
Geneset3 =[]
for i in f3:
    if i[0:3]!="Gene":
        dd3 = i.strip()
        Geneset3.append(dd3)
        

#f4 MAGMA analysis of GWAS on null trait (MAGMA-based P value <0.05)
Geneset4=[]
for j in f4:
    if j[0:3]!="Gene":
        dd4 = j.strip()
        Geneset4.append(dd4)



#Part III Calculate overlapped gene rates at 3 thresholds of P values 


#Establish a novel function for calculate the overlapped gene rates
def fun_rate(x,y,z):
    sum1 = len(x)
    Combined = []
    for num in x:
        if num in y:
            Combined.append(num)
    sum2 = len(Combined)
    rate = sum2/sum1
    
    Combined_null = []
    for num2 in x:
        if num2 in z:
            Combined_null.append(num2)
    sum_null_1 = len(Combined_null)
    if sum_null_1 == 0:
        sum_null_1 = 0.0001
    rate2 = sum_null_1/sum1
    return(rate,rate2)
 
           
#For Sherlock analysis of Zeller et al eQTL vs. MAGMA  
#At the threshold of P value = 0.05
Rate1 = fun_rate(Geneset_1,Geneset3,Geneset4)
#At the threshold of P value = 0.01
Rate2 = fun_rate(Geneset_2,Geneset3,Geneset4)
#At the threshold of P value = 0.001
Rate3 = fun_rate(Geneset_3,Geneset3,Geneset4)
CAD_group = [Rate1[0],Rate2[0],Rate3[0]]
Null_group =  [Rate1[1],Rate2[1],Rate3[1]]

#For Sherlock analysis of Dixon et al eQTL vs. MAGMA
#At the threshold of P value = 0.05
Rate4 = fun_rate(Geneset2_1,Geneset3,Geneset4)
#At the threshold of P value = 0.01
Rate5 = fun_rate(Geneset2_2,Geneset3,Geneset4)
#At the threshold of P value = 0.001
Rate6 = fun_rate(Geneset2_3,Geneset3,Geneset4)
CAD_group_2 = [Rate4[0],Rate5[0],Rate6[0]]
Null_group_2 =  [Rate4[1],Rate5[1],Rate6[1]]



#paired t-test for Sherlock and MAGMA on CAD VS. Sherlock and MAGMA on null trait
from scipy.stats import ttest_rel
ttest_rel(CAD_group,Null_group)



if __name__ == "__main__":
    print(ttest_rel(CAD_group,Null_group))

#End



