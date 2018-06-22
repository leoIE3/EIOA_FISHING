#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 11:23:33 2018

@author: leoie
"""
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from numpy import inf 
idx=pd.IndexSlice
print(os.getcwd())
wd=os.getcwd()
#%%
#raw_data_dir=os.getcwd()
#iot_filename="mrIot_3.3_2011.txt" #this is how the transaction matrix Z is called in exiobase
#fd_filename="mrFinalDemand_3.3_2011.txt" #This is how the Y final demand is called in exiobase
#em_filename="mrEmissions_3.3_2011.txt"
#fd_em_filename="mrFDEmissions_3.3_2011.txt"
#fi=pd.read_csv(fi_filename,sep='\t',index_col=[0,1],header=[0,1])
def pickling():
    z=pd.read_pickle('z.pkl')
    y=pd.read_pickle('y.pkl')
    fi=pd.read_pickle('fi.pkl')
    em=pd.read_pickle('em.pkl')
    em_y=pd.read_pickle('em_y.pkl')
    return (z,y,fi,em,em_y)

z,y,fi,em,em_y=pickling()
va_index=list(range(9))
co2_index=[0,28,77,78]
employment=[9,10,11,12,13,14]

products=[]
for i in range(len(z)):
    products.append(z.index.values[i][1])
products=list(set(products))


va=fi.iloc[va_index]
va_jobs=fi.iloc[employment].sum(0)
em_co2=em.iloc[co2_index]
em_co2_y=em_y.iloc[co2_index]
em_co2_air=em_co2.iloc[0]
y_categories=[]
for i in range(len(y.columns)):
    y_categories.append(y.columns[i][1])
y_categories=list(set(y_categories))

 
#sum in rows .sum(0) sum in columns .sum(1)
x_out=z.sum(1)+y.sum(1) #total gross output per FD
x_inp=z.sum(0)+va.sum(0) #total gross output per VA
#x_inp_=z.sum(0)+fi.sum(0)
#%%

"Calculate Inverse and A matrix"
inv_x=1/x_out
inv_x[inv_x==inf]=0
diag=np.diag(inv_x)
A=z.dot(diag)
A.columns=z.columns
#A_NO=A.loc[:,'NO']

"Calculate B vector and air emissions vector"
B=pd.DataFrame(np.dot(em_co2,diag))
B.columns=A.columns
B.index=em_co2.index
B_co2_air=pd.DataFrame(B.loc[idx[B.index.values[0]],])

"Calculate I matrix vector and Leontief inverse"
I=np.identity(len(z))
L=np.subtract(I,A)
L=np.linalg.inv(np.array(L))
L=pd.DataFrame(L)
L.columns=A.columns
L.index=A.index
#%%
"Calculate M matrix and can be either done with b' or diag(b)"
#M=pd.DataFrame(np.dot(B,L)) #because of hotspot analysis
M=pd.DataFrame(np.dot(np.transpose(B_co2_air),L)) #because of contribution analysys
M=M.transpose()
M.columns=B_co2_air.columns
M.index=B_co2_air.index

#%%
"Calculate Delta X matrix"
Delta_x=pd.DataFrame(np.dot(L,y))
Delta_x.columns=y.columns
Delta_x.index=L.index

"Calculate Delta R matrix can be either done with b' or diag(b)"
#Delta_r=pd.DataFrame(np.dot(np.dot(B,L),y))
#Delta_r.columns=y.columns
#Delta_r.index=B.index

Delta_r=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),y))
Delta_r.columns=y.columns
Delta_r.index=np.transpose(B_co2_air).index
#%%
"Calculate total impact"
HH=int(y_categories.index('Final consumption expenditure by households'))
y_sum=pd.DataFrame()
for i in range(len(y_categories)):
    y_sum[y_categories[i]]=y.loc[:,idx[:,y_categories[i],:]].sum(1)
y_sum_HH=y_sum.iloc[:,idx[HH]]
y_sum_HH.index=y.index
Delta_r=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),y_sum_HH))

y_sum=y_sum.sum(1)
contribution=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),np.diag(y_sum)))
contribution.columns=z.columns
contribution=np.transpose(contribution)
#jobs=pd.DataFrame(np.dot(np.dot(va_jobs,L),y_sum_HH))
#jobs.index=va_jobs.index

jobs_diag=pd.DataFrame(np.dot(np.dot(np.diag(va_jobs),L),y_sum_HH))
jobs_diag.index=y_sum_HH.index


B_co2_air1=B.loc[idx[B.index.values[0]],]
hotspot=pd.DataFrame(np.dot(np.dot(np.diag(B_co2_air1),L),y_sum_HH))
#contribution.columns=z.columns


#%%
"Specific CO2 emissions for fish sectors"
fish_index=['Fish and other fishing products; services incidental of fishing','Fish products']
CO2_fish_sector=[]
codes=['PT','NO','ES']
for i in range(len(codes)):
    for j in range(len(fish_index)):
        CO2_fish_sector.append(np.transpose(B_co2_air).loc[:,idx[codes[i],fish_index[j],]])


#%%
#FISH
#fish_co2.transpose()
#fish_co21=fish_co2.transpose()


#hotspot_fish=hotspot.loc[idx[:,['Fish and other fishing products; services incidental of fishing','Fish products'],:],:]    
#hotspot_fish_sum=pd.DataFrame()
#for i in range(len(y_categories)):
#    hotspot_fish_sum[y_categories[i]]=hotspot_fish.loc[:,idx[:,y_categories[i],:]].sum(1)
    
"increase demand for Fish and other fishing products; services incidental of fishing and Fish products"    

#increase by 50% 
delta_fish=[]
for i in range(len(codes)):
    delta_fish.append(list(y.loc[idx[codes[i],[fish_index[0],fish_index[1]],:],\
                  idx['PT',y_categories[HH],:]].values))


#assing change in FD to other categories    
PT1=('PT',fish_index[0],'M.EUR') #fish and others
PT2=('PT',fish_index[1],'M.EUR') #fish products
NO1=('NO',fish_index[0],'M.EUR')
NO2=('NO',fish_index[1],'M.EUR')
ES1=('ES',fish_index[0],'M.EUR')
ES2=('ES',fish_index[1],'M.EUR')
list_fish=[PT1,PT2,NO1,NO2,ES1,ES2]    

"Halve values for PT y fish related" 
y.at[PT1,('PT',y_categories[HH])]=float(delta_fish[0][0])/2   #fish and other
y.at[PT2,('PT',y_categories[HH])]=float(delta_fish[0][1])/2   #fish products

#case 1 - to NO     
y.at[NO1,('PT',y_categories[HH])]=float(delta_fish[0][0])/2+y.at[NO1,('PT',y_categories[HH])]   
y.at[NO2,('PT',y_categories[HH])]=float(delta_fish[0][1])/2+y.at[NO2,('PT',y_categories[HH])]

y_sum_NO=pd.DataFrame()
for i in range(len(y_categories)):
    y_sum_NO[y_categories[i]]=y.loc[:,idx[:,y_categories[i],:]].sum(1)
#y_sum_HH_NO=y_sum_NO.iloc[:,idx[y_categories[HH]]]
y_sum_HH_NO=y_sum_NO.iloc[:,idx[HH]]
Delta_r_NO=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),y_sum_HH_NO))
y_sum_NO=y_sum_NO.sum(1)
contribution_NO=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),np.diag(y_sum_NO)))
contribution_NO.columns=z.columns
contribution_NO=np.transpose(contribution_NO)
hotspot_NO=pd.DataFrame(np.dot(np.dot(np.diag(B_co2_air1),L),y_sum_HH_NO))
jobs_diag_NO=pd.DataFrame(np.dot(np.dot(np.diag(va_jobs),L),y_sum_HH_NO))
jobs_diag_NO.index=y_sum_HH.index
"reset values for NO y fish related"

#case 1 - to Es  
#reinitialize NO
y.at[NO1,('PT',y_categories[HH])]=float(delta_fish[1][0])  
y.at[NO2,('PT',y_categories[HH])]=float(delta_fish[1][1])

y.at[ES1,('PT',y_categories[HH])]=float(delta_fish[0][0])/2+y.at[ES1,('PT',y_categories[HH])]   
y.at[ES2,('PT',y_categories[HH])]=float(delta_fish[0][1])/2+y.at[ES2,('PT',y_categories[HH])]

y_sum_ES=pd.DataFrame()
for i in range(len(y_categories)):
    y_sum_ES[y_categories[i]]=y.loc[:,idx[:,y_categories[i],:]].sum(1)
#y_sum_HH_ES=y_sum_ES.iloc[:,idx[y_categories[HH]]]
y_sum_HH_ES=y_sum_ES.iloc[:,idx[HH]]
Delta_r_ES=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),y_sum_HH_ES))
y_sum_ES=y_sum_ES.sum(1)
contribution_ES=pd.DataFrame(np.dot(np.dot(np.transpose(B_co2_air),L),np.diag(y_sum_ES)))
contribution_ES.columns=z.columns
contribution_ES=np.transpose(contribution_ES)
hotspot_ES=pd.DataFrame(np.dot(np.dot(np.diag(B_co2_air1),L),y_sum_HH_ES))
jobs_diag_ES=pd.DataFrame(np.dot(np.dot(np.diag(va_jobs),L),y_sum_HH_ES))
jobs_diag_ES.index=y_sum_HH.index

contribution.to_csv('contribution.csv',sep='\t')
hotspot.to_csv('hotspot.csv',sep='\t')

contribution_NO.to_csv('contributionNO.csv',sep='\t')
hotspot_NO.to_csv('hotspotNO.csv',sep='\t')
contribution_ES.to_csv('contributionES.csv',sep='\t')
hotspot_ES.to_csv('hotspotES.csv',sep='\t')
print(jobs_diag.sum(0),jobs_diag_NO.sum(0),jobs_diag_ES.sum(0))

jobs_plot=[]
for i in range(len(codes)):
    jobs_plot.append(list(jobs_diag.loc[idx[codes[i],[fish_index[0],fish_index[1]],:]].values))
    jobs_plot.append(list(jobs_diag_NO.loc[idx[codes[i],[fish_index[0],fish_index[1]],:]].values))                
    jobs_plot.append(list(jobs_diag_ES.loc[idx[codes[i],[fish_index[0],fish_index[1]],:]].values))





width = 0.35  # the width of the bars

fig, ax = plt.subplots()
#rects1 = ax.bar([int(jobs_plot[0][0]),int(jobs_plot[3][0]),int(jobs_plot[6][0])], width,
#                color='SkyBlue', label='PT C_FISH')
bar1 = ax.bar('PT',[2,2,2] ,width,
                color='SkyBlue', label='PT C_FISH')
bar2 = ax.bar([4, 4, 4],, width,
                color='Red', label='PT C_FISH')
#rects2 = ax.bar(int(jobs_plot[1][0]), int(jobs_plot[4][0]),int(jobs_plot[7][0]), width,
                #color='IndianRed', label='NO C_FISH')
#rects3 = ax.bar(int(jobs_plot[2][0]), int(jobs_plot[4][0]),int(jobs_plot[8][0]), width, 
                #color='Black', label='ES C_FISH')
# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Thousand jobs')
ax.set_title('Job creation')
ax.set_xticks(ind)
ax.set_xticklabels(('PT', 'NO', 'ES'))
ax.legend()
plt.show()
