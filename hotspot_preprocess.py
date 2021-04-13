%%capture
#libraries
import csv
import pytz
import time
%pip install pandas
import pandas as pd
from datetime import datetime
import os
os.environ['R_HOME'] = '/path/to/R'

#analysis of dataset collected by two low-cost OPC-N2s on two trash-trucks in NYC

#set timezone
tz = pytz.timezone('US/Eastern')

#read in OPC data
aq = "AQ_orgfid.csv"
data = pd.read_csv(aq, engine='python')
print("aq read")

#format time
times = data.loc[:,['time']]
df = pd.DataFrame(times)
df['time'] = pd.to_datetime(df['time'], unit='s')
#print(df)
#print(df.dtypes)


#sfr
sfr_in = data.loc[:,['sfr']]
sfr_df = pd.DataFrame(sfr_in)
#r code (below) says remove - investigate paste0 
# /// aq$sfr<-as.numeric(as.character(paste0(aq$sfr)))

#add up particles from bins
#Number concentration of particles between 0.38 and 1 micrometer
n1_out = data['n1_out'] = data.loc[:,['bin0','bin1', 'bin2']].sum(axis=1)
#print("N1")
#print(n1_out)

#Number concentration of particles between 1 and 3 micrometer
n3_out = data['n3_out'] = data.loc[:,['bin3','bin4', 'bin5', 'bin6']].sum(axis=1)
#print("N3")
#print(n3_out)

#Number concentration of particles between 1 and 10 micrometer
n10_out = data['n10_out'] = data.loc[:,['bin7','bin8', 'bin9', 'bin10', 'bin11']].sum(axis=1)
#print("N10")
#print(n10_out)

#something about NAs
#aq<-aq[!is.na(aq$n1),]

#divide by sfr 
'''
data['n1_out1'] = data[n1_out]/data[sfr_df]
print(n10_out)
n3_out1 = (n3_out)/sfr_df
print(n10_out)
n10_out1 = (n10_out)/sfr_df
print(n10_out)
'''
