#hotspot background processing
#libraries
import csv
import pytz
import time
%pip install pandas
import pandas as pd
from datetime import datetime

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
print(df)
print(df.dtypes)

#remove sfr
# /// aq$sfr<-as.numeric(as.character(paste0(aq$sfr)))

#add up particles from bins
#Number concentration of particles between 0.38 and 1 micrometer
n1_out = data['n1_out'] = data.loc[:,['bin0','bin1', 'bin2']].sum(axis=1)
print("N1")
print(n1_out)

#Number concentration of particles between 1 and 3 micrometer
n3_out = data['n3_out'] = data.loc[:,['bin3','bin4', 'bin5', 'bin6']].sum(axis=1)
print("N3")
print(n3_out)

#Number concentration of particles between 1 and 10 micrometer
n10_out = data['n10_out'] = data.loc[:,['bin7','bin8', 'bin9', 'bin10', 'bin11']].sum(axis=1)
print("N10")
print(n10_out)

#something about NAs
#aq<-aq[!is.na(aq$n1),]

#divide by sfr
'''
n1_out1 = (n1/out)/sfr
n3_out1 = (n3/out)/sfr
n10_out1 = (n10/out)/sfr

in essence this is it - need to define sfr and perform division through column though
'''


# this whole thing
'''
aq$lat<-as.numeric(aq$lat)
aq$long<-as.numeric(as.character(paste(aq$long)))
aq<-subset(aq,aq$lat!=0)
aq<-subset(aq, !is.na(aq$pm25))
aq$bin13<-as.numeric(as.character(paste0(aq$bin13)))
aq$bin12<-as.numeric(as.character(paste0(aq$bin12)))
aq$bin14<-as.numeric(as.character(paste0(aq$bin14)))
aq$bin15<-as.numeric(as.character(paste0(aq$bin15)))
aq$bin16<-as.numeric(as.character(paste0(aq$bin16)))
aq$bin17<-as.numeric(as.character(paste0(aq$bin17)))
aq$bin13<-as.numeric(as.character(paste0(aq$bin13)))
aq$bin19<-as.numeric(as.character(paste0(aq$bin19)))
aq$bin20<-as.numeric(as.character(paste0(aq$bin20)))
aq$bin21<-as.numeric(as.character(paste0(aq$bin21)))
aq$bin22<-as.numeric(as.character(paste0(aq$bin22)))
aq$bin23<-as.numeric(as.character(paste0(aq$bin23)))
'''

from rpy2.robjects.packages import importr
    openair = importr('openair')


#and then this whole thing

'''
#Creating date
aq$date<-substring(aq$time,1,10)
aq$date<-as.POSIXct(aq$date, format="%Y-%m-%d")
aq$date<-format(aq$date, "%Y-%m-%d")
aq$date<-as.Date(aq$date," %Y-%m-%d")

aq<-aq[!is.na(aq$n3),]
aq$x<- aq$join_x
aq$y<-aq$join_y
aq$pmw10b_NYCCAS<-aq$join_pmw10b
aq<-aq[ , !grepl( "join" , names(aq) ) ]
'''
