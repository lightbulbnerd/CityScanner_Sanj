#hotspot background processing

import csv
import pandas as pd
from datetime import datetime

#analysis of dataset collected by two low-cost OPC-N2s on two trash-trucks in NYC
date_default_timezone_set("America/NewYork");
#Reading in the data

#read in csv
with open("/Users/sanjanapaul/Downloads/AQ_orgfid.csv") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            print('Line read!')
            line_count += 1
    print(f'Processed {line_count} lines.')
#aq<-read.csv("nyc_aq_joined.csv")
print("aq read")

# just ran chunk 3 in the r code, which is what this is. compare actual and
# desired results and adjust as necessary

# cleaning data below -
# subset=None:every column is used to determine if two rows are different;
# to change that specify columns as an array
# inplace=True:the data structure is changed and duplicate rows are gone
'''
aq_in = "AQ_orgfid.csv"
aq = "AQ_orgfid_clean.csv"
df = pd.read_csv(file_in, sep="\t or ,")
df.drop_duplicates(subset=None, inplace=True)
df.to_csv(file_name_output, index=False)

#Formatting string from CSV to time
pd.to_datetime(aq$time, origin="1970-01-01 00:00")
pd.to_datetime(aq$date, format'"%Y-%m-%d")
'''
