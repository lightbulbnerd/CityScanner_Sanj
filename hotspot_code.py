

'''
Python version of City Scanner hotspot code, based on Priyanka's original
R analysis
'''
#what all libraries should we get in here?
import csv
import pandas as pd
from datetime import datetime

#analysis of dataset collected by two low-cost OPC-N2s on two trash-trucks in NYC
date_default_timezone_set("America/NewYork");
#Reading in the data

#comment everything else out and try this chunk first
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
aq_in = "AQ_orgfid.csv"
aq = "AQ_orgfid_clean.csv"
df = pd.read_csv(file_in, sep="\t or ,")
df.drop_duplicates(subset=None, inplace=True)
df.to_csv(file_name_output, index=False)

#Formatting string from CSV to time
pd.to_datetime(aq$time, origin="1970-01-01 00:00")
pd.to_datetime(aq$date, format'"%Y-%m-%d")

'''
#Creating date
aq$date<-substring(aq$time,1,10)
aq$date<-as.POSIXct(aq$date, format="%Y-%m-%d") #
aq$date<-format(aq$date, "%Y-%m-%d")
aq$date<-as.Date(aq$date," %Y-%m-%d")
#Creating time
aq$time<-as.POSIXct(aq$time,"%H:%M:%S" )
#Creating hour
aq$h <-substring(aq$time, 12,13)
aq$h<-as.integer(aq$h)
#Remove SFR
aq$sfr<-as.numeric(as.character(paste0(aq$sfr)))
aq$bin4<-as.numeric(as.character(paste0(aq$bin4)))
aq$bin5<-as.numeric(as.character(paste0(aq$bin5)))
aq$bin6<-as.numeric(as.character(paste0(aq$bin6)))
aq$bin11<-as.numeric(as.character(paste0(aq$bin11)))
aq<-aq[!is.na(aq$sfr),]
aq<-aq[aq$sfr!=0,]
#Number concentration of particles between 0.38 and 1 micrometer
aq$n1<-aq$bin0+aq$bin1+aq$bin2
#Number concentration of particles between 1 and 3 micrometer
aq$n3<-aq$bin3+aq$bin4+aq$bin5+aq$bin6
#Number concentration of particles between 1 and 10 micrometer
aq$n10<-aq$bin3+aq$bin4+aq$bin5+aq$bin6+aq$bin7+aq$bin8+aq$bin9+aq$bin10+aq$bin11
aq<-aq[!is.na(aq$n1),]

#Dividing by SFR
aq$n1<-aq$n1/aq$sfr
aq$n3<-aq$n3/aq$sfr
aq$n10<-aq$n10/aq$sfr

#Making sure the dates are correct
aq<-subset(aq,aq$date>="2020-01-01")


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


require(openair)
aq$date<-aq$time
nyc_aq<-subset(aq, select=c(date, pm1, pm25, pm10, n1, n3, n10))
aq_hourly<-timeAverage(nyc_aq, avg.time="hour")
summary(aq_hourly)

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


#Cleaning AQ data
{r}
nrow(aq)
#aq$long<-as.numeric(as.character(paste(aq$long)))
#aq$lat<-as.numeric(as.character(paste(aq$lat)))
aq1<-aq[aq$long< -72,]
aq1<-aq1[aq1$lat>40,]
aq<-aq1
rm(aq1)
nrow(aq)

#Function for Splitting line into 30 m segments
#https://stackoverflow.com/questions/38700246/how-do-i-split-divide-polyline-shapefiles-into-equally-length-smaller-segments

#Creating 30 meter segments
#https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html
#https://github.com/r-spatial/lwgeom/issues/16

## Reading in air quality data from ropenaq
https://openaq.org/#/location/Bronx%20-%20IS74?_k=0z7fjl
{r}
require(ropenaq)
#bronx_is74<-aq_measurements( country="US", city="New York-Northern New Jersey-Long Island", location="Bronx - IS74", date_from = '2018-12-31', date_to = '2020-02-06', limit = 10000)
bronx_is74<-read.csv("/Users/sanjanapaul/Downloads/bronxis74.csv")
bronx_is74$date<-as.character(paste(bronx_is74$local))
bronx_is74$date<-paste(substr(bronx_is74$date,1,10 ),substr(bronx_is74$date, 12, 19))
bronx_is74$date<-as.POSIXct(bronx_is74$date, "%Y-%m-%d %H:%M:%S")

bronxis74_pm25<-subset(bronx_is74, parameter="pm25")
bronxis74_pm25<-subset(bronxis74_pm25, bronxis74_pm25$date>="2020-01-20 18:00:00" & bronxis74_pm25$date <="2020-02-14 15:00:00")
summary(bronxis74_pm25$value)


#Calculating background air pollution concentrations
{r}
require(rlist)
#Calculating background
#Method 1: For each day take the lower percentile of the values (using 10%)
#Splitting the data by date and hour of the day
out<-split(aq, f=list(aq$date, aq$h))
require(rlist)
try<-c()
for(j in 1:length(out)){
  if(nrow(out[[j]])==0){
    try<-list.append(try, FALSE)
  } else { try<-list.append(try, TRUE)}
}
cnames<-c("date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")
#Subsetting the dataframe
aq_sub<-aq[cnames]
#Initializing the background percentile for Day 1
bg_percentile<-aq_sub[1,]
for (i in 1:length(out)){
  #print(i)
  if(nrow(out[[i]])==0){next}
  else{
    #print(i)
  bg_percentile[i,]<-(aggregate(cbind(n1,n3,n10,bin0, bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12, bin13, bin14, bin15, bin16, bin17, bin18, bin19, bin20, bin21, bin22, bin23,  pm1, pm25, pm10)~date, data=out[[i]], FUN=function(x){quantile(x, 0.1)}))}
  bg_percentile$hour[i]<-unique(out[[i]]$h)[1]
  }'''

print("Calculating lower 10 percentile background")


#Method 2: Background from the NYC Bronx reference air quality monitoring station
hourly_PMdata<-bronxis74_pm25
hourly_PMdata$date<-as.character(paste0(substr(bronxis74_pm25$local, 1, 10), " ", substr(bronxis74_pm25$local, 12, 19)))
hourly_PMdata$date<-as.POSIXct(hourly_PMdata$date, format="%Y-%m-%d %H:%M:%S")
hourly_PMdata<-hourly_PMdata[order(hourly_PMdata$date),]
hourly_PMdata$h<-as.integer(substr(hourly_PMdata$date, 12, 13))
hourly_PMdata<-subset(hourly_PMdata,select=c(date, value, h) )
hourly_PMdata<-subset(hourly_PMdata, hourly_PMdata$value>=0)


#Method 3: Spline of minimums
#Calculate index corresponding to every 30 sec intervals
#Subsetting to the required number of columns
cnames<-c("time","date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")
aq_sub<-aq[cnames]
#Filtering out points > 100 micrograms/meter cube
aq_sub<-subset(aq_sub, aq_sub$PM2.5<50)
#Initializing the 30 second averages for each day
aq30sec<-list(rep(aq_sub[1,], length(unique(aq_sub$date))))

out<-split(aq, f=aq$date)
for(i in (1:length(out))){
cnames<-c("time","date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")
aq_sub<-out[[i]][cnames]
  aq30sec[[i]]<- aq_sub %>% group_by(time =cut(time, breaks="30 sec")) %>% summarize_all(funs(mean), na.rm=TRUE)
}
print("Calculated average every 30 seconds")

#Finding minimum every 10 minutes
cnames<-c("time","date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")

aq_sub<-aq[cnames]
#Initializing 10 minute aq minimums for each day
aq10min<-list(rep(aq_sub[1,], length(unique(aq_sub$date))))
for(i in (1:length(out))){
cnames<-c("time","date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")

  aq_sub<-aq30sec[[i]][cnames]
  aq_sub$time<-as.POSIXct(aq_sub$time,  origin="1970-01-01 00:00")
  aq10min[[i]]<- aq_sub %>% group_by(time = cut(time, breaks="10 min")) %>% summarize_all(funs(min), na.rm=TRUE)
}
print("Calculated min every 10 minutes")


#Fitting a spline
#k=15 seems to be the maximum number of dimensions to construct a smooth spline fit
out<-split(aq, f=aq$date)
for(i in (1:length(out))){
 cnames<-c("time","date","n1", "n3","n10", "bin0","bin1", "bin2","bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15","bin16", "bin17", "bin18", "bin19", "bin20", "bin21", "bin22", "bin23", "pm1", "pm25", "pm10")

  aq_sub<-aq10min[[i]][cnames]
  aq_sub$time<-as.POSIXct(aq_sub$time,  origin="1970-01-01 00:00")
  aq_sub$time<-as.integer(aq_sub$time)

  out[[i]]$time<-as.integer(as.POSIXct(out[[i]]$time, format="%H:%M:%S"))
  if(i!=29){
    n=3

  bn1<-gam(n1 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgn1<-predict.gam(bn1, data.frame(time=out[[i]]$time))
  print("Predicted for n1")

  bn3<-gam(n3 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgn3<-predict.gam(bn3, data.frame(time=out[[i]]$time))
  print("Predicted for n3")

  bn10<-gam(n10 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgn10<-predict.gam(bn10, data.frame(time=out[[i]]$time))
  print("Predicted for n10")

bpm1<-gam(pm1 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgpm1<-predict.gam(bpm1, data.frame(time=out[[i]]$time))
  print("Predicted for PM1")

  bpm25<-gam(pm25 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgpm25<-predict.gam(bpm25, data.frame(time=out[[i]]$time))
  print("Predicted for PM2.5")

  bpm10<-gam(pm10 ~ s(time, bs = "ts", k= n), data=aq_sub)
  out[[i]]$bgpm10<-predict.gam(bpm10, data.frame(time=out[[i]]$time))
  print("Predicted for PM10")
  print(i)
}

}
print("Finished prediction of background air pollution")
#rm(bpm1, bpm25, bpm10, bn1, bn3, bn10, bbin0, bbin1, bbin2, bbin3, bbin4, bbin5, bbin6, bbin7, bbin8, bbin9, bbin10, bbin11, bbin12, bbin13, bbin14, bbin15)
# Combine all the files
out_merge<-out[[1]]
for (i in 2:length(out)){
  if(ncol(out[[i]])!=ncol(out_merge)){next}
  print(i)
  out_merge<-rbind(out_merge, out[[i]])
}
#rm(out)
#Converting all the negative background air pollutants to 0

#out_merge[out_merge<0]<-0
out_merge$Longitude<-out_merge$long
out_merge$Latitude<-out_merge$lat
```
#Comparing background values
```{r}
#Plotting the three values of background air pollution
temp<-subset(out_merge, select=c(bgpm1, bgpm25, bgpm10, bgn1, bgn3, bgn10,  date, h, time))
temp$date<-as.Date(temp$date, "%Y-%m-%d")
#Check time zone
temp$time<-as.POSIXct(temp$time,  origin="1970-01-01 00:00")
#Creating time
#temp$time<-as.POSIXct(aq$time,"%H:%M:%S" )

hourly_PMdata$date <- as.Date(hourly_PMdata$date, "%Y-%m-%d")

background<-merge(temp, hourly_PMdata, by=c("date", "h"), all.x=TRUE)
rm(temp)

background<-background[!duplicated(background),]

bg_percentile$date<-as.Date(bg_percentile$date, format="%Y-%m-%d")
#bg_percentile$date<-format(bg_percentile$date, "%Y-%m-%d")
background<-merge(background, bg_percentile, by.x=c("date", "h"), by.y=c("date", "hour"), all.x=TRUE)
background<-background[!duplicated(background),]

background$bgn1<-as.numeric(background$bgn1)
background$bgn3<-as.numeric(background$bgn3)
background$bgn10<-as.numeric(background$bgn10)

background$id<-seq.int(nrow(background))
background$bgpm25<-as.numeric(background$bgpm25)
background$value<-as.numeric(background$value)
background$pm25<-as.numeric(background$pm25)

a<-ggplot(background, aes(x=id) )+geom_line(aes(y=bgpm25, color="Spline of minimums"))+geom_line(aes(y=pm25, color="Lowest tenth Percentile"))+
  geom_line(aes(y=value, color="Reference monitor"))+xlab("")+ylab("Background (micrograms/cubic meter")+scale_color_manual(name="Methods", values=c("Spline of minimums"="red", "Lowest tenth Percentile"="blue", "Reference monitor"="yellow"))
a

b<-ggplot(background, aes(x=id) )+geom_line(aes(y=bgn1, color="Spline of minimums"))+geom_line(aes(y=n1, color="Lowest tenth Percentile"))+scale_color_manual(name="Methods", values=c("Spline of minimums"="red", "Lowest tenth Percentile"="blue"))+xlab("")+ylab("Background N1")
b

c<-ggplot(background, aes(x=id) )+geom_line(aes(y=bgn3, color="Spline of minimums"))+geom_line(aes(y=n3, color="Lowest tenth Percentile"))+scale_color_manual(name="Methods", values=c("Spline of minimums"="red", "Lowest tenth Percentile"="blue"))+xlab("")+ylab("Background N3")
c

d<-ggplot(background, aes(x=id) )+geom_line(aes(y=bgn10, color="Spline of minimums"))+geom_line(aes(y=n10, color="Lowest tenth Percentile"))+scale_color_manual(name="Methods", values=c("Spline of minimums"="red", "Lowest tenth Percentile"="blue"))+xlab("")+ylab("Background N10")
d

#timePlot(background, pollutant=c("bgpm25","pm25", "value" ), date.pad=TRUE, group=TRUE)


```

#Correcting for background pollution
```{r}
#Correcting for background air pollution
background$bgpm25_ref<-background$value

background$bgpm1_percentile<-background$pm1
background$bgpm25_percentile<-background$pm25
background$bgpm10_percentile<-background$pm10

background$bgn1_percentile<-background$n1
background$bgn3_percentile<-background$n3
background$bgn10_percentile<-background$n10

temp<-subset(background, select=c(time, bgpm25_ref, bgpm25_percentile, bgpm1_percentile, bgpm10_percentile, bgn1_percentile, bgn3_percentile, bgn10_percentile))
out_merge$time<-as.POSIXct(out_merge$time,  origin="1970-01-01 00:00")
out_merge<-merge(out_merge, temp, by="time", all.x=TRUE)

#out_merge$bgpm25_ref<-background$value
#out_merge$bgpm25_percentile<-background$pm25


meanpollutants<-aggregate(cbind(bgpm1, bgpm25, bgpm10, bgn1, bgn3, bgn10,  bgpm25_ref, bgpm25_percentile, bgn1_percentile, bgn3_percentile, bgn10_percentile, bgpm1_percentile, bgpm10_percentile)~ date , out_merge, FUN=function(x){mean(x)})
colnames(meanpollutants)<-c("date","meanbgpm1", "meanbgpm25", "meanbgpm10","meanbgn1", "meanbgn3", "meanbgn10", "meanbgpm25_ref", "meanbgpm25_percentile", "meanbgn1_percentile", "meanbgn3_percentile", "meanbgn10_percentile", "meanbgpm1_percentile", "meanbgpm10_percentile" )

#medianpollutants<-aggregate(cbind(bgpm1, bgpm25, bgpm10, bgn1, bgn3, bgn10,  bgpm25_ref, bgpm25_percentile)~ date , out_merge, FUN=function(x){median(x)})
#colnames(medianpollutants)<-c("date","medianbgpm1", "medianbgpm25", "medianbgpm10","medianbgn1", "medianbgn3", "medianbgn10", "medianbgpm25_ref", "medianbgpm25_percentile" )

out_merge<-merge(out_merge, meanpollutants,
                 by="date", all=TRUE)
#out_merge<-merge(out_merge, medianpollutants, by="date", all=TRUE)
rm(meanpollutants)

out_merge$spm1<-ifelse(out_merge$bgpm1<= out_merge$pm1, out_merge$pm1-out_merge$bgpm1 + out_merge$meanbgpm1, out_merge$pm1*out_merge$meanbgpm1/out_merge$bgpm1)

out_merge$spm25<-ifelse(out_merge$bgpm25<= out_merge$pm25, out_merge$pm25-out_merge$bgpm25+ out_merge$meanbgpm25, out_merge$pm25*out_merge$meanbgpm25/out_merge$bgpm25)

out_merge$spm10<-ifelse(out_merge$bgpm10<= out_merge$pm10, out_merge$pm10-out_merge$bgpm10+ out_merge$meanbgpm10, out_merge$pm10*out_merge$meanbgpm10/out_merge$bgpm10)

out_merge$sn1<-ifelse(out_merge$bgn1<= out_merge$n1, out_merge$n1-out_merge$n1 + out_merge$meanbgn1, out_merge$n1*out_merge$meanbgn1/out_merge$bgn1)
out_merge$sn3<-ifelse(out_merge$bgn3<= out_merge$n3, out_merge$n3-out_merge$bgn3 + out_merge$meanbgn3, out_merge$n3*out_merge$meanbgn3/out_merge$bgn3)
out_merge$sn10<-ifelse(out_merge$bgn10<= out_merge$n10, out_merge$n10-out_merge$bgn10  + out_merge$meanbgn10, out_merge$n10*out_merge$meanbgn10/out_merge$bgn10)

#Correcting usingreference monitor
out_merge$spm25_ref<-ifelse(out_merge$bgpm25_ref<= out_merge$pm25, out_merge$pm25-out_merge$bgpm25_ref + out_merge$meanbgpm25_ref, out_merge$pm25*out_merge$meanbgpm25_ref/out_merge$bgpm25_ref)

#Correcting using percentile method
out_merge$spm25_percentile<-ifelse(out_merge$bgpm25_percentile<= out_merge$pm25, out_merge$pm25-out_merge$bgpm25_percentile+ out_merge$meanbgpm25_percentile, out_merge$pm25*out_merge$meanbgpm25_percentile/out_merge$bgpm25_percentile)

out_merge$spm10_percentile<-ifelse(out_merge$bgpm10_percentile<= out_merge$pm10, out_merge$pm10-out_merge$bgpm10_percentile+ out_merge$meanbgpm10_percentile, out_merge$pm10*out_merge$meanbgpm10_percentile/out_merge$bgpm10_percentile)

out_merge$spm1_percentile<-ifelse(out_merge$bgpm1_percentile<= out_merge$pm1, out_merge$pm1-out_merge$bgpm1_percentile+ out_merge$meanbgpm1_percentile, out_merge$pm1*out_merge$meanbgpm1_percentile/out_merge$bgpm1_percentile)

out_merge$sn1_percentile<-ifelse(out_merge$bgn1_percentile<= out_merge$n1, out_merge$n1-out_merge$bgn1_percentile+ out_merge$meanbgn1_percentile, out_merge$n1*out_merge$meanbgn1_percentile/out_merge$bgn1_percentile)

out_merge$sn3_percentile<-ifelse(out_merge$bgn3_percentile<= out_merge$n3, out_merge$n3-out_merge$bgn3_percentile+ out_merge$meanbgn3_percentile, out_merge$n3*out_merge$meanbgn3_percentile/out_merge$bgn3_percentile)

out_merge$sn10_percentile<-ifelse(out_merge$bgn10_percentile<= out_merge$n10, out_merge$n10-out_merge$bgn10_percentile+ out_merge$meanbgn10_percentile, out_merge$n10*out_merge$meanbgn10_percentile/out_merge$bgn10_percentile)

out_merge$id<-seq.int(nrow(out_merge))
a1<-ggplot(out_merge, aes(x=id) )+geom_line(aes(y=spm25_ref, color="Reference monitor"))+geom_line(aes(y=spm25, color="Spline of minimums"))+geom_line(aes(y=spm25_percentile, color="Lowest tenth Percentile"), position=position_jitter(w=0.02, h=0))+xlab("")+ylab("Background-corrected PM2.5")+scale_color_manual(name="Methods", values=c("Spline of minimums"="red", "Lowest tenth Percentile"="blue", "Reference monitor"="yellow"))
a1

#Checking difference between measured and corrected PM2.5
b1<-ifelse(out_merge$pm25>0, (out_merge$pm25-out_merge$spm25)/out_merge$pm25, 0)
print(summary(b1[is.finite(b1)]))
#Mean difference is -2.1%. Median diff = 0

b2<-ifelse(out_merge$pm25>0, (out_merge$pm25-out_merge$spm25_ref)/out_merge$pm25, 0)
print(summary(b2[is.finite(b2)]))
#Mean diff is -5.5%, Median = -1.9%

b3<-ifelse(out_merge$pm25>0, (out_merge$pm25-out_merge$spm25_percentile)/out_merge$pm25, 0)
print(summary(b3[is.finite(b3)]))
#Mean = -14.6% Median = -1.1%

 b4<-(out_merge$spm25_ref - out_merge$spm25)/out_merge$spm25
print(summary(b4[is.finite(b4)]))
#Mean is 4.9% Median = 1.1%

 b5<-(out_merge$spm25_percentile - out_merge$spm25)/out_merge$spm25
print(summary(b5[is.finite(b5)]))
#Mean is 13.1%, Median = 0%

#Checking the difference for N1
b6<-ifelse(out_merge$n1>0, (out_merge$n1-out_merge$sn1)/out_merge$n1, 0)
print(summary(b6[is.finite(b6)]))
#Mean = 29.9%
#Median = 32.7%

b7<-ifelse(out_merge$n1>0, (out_merge$n1-out_merge$sn1_percentile)/out_merge$n1, 0)
print(summary(b7[is.finite(b7)]))
#Mean = -9.9%
#Median =-1.5%

b8<-ifelse(out_merge$n1>0, (out_merge$sn1-out_merge$sn1_percentile)/out_merge$sn1, 0)
print(summary(b8[is.finite(b8)]))
#Mean=-129.7%
#Median=-57.1%

#Checking the difference for N10
b9<-ifelse(out_merge$n10>0, (out_merge$n10-out_merge$sn10)/out_merge$n10, 0)
print(summary(b9[is.finite(b9)]))
#Mean = -32.7%
#Median = 0%

b10<-ifelse(out_merge$n10>0, (out_merge$n10-out_merge$sn10_percentile)/out_merge$n10, 0)
print(summary(b10[is.finite(b10)]))
#Mean = -0.7%
#Median = 0%

b11<-ifelse(out_merge$n10>0, (out_merge$sn10-out_merge$sn10_percentile)/out_merge$sn10, 0)
print(summary(b11[is.finite(b11)]))
#Mean =-15.9%
#Median =0%

#Checking the difference for N3
b12<-ifelse(out_merge$n3>0, (out_merge$n3-out_merge$sn3)/out_merge$n3, 0)
print(summary(b12[is.finite(b12)]))
#Mean = -13.9%
#Median = 0%

b13<-ifelse(out_merge$n3>0, (out_merge$n3-out_merge$sn3_percentile)/out_merge$n3, 0)
print(summary(b13[is.finite(b13)]))
#Mean = -0.5%
#Median = 0%

b14<-ifelse(out_merge$n3>0, (out_merge$sn3-out_merge$sn3_percentile)/out_merge$sn3, 0)
print(summary(b14[is.finite(b14)]))
#Mean =-3.6%
#Median =0%


#rm(b1, b2, aq10min, aq30sec)

```
#Calculate the number of stationary points
```{r}
#Calculate how much data is stationary
temp<-subset(out_merge, select=c(Latitude, Longitude))
temp<-temp[!duplicated(temp),]
print(nrow(temp))
#This is the number of unique lat and longs which means its the number ofnon-stationary data
print(nrow(temp)/nrow(aq))
```

#Calculating velocity
```{r}
require(geosphere)
out<-split(out_merge, f=out_merge$date)
for(i in 1:length(out)){
  print(i)
  out[[i]]$distance_trav<-NA
  out[[i]]$velocity<-NA
  for(j in 2:nrow(out[[i]])){

    out[[i]]$distance_trav[j]<-distm(c(out[[i]]$Longitude[j], out[[i]]$Latitude[j]),c(out[[i]]$Longitude[j-1], out[[i]]$Latitude[j-1]), distHaversine)[1]
    difftime1<-as.double(out[[i]]$time[j]-out[[i]]$time[j-1])

    out[[i]]$velocity[j]<-out[[i]]$distance_trav[j]/difftime1
 #   print(out[[i]]$velocity[j])
  }
}

#Combine to get velocity
out_merge<-out[[1]]
for (i in 2:length(out)){
  if(ncol(out[[i]])!=ncol(out_merge)){next}
  print(i)
  out_merge<-rbind(out_merge, out[[i]])
}

```

#velocity calculations
```{r}
summary(out_merge$velocity[is.finite(out_merge$velocity)])
#Average velocity when the truck is moving
summary(out_merge[is.finite(out_merge$velocity)!=0,]$velocity)
#The percentage of points < 5m/s=18 kmph
 sum(abs(out_merge$velocity)<5, na.rm=TRUE)/nrow(out_merge)
 cor(abs(out_merge$velocity), out_merge$pm25, use="complete.obs")
 no_vel<-out_merge[out_merge$velocity!=0,]
 cor(abs(no_vel$velocity), no_vel$pm25, use="complete.obs")
```

#join with raster from NYCCAS
```{r}
rm(out,  temp, hourly_PMdata, a, b, c, d, bpm25, bpm10, bn1, bn3, bn10, bg_percentile, aq10min, aq30sec, aq_sub, a1)
rm(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12,b13, b14, bpm1)
out_merge<-out_merge[!duplicated(out_merge),]
NYCCAS<-raster("/Users/sanjanapaul/Downloads/pmw10b/pmw10b")
pts<-subset(out_merge, select=c(Longitude, Latitude))
coordinates(pts) <- c("Longitude", "Latitude")
crs(pts)<-CRS("+init=epsg:4326")
NYCCAS1<-projectRaster(NYCCAS, crs=crs(pts))
NYCCAS_grid<-rasterToPoints(NYCCAS1)
head(NYCCAS_grid)
aq_NYCCAS<-raster::extract(NYCCAS1, pts, sp=TRUE)
aq_NYCCAS<-as.data.frame(aq_NYCCAS)
out_merge$pmw10b<-aq_NYCCAS$pmw10b
aq_NYCCAS_pmw10b<-aq_NYCCAS

#Joining with No2 and BC
NYCCAS_no2<-raster("/Users/sanjanapaul/Downloads/no2w10/no2w10")
pts<-subset(out_merge, select=c(Longitude, Latitude))
coordinates(pts) <- c("Longitude", "Latitude")
crs(pts)<-CRS("+init=epsg:4326")
NYCCAS1<-projectRaster(NYCCAS_no2, crs=crs(pts))
NYCCAS_grid<-rasterToPoints(NYCCAS1)
head(NYCCAS_grid)
aq_NYCCAS<-raster::extract(NYCCAS1, pts, sp=TRUE)
aq_NYCCAS<-as.data.frame(aq_NYCCAS)
out_merge$no2<-aq_NYCCAS$no2w10
aq_nyccas_no2<-aq_NYCCAS

NYCCAS_bc<-raster("/Users/sanjanapaul/Downloads/no2w10/bcw10")
pts<-subset(out_merge, select=c(Longitude, Latitude))
coordinates(pts) <- c("Longitude", "Latitude")
crs(pts)<-CRS("+init=epsg:4326")
NYCCAS1<-projectRaster(NYCCAS_bc, crs=crs(pts))
NYCCAS_grid<-rasterToPoints(NYCCAS1)
head(NYCCAS_grid)
aq_NYCCAS<-raster::extract(NYCCAS1, pts, sp=TRUE)
aq_NYCCAS<-as.data.frame(aq_NYCCAS)
out_merge$bc<-aq_NYCCAS$bcw10
```

#Finding NYCCAS grid cell
(Alternative for the future:
https://stackoverflow.com/questions/34242154/merging-two-data-frames-both-with-coordinates-based-on-the-closest-location
)
```{r}
require(sf)
df1<-subset(out_merge, select=c(Longitude, Latitude))
colnames(df1)<-c("X", "Y")
my_spdf<- SpatialPoints(df1, proj4string=CRS("+init=epsg:4326"))

NYCCAS_grid<-as.data.frame(NYCCAS_grid)
df2<-subset(NYCCAS_grid, select=c(x,y))
colnames(df2)<-c("X", "Y")
my_spdf2<- SpatialPoints(df2, proj4string=CRS("+init=epsg:4326"))
my_spdf1<- SpatialPoints(df1, proj4string=CRS("+init=epsg:4326"))
my_sfdf1 = st_sf(id_pt = 1:length(my_spdf1), geom = st_as_sfc(my_spdf1))
my_sfdf2 = st_sf(id_pt = 1:length(my_spdf2), geom = st_as_sfc(my_spdf2))

shapefile_roads_sf_w_pts = st_join(my_sfdf1,
                                    my_sfdf2,
                                    join = st_nearest_feature)
colnames(shapefile_roads_sf_w_pts)<-c("id.x", "ID", "geom")
```

#Joining out_merge with NYCCAS
```{r}
out_merge$id.x<-seq(1:nrow(out_merge))
out_merge<-merge(out_merge, shapefile_roads_sf_w_pts)
```

#The heavy lifting
```{r}
temp<-subset(out_merge, select=c(ID, pmw10b, no2, bc))
temp<-temp[!duplicated(temp),]
out30<-out_merge
out30$time<-as.POSIXct(out30$time, format="%Y-%m-%d %H:%M:%S")

#Find number of days for each segment
days30 <- aggregate(cbind(time, date) ~ ID, data=out30, FUN=function(x){length(unique(x))} )

colnames(days30)<-c("ID", "number_obs", "number_days")

#Plot number of days each segment is sampled for
day30<-aggregate(cbind(time, date)~ID, data=out30, FUN=function(x){length(unique(x))})
colnames(day30)<-c("ID", "number_obs", "number_days")
hist(day30$number_days, nclass=24, xlab="Number of days", ylab="Frequency", main="")

#Aggregateby road segment median
segment_median30<-aggregate(cbind(sn1, sn1_percentile, sn3, sn3_percentile, sn10, sn10_percentile, spm1, spm25, spm25_percentile, spm25_ref, spm10,  pm1, pm25, pm10, n1, n3, n10)~ID, data=out30,FUN=function(x){median(x, na.rm=TRUE)})

colnames(segment_median30)<-c("ID", "sn1median","sn1median_percentile", "sn3median", "sn3median_percentile", "sn10median","sn10median_percentile", "spm1median", "spm25median", "spm25median_percentile", "spm25median_ref", "spm10median",  "pm1median", "pm25median", "pm10median", "n1median", "n3median", "n10median")

#Aggregateby road segment mean
segment_mean30<-aggregate(cbind(sn1, sn1_percentile, sn3, sn3_percentile, sn10, sn10_percentile, spm1, spm25, spm25_percentile, spm25_ref, spm10,  pm1, pm25, pm10, n1, n3, n10)~ID, data=out30,FUN=function(x){mean(x, na.rm=TRUE)})

colnames(segment_mean30)<-c("ID", "sn1mean","sn1mean_percentile", "sn3mean", "sn3mean_percentile", "sn10mean","sn10mean_percentile", "spm1mean", "spm25mean", "spm25mean_percentile", "spm25mean_ref", "spm10mean",  "pm1mean", "pm25mean", "pm10mean", "n1mean", "n3mean", "n10mean")

#Aggregateby road segment var
segment_var30<-aggregate(cbind(sn1, sn1_percentile, sn3, sn3_percentile, sn10, sn10_percentile, spm1, spm25, spm25_percentile, spm25_ref, spm10,  pm1, pm25, pm10, n1, n3, n10)~ID, data=out30,FUN=function(x){var(x, na.rm=TRUE)})

colnames(segment_var30)<-c("ID", "sn1var","sn1var_percentile", "sn3var", "sn3var_percentile", "sn10var","sn10var_percentile", "spm1var", "spm25var", "spm25var_percentile", "spm25var_ref", "spm10var",  "pm1var", "pm25var", "pm10var", "n1var", "n3var", "n10var")

#Aggregateby road segment skewness
segment_skewness30<-aggregate(cbind(sn1, sn1_percentile, sn3, sn3_percentile, sn10, sn10_percentile, spm1, spm25, spm25_percentile, spm25_ref, spm10,  pm1, pm25, pm10, n1, n3, n10)~ID, data=out30,FUN=function(x){skewness(x, na.rm=TRUE)})

colnames(segment_skewness30)<-c("ID", "sn1skew","sn1skew_percentile", "sn3skew", "sn3skew_percentile", "sn10skew","sn10skew_percentile", "spm1skew", "spm25skew", "spm25skew_percentile", "spm25skew_ref", "spm10skew",  "pm1skew", "pm25skew", "pm10skew", "n1skew", "n3skew", "n10skew")

#Combine mean and median
segment30<-merge(segment_median30, segment_mean30, by=c("ID"))
#Combine with skewness
segment30<-merge(segment30, segment_skewness30, by=c( "ID"))
#combine with number of days
segment30<-merge(segment30, days30, by="ID")
#Combine with variance
segment30<-merge(segment30, segment_var30, by="ID")

temp<-subset(out_merge, select=c(ID, pmw10b, no2, bc))
temp<-temp[order(temp$pmw10b),]
temp<-temp[!duplicated(temp$ID),]
temp<-temp[temp$ID %in% segment30$ID,]
segment30<-merge(segment30, temp, by="ID", all.x=TRUE)

cor(segment30$spm25median, segment30$pmw10b, use="pairwise.complete.obs")

cor(segment30$spm25median, segment30$no2, use="pairwise.complete.obs")

cor(segment30$spm25median, segment30$bc, use="pairwise.complete.obs")
#Joining with latitude and longitude of x and y
NYCCAS_grid$ID<-seq(1:nrow(NYCCAS_grid))
temp<-NYCCAS_grid
temp<-temp[!duplicated(temp),]
temp<-temp[temp$ID %in% segment30$ID,]
segment30<-merge(segment30, temp, by="ID", all.x=TRUE)
write.csv(segment30, "segment30_nyc3.csv")

```

#Calculating variance in median
```{r}
#Use bootstrap to calculate the  mad
#PM25

out30$join_ID<-out30$ID

segment30$madpm25<- -10
segment30$madpm25_ref<- -10
segment30$madpm25_percentile<- -10

#PM1
segment30$madpm1<- -10

#PM10
segment30$madpm10<- -10

#n1
segment30$madn1<- -10
segment30$madn1_percentile<- -10

#n3

segment30$madn3<- -10
segment30$madn3_percentile<- -10

#n10

segment30$madn10<- -10
segment30$madn10_percentile<- -10

#uncorrected data
segment30$madpm25_u<- -10
segment30$madn1_u<- -10
segment30$madn3_u<- -10
segment30$madn10_u<- -10

#30 meter
#fn1n3n10

#segment30$madfn1n10<- -10
#segment30$madsfns1n10_percentile<- -10
#segment30$madsfn1sn10<- -10



samplemean <- function(x, d) {
  return(mean(x[d]))
}

samplemedian <- function(x, d) {
  return(median(x[d]))
}

#print(Sys.time())
for(i in 1:nrow(segment30)){
  print(i)
  sub<-out30[out30$ID %in% segment30$ID[i],]
  # sub<-subset(sub, select=spm25)


  bootmedianpm25<-boot(data=sub$spm25, statistic=samplemedian, R=1000)
  segment30$bootmadpm25[i]<-sd(bootmedianpm25$t)

  #bootmedianpm25_ref<-boot(data=sub$spm25_ref, statistic=samplemedian, R=1000)
  #segment30$bootmadpm25_ref[i]<-sd(bootmedianpm25_ref$t)

  bootmedianpm25_percentile<-boot(data=sub$spm25_percentile, statistic=samplemedian, R=1000)
  segment30$bootmadpm25_percentile[i]<-sd(bootmedianpm25_percentile$t)

  bootmedianpm10<-boot(data=sub$spm10, statistic=samplemedian, R=1000)
   segment30$bootmadpm10[i]<-sd(bootmedianpm10$t)

  bootmedianpm1<-boot(data=sub$spm1, statistic=samplemedian, R=1000)
  segment30$bootmadpm1[i]<-sd(bootmedianpm1$t)

  bootmediann1<-boot(data=sub$sn1, statistic=samplemedian, R=1000)
  segment30$bootmadn1[i]<-sd(bootmediann1$t)

  bootmediann3<-boot(data=sub$sn3, statistic=samplemedian, R=1000)
  segment30$bootmadn3[i]<-sd(bootmediann3$t)

  bootmediann10<-boot(data=sub$sn10, statistic=samplemedian, R=1000)
  segment30$bootmadn10[i]<-sd(bootmediann10$t)


  bootmediann1_percentile<-boot(data=sub$sn1_percentile, statistic=samplemedian, R=1000)
  bootmediann3_percentile<-boot(data=sub$sn3_percentile, statistic=samplemedian, R=1000)
  bootmediann10_percentile<-boot(data=sub$sn10_percentile, statistic=samplemedian, R=1000)

  segment30$bootmadn3_percentile[i]<-sd(bootmediann3_percentile$t)
  segment30$bootmadn1_percentile[i]<-sd(bootmediann1_percentile$t)
  segment30$bootmadn10_percentile[i]<-sd(bootmediann10_percentile$t)

  bootmedianpm25_u<-boot(data=sub$pm25, statistic=samplemedian, R=1000)
  segment30$bootmadpm25_u[i]<-sd(bootmedianpm25_u$t)

  bootmedianpm10_u<-boot(data=sub$pm10, statistic=samplemedian, R=1000)
  segment30$bootmadpm10_u[i]<-sd(bootmedianpm10_u$t)

  bootmedianpm1_u<-boot(data=sub$pm1, statistic=samplemedian, R=1000)
  segment30$bootmadpm1_u[i]<-sd(bootmedianpm1_u$t)

  bootmediann1_u<-boot(data=sub$n1, statistic=samplemedian, R=1000)
  segment30$bootmadn1_u[i]<-sd(bootmediann1_u$t)

  bootmediann3_u<-boot(data=sub$n3, statistic=samplemedian, R=1000)
  segment30$bootmadn3_u[i]<-sd(bootmediann3_u$t)

  bootmediann10_u<-boot(data=sub$n10, statistic=samplemedian, R=1000)
  segment30$bootmadn10_u[i]<-sd(bootmediann10_u$t)
}

#Final file
write.csv(segment30, "segment30_nyc2.csv")
```

#Finishing the calculation
```{r}
segment30<-segment30[!duplicated(segment30),]
#Calculating normalised errors
#30m
segment30$normmadpm25<-ifelse(segment30$spm25median!=0,segment30$bootmadpm25/segment30$spm25median, 0)
segment30$normvarpm25<-ifelse(segment30$spm25mean!=0, segment30$spm25var/segment30$spm25mean, 0)

segment30$normmadpm25_percentile<-ifelse(segment30$spm25median_percentile!=0, segment30$bootmadpm25_percentile/segment30$spm25median_percentile, 0)

segment30$normvarpm25_percentile<-ifelse(segment30$spm25mean_percentile!=0, segment30$spm25var_percentile/segment30$spm25mean_percentile, 0)

segment30$normmadpm1<-ifelse(segment30$spm1median!=0, segment30$bootmadpm1/segment30$spm1median, 0)
segment30$normvarpm1<-ifelse(segment30$spm1mean!=0, segment30$spm1var/segment30$spm1mean, 0)

segment30$normmadpm10<-ifelse(segment30$spm10median!=0, segment30$bootmadpm10/segment30$spm10median, 0)
segment30$normvarpm10<-ifelse(segment30$spm10mean!=0, segment30$spm10var/segment30$spm10mean, 0)

segment30$normmadn1<-ifelse(segment30$sn1median!=0, segment30$bootmadn1/segment30$sn1median, 0)
segment30$normvarn1<-ifelse(segment30$sn1mean!=0, segment30$sn1var/segment30$sn1mean, 0)

segment30$normmadn3<-ifelse(segment30$sn3median!=0, segment30$madn3/segment30$sn3median, 0)
segment30$normvarn3<-ifelse(segment30$sn3mean!=0, segment30$sn3var/segment30$sn3mean, 0)

segment30$normmadn10<-ifelse(segment30$sn10median!=0, segment30$bootmadn10/segment30$sn10median, 0)
segment30$normvarn10<-ifelse(segment30$sn10mean!=0, segment30$sn10var/segment30$sn10mean, 0)

segment30$normmadn1_percentile<-ifelse(segment30$sn1median_percentile!=0, segment30$bootmadn1_percentile/segment30$sn1median_percentile, 0)

segment30$normvarn1_percentile<-ifelse(segment30$sn1mean_percentile!=0, segment30$sn1var_percentile/segment30$sn1mean_percentile, 0)

segment30$normmadn3_percentile<-ifelse(segment30$sn3median_percentile!=0, segment30$bootmadn3_percentile/segment30$sn3median_percentile, 0)
segment30$normvarn3_percentile<-ifelse(segment30$sn3mean_percentile!=0, segment30$sn3var_percentile/segment30$sn3mean_percentile, 0)

segment30$normmadn10_percentile<-ifelse(segment30$sn10median_percentile!=0, segment30$bootmadn10_percentile/segment30$sn10median_percentile, 0)
segment30$normvarn10_percentile<-ifelse(segment30$sn10mean_percentile!=0, segment30$sn10var_percentile/segment30$sn10mean_percentile, 0)


segment30$normmadpm25_u<-ifelse(segment30$pm25median!=0, segment30$bootmadpm25_u/segment30$pm25median, 0)
segment30$normvarpm25_u<-ifelse(segment30$pm25mean!=0, segment30$pm25var/segment30$pm25mean, 0)

#Write the csv sheets
require(data.table)
write.csv(segment30, file="segment30_nyc3_complete.csv")

#Find percentage of segments for which mean error in normalise dmedian <=0.2
a<-segment30[is.finite(segment30$normmadpm25),]
sum(a$normmadpm25<=0.2)/nrow(a)
sum(a$normmadpm25<=0.2 & a$number_days>1)/nrow(a)

a1<-segment30[is.finite(segment30$normmadpm25_u),]
sum(a1$normmadpm25_u<=0.2)/nrow(a1)
sum(a1$normmadpm25_u<=0.2 & a1$number_days>1)/nrow(a1)

b<-segment30[is.finite(segment30$normmadn1),]
sum(b$normmadn1<=0.2)/nrow(b)
sum(b$normmadn1<=0.2 & b$number_days>1)/nrow(b)

c<-segment30[is.finite(segment30$normmadn3),]
sum(c$normmadn3<=0.2)/nrow(c)
sum(c$normmadn3<=0.2 & c$number_days>1)/nrow(c)

d<-segment30[is.finite(segment30$normmadn10),]
sum(d$normmadn10<=0.2)/nrow(d)
sum(d$normmadn10<=0.2 & d$number_days>1)/nrow(d)

c1<-segment30[is.finite(segment30$normmadn3_u),]
sum(c1$normmadn3_u<=0.2)/nrow(c1)
sum(c1$normmadn3_u<=0.2 & c1$number_days>1)/nrow(c1)

d1<-segment30[is.finite(segment30$normmadn10_u),]
sum(d1$normmadn10_u<=0.2)/nrow(d1)
sum(d$normmadn10_u<=0.2 & d1$number_days>1)/nrow(d1)

b<-segment30[is.finite(segment30$normmadn1_percentile),]
sum(b$normmadn1_percentile<=0.2)/nrow(b)
sum(b$normmadn1_percentile<=0.2 & b$number_days>1)/nrow(b)

c<-segment30[is.finite(segment30$normmadn3_percentile),]
sum(c$normmadn3_percentile<=0.2)/nrow(c)
sum(c$normmadn3_percentile<=0.2 & c$number_days>1)/nrow(c)

d<-segment30[is.finite(segment30$normmadn10_percentile),]
sum(d$normmadn10_percentile<=0.2)/nrow(d)
sum(d$normmadn10_percentile<=0.2 & d$number_days>1)/nrow(d)

#Checking diff background PM2.5
temp2<-(segment30$bootmedianpm25-segment30$bootmedianpm25_percentile)/segment30$bootmedianpm25_percentile
summary(temp2)

temp1<-(segment30$bootmedianpm25-segment30$bootmedianpm25_ref)/segment30$bootmedianpm25_ref
summary(temp1)
```

#Calculating zscore for NYCCAS and CityScanner datasets
```{r}
segment30$zscore_spm25median<-scale(segment30$spm25median)
segment30$zscore_spm25median_percentile<-scale(segment30$spm25median_percentile)
segment30$zscore_spm25median_ref<-scale(segment30$spm25median_ref)

segment30$zscore_spm25mean<-scale(segment30$spm25mean)
segment30$zscore_spm25mean_percentile<-scale(segment30$spm25mean_percentile)
segment30$zscore_spm25mean_ref<-scale(segment30$spm25mean_ref)

segment30$zscore_pm25mean<-scale(segment30$pm25mean)
segment30$zscore_pm25median<-scale(segment30$pm25median)

segment30$zscore_pmw10b<-scale(segment30$pmw10b)

#Final file 2

cor(segment30$zscore_spm25median, segment30$zscore_pmw10b, use="pairwise.complete.obs")
cor(segment30$spm25median, segment30$pmw10b, use="pairwise.complete.obs")


```
#Getting to this point
```{r}
int<-read.csv("/Users/sanjanapaul/Downloads/Bronx.csv")

df1<-subset(segment30, select=c(x, y))
colnames(df1)<-c("X", "Y")
my_spdf<- SpatialPoints(df1, proj4string=CRS("+init=epsg:4326"))
gc()

df2<-subset(int, select=c(Lon, Lat))
colnames(df2)<-c("X", "Y")
my_spdf2<- SpatialPoints(df2, proj4string=CRS("+init=epsg:4326"))
my_spdf1<- SpatialPoints(df1, proj4string=CRS("+init=epsg:4326"))
my_sfdf1 = st_sf(id_pt = 1:length(my_spdf1), geom = st_as_sfc(my_spdf1))
my_sfdf2 = st_sf(id_pt = 1:length(my_spdf2), geom = st_as_sfc(my_spdf2))

#rm(my_spdf1, my_spdf2, NYCCAS, NYCCAS1, NYCCAS_grid, NYCCAS_bc, NYCCAS_no2, pts, temp, segment_mean30, segment_median30, segment_skewness30, segment_var30, bronx_is74, bronxis74_pm25, aq_NYCCAS, aq_hourly, shapefile_roads_sf_w_pts)
#rm(out30, background, try, difftime1, day30, days30)
#gc()

```
#Distance
```{r}
gc()
dist = pointDistance(my_sfdf1, my_sfdf2, lonlat=TRUE, allpairs = TRUE)
colnames(dist)<-int$int_ID

```
#Limiting to intersections close by
https://stackoverflow.com/questions/31848156/delete-columns-rows-with-more-than-x-missing
```{r}
# delete columns with min>200
miss <- c()
for(i in 1:ncol(dist)) {
  if((min(dist[,i])) > 200) miss <- append(miss,i)
}
dist1 <- dist[,-miss]

dist2<-dist1
dist2[dist2>200]<-NA
#out_merge1<-cbind(out_merge, dist1)
#temp<-subset(out_merge, select=c(id.x, date, time, lat, long, Latitude, Longitude, pm1, pm25, pm10, spm1, spm25, spm10, n1, n3, n10, sn1, sn3, sn10, distance_trav, velocity, ID, pmw10b, no2, bc ))
segment30_1<-cbind(segment30, dist, dist2)
```


#Wide to long
```{r}
require(tidyverse)
require(reshape2)
require(data.table)
#out_merge2 = out_merge1 %>%
#  gather(colnames(dist2), key = intersection, value = distance, na.rm=TRUE)
segment30_2<-melt(as.data.table(segment30_1),id.var=colnames(segment30), measure.vars= colnames(dist2), na.rm=TRUE)

segment30_2<-data.table(segment30_2)

int<-data.table(int)
int$int_ID<-as.character(int$int_ID)
segment30_3<-merge(segment30_2, int, by.x="variable", by.y="int_ID", all.x=TRUE)
save.image("~/NYC_CS_uptointersectionregression.RData")
```

#Regression
```{r}
load("~/NYC_CS_uptointersectionregression.RData")
require(lme4)
segment30_3$total<-segment30_3$Incoming +segment30_3$Outgoing
segment30_3<-segment30_3[!is.na(segment30_3$variable),]
#Filter high PM2.5 values
#out_merge4<-out_merge3[out_merge3$spm25<30,]
#out_merge4<-out_merge4[!is.na(out_merge4$variable),]

v<-unique(segment30_3$variable)
a<-segment30_3[segment30_3$variable==v[2],]
plot(a$value, a$spm25)
p<-ggplot(a, aes(x=value, y=spm25median))+geom_line()
require(ggplot2)
p <- ggplot(segment30_3, aes(factor(total), spm25)) + geom_boxplot()
kruskal.test(spm25median ~ factor(total), segment30_3)
#No Significant differences

#fit<-lmer(spm25median ~ (value|total), segment30_3)
#summary(fit)
fit_n1<-lm(sn1median~ value:total, segment30_3)
summary(fit_n1)

fit_n10<-lm(sn10median~ value:total, segment30_3)
summary(fit_n10)

fit_pm25<-lm(spm25median~ value:total, segment30_3)
summary(fit_pm25)

fit_pm25<-lm(spm25median~ value:factor(total), segment30_3)
summary(fit_pm25)

fit<-lm(sn1median~ value*factor(total), segment30_3)
summary(fit)

#m1 = lmer(sn1~value+(1+value|total), data = out_merge4, REML = FALSE,  control = lmerControl( optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

#m2 = lmer(sn1~value+(value|total), data = out_merge4,  REML = FALSE, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))



```
#Effect of landuse
```{r}
landuse<-sf::st_read("/Users/priyankadesouza/Downloads/nyc_mappluto_20v6_fgdb/MapPLUTO_20v6.gdb")
pts<-st_as_sf(segment30, coords=c('x', 'y'))
st_crs(pts)<- CRS("+init=epsg:4326")
pts<-st_transform(pts, crs=crs(landuse))
segment30_landuse<-over(as(pts, "Spatial"), as(landuse, "Spatial"))
#Lots of NA's because roads are not included in PLUTO
segment30<-cbind(segment30, segment30_landuse)
a<-segment30[segment30$spm25median< 50,]
p<-ggplot(segment30, aes(x=LandUse, y=spm25median))+geom_boxplot()+ylim(0, 20)

kruskal.test(spm25median ~ factor(LandUse), segment30)
```
#CDC 500 indicators
```{r}
cdc<-st_read("/Users/priyankadesouza/Downloads/NewYork_500cities")
pts<-st_as_sf(segment30, coords=c('x', 'y'))
st_crs(pts)<- CRS("+init=epsg:4326")
pts<-st_transform(pts, crs=crs(cdc))
segment30_cdc<-over(as(pts, "Spatial"), as(cdc, "Spatial"))
segment30<-cbind(segment30, segment30_cdc)

cor(segment30$spm25median, segment30$E_POV)
fit<-lm(spm25median ~ EP_MINRTY, segment30)
fit<-lm(sn1median ~ EP_MINRTY, segment30)

#Better comparison with NYCCAS
#Disadvantaged communities
#How to do hyperlocal measurements: quantifying error between NYCCAS
#Socioeconomic variables of the city
#Average as the ground truth

#Steve Barrett
#Signal processing
#p<-ggplot(segment30, aes(x=E_POV, y=spm25median))+geom_boxplot()
```

#Creating hotspot
```{r}
pm25_100<-subset(aq, aq$pm25>100)

coordinates(pm25_100)<- c("long", "lat")
proj4string(pm25_100)<-CRS("+init=epsg:4326")
require(geosphere)
mdist <- distm(pm25_100)
hc <- hclust(as.dist(mdist), method="complete")

# define the distance threshold, in this case 100 m
d=100
# define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
pm25_100$clust <- cutree(hc, h=d)
cent <- matrix(ncol=2, nrow=max(pm25_100$clust))
require(rgeos)
for (i in 1:max(pm25_100$clust))
    # gCentroid from the rgeos package
    cent[i,] <- gCentroid(subset(pm25_100, clust == i))@coords

cent<-as.data.frame(cent)
colnames(cent)<-c("Longitude", "Latitude")

a<-pm25_100$clust
a<-as.data.frame(a)

n<-a %>%
  dplyr::group_by(a) %>%
  dplyr::summarise(n = n())

n_date<-pm25_100@data %>%
  dplyr::group_by(clust) %>%
  dplyr::summarise(length(unique(date)))

colnames(n_date)<-c("clust", "days")

n<-a %>%
    dplyr::group_by(a) %>%
    dplyr::summarise(n = n())

cent$n<-n$n
cent$cluster<-n$a
cent$n_date<-n_date$days
write.csv(cent, "cent.csv")


```

#Not run yet
#Joining with Roads

#Calculate the number of times a road has been sampled
```{r}
road_day<-aggregate(date~join_stname_lab, data=out_merge, FUN=function(x){length(unique(x))})
hist(road_day$date, xlab="Number of Roads sampled", main="")
#Average number of times a road has been sampled
mean(road_day$date)
#Number ofroads sampled >= 15 times
sum(road_day$date>=15)
#Number of roads sampled >=10 and < 15 times
sum(road_day$date>=10 & road_day$date<15)
#Number of roads sampled < 10 times
sum(road_day$date<10)
#Max number of times a road has been sampled
max(road_day$date)

```

#Clustering Analysis
#Sizes: http://www.alphasense.com/WEB1213/wp-content/uploads/2019/03/OPC-N3.pdf
```{r}
set.seed(1234)
bins<-subset(out30, select=c(bin0, bin1, bin2, bin3, bin4, bin5, bin6, bin7, bin8, bin9, bin10, bin11, bin12, bin13, bin14, bin15, bin16, bin17, bin18, bin19, bin20, bin21, bin22, bin23))


bins$dN0<-bins$bin0/log10(0.46/0.35)
bins$dN1<-bins$bin1/log10(0.66/0.46)
bins$dN2<-bins$bin2/log10(1/0.66)
bins$dN3<-bins$bin3/log10(1.3/1)
bins$dN4<-bins$bin4/log10(1.7/1.3)
bins$dN5<-bins$bin5/log10(2.3/1.7)
bins$dN6<-bins$bin6/log10(3/2.3)
bins$dN7<-bins$bin7/log10(4/3)
bins$dN8<-bins$bin8/log10(5.2/4)
bins$dN9<-bins$bin9/log10(6.5/5.2)
bins$dN10<-bins$bin10/log10(8/6.5)
bins$dN11<-bins$bin11/log10(10/8)
bins$dN12<-bins$bin12/log10(12/10)
bins$dN13<-bins$bin13/log10(14/12)
bins$dN14<-bins$bin14/log10(16/14)
bins$dN15<-bins$bin15/log10(18/16)
bins$dN16<-bins$bin16/log10(20/18)
bins$dN17<-bins$bin17/log10(22/20)
bins$dN18<-bins$bin18/log10(25/22)
bins$dN19<-bins$bin19/log10(28/25)
bins$dN20<-bins$bin20/log10(31/28)
bins$dN21<-bins$bin21/log10(34/31)
bins$dN22<-bins$bin22/log10(37/34)
bins$dN23<-bins$bin23/log10(40/37)

bins<-subset(bins, select=c(dN0, dN1, dN2, dN3, dN4, dN5,
                            dN6, dN7, dN8, dN9, dN10, dN11, dN12, dN13,dN14, dN15, dN16, dN17, dN18, dN19, dN20, dN21, dN22, dN23))
bins<-bins[!is.na(bins$dN21),]
bins<-bins[!is.na(bins$dN20),]
bins<-bins[!is.na(bins$dN19),]
bins<-bins[!is.na(bins$dN13),]
bins<-bins[!is.na(bins$dN12),]
bins<-bins[!is.na(bins$dN14),]

#Evaluation cluster tendency
#Hopkins statistic just assesses if a sample is not from a uniform distribution
#h<-hopkins(data=bins, byrow=TRUE, n=100)
#kmeans
#k<-as.data.frame(matrix(ncol=20, nrow=0))
#colnames(k)<-c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")


```
#Determine number of clusters
```{r}
wss <- (nrow(bins)-1)*sum(apply(bins,2,var))
for (i in 2:30) {
  print(i)
  wss[i] <- sum(kmeans(bins, algorithm="MacQueen", i, iter.max=1000)$withinss)
}
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
#Minimum cluster number: 16

gc()
```


```{r}
k3<-kmeans(subset(bins, select=c(dN0, dN1, dN2, dN3, dN4, dN5,
dN6, dN7, dN8, dN9, dN10, dN11, dN12, dN13,
dN14, dN15, dN16, dN17, dN18, dN19, dN20, dN21, dN22, dN23)), algorithm="MacQueen", 3, iter.max=100000)
#be careful here of Warning message QUICK-TRANSfer stage steps exceed maximum
out30<-out30[!is.na(out30$bin21),]
out30<-out30[!is.na(out30$bin20),]
out30<-out30[!is.na(out30$bin19),]
out30<-out30[!is.na(out30$bin13),]
out30<-out30[!is.na(out30$bin12),]
out30<-out30[!is.na(out30$bin14),]

bins$cluster3<-k3$cluster
print(summary(bins$cluster3))
rm(k3)
out30$cluster3<-bins$cluster3
out_merge<-out_merge[!is.na(out_merge$bin21),]
out_merge<-out_merge[!is.na(out_merge$bin20),]
out_merge<-out_merge[!is.na(out_merge$bin19),]
out_merge<-out_merge[!is.na(out_merge$bin13),]
out_merge<-out_merge[!is.na(out_merge$bin12),]
out_merge<-out_merge[!is.na(out_merge$bin14),]
out_merge$cluster3<-bins$cluster3
#out_merge1$cluster3<-bins$cluster3
print(summary(out_merge$cluster3))
```

#Plotting cluster properties
```{r}
cluster3<-aggregate(cbind(dN0, dN1, dN2, dN3, dN4, dN5, dN6, dN7,dN8, dN9, dN10, dN11, dN12, dN13,dN14, dN15, dN16, dN17,dN18, dN19, dN20, dN21, dN22, dN23) ~ cluster3, data=bins, FUN=function(x){mean(x)})
cluster3<-t(cluster3)
colnames(cluster3)<-cluster3[1,]
cluster3<-cluster3[-1,]
df <- melt(cluster3 ,  id.vars = 'cluster', variable.name = 'series')
diameter<-c(0.35, 0.46, 0.66, 1, 1.3, 1.7, 2.3, 3, 4, 5.2, 6.5, 8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 31, 34, 37)
#40
df$diameter<-mapvalues(df$Var1, from=c("dN0","dN1","dN2","dN3","dN4", "dN5", "dN6", "dN7", "dN8",
  "dN9","dN10", "dN11", "dN12", "dN13", "dN14", "dN15", "dN16","dN17", "dN18", "dN19", "dN20", "dN21", "dN22", "dN23"),
                       to=c(0.35, 0.46, 0.66, 1, 1.3, 1.7, 2.3, 3, 4, 5.2, 6.5, 8, 10, 12, 14, 16, 18, 20, 22, 25, 28, 31, 34, 37))
cluster3<-as.data.frame(cluster3)
#Printing log-normal distribution

breaksx<-c(0.38, 1, 3, 10, 20, 30, 40)
cluster3$diameter<-diameter
df$diameter<-as.numeric(paste(df$diameter))
a<-ggplot(df, aes(x=diameter, y=value, group=factor(Var2), colour=factor(Var2)))+
  geom_line()+xlab("Diameter")+ylab("dN/dlogD")+labs(colour="Cluster number")+scale_x_continuous(trans='log', breaks=breaksx, labels=breaksx)
a

bins$h<-out_merge$h
#Plot cluster counts over hours of the day

#Calculate number of times each cluster appears for each hour
c1<-sapply(split(bins,bins$h),function(x){sum(x$cluster3==1)})
c2<-sapply(split(bins,bins$h),function(x){sum(x$cluster3==2)})
c3<-sapply(split(bins,bins$h),function(x){sum(x$cluster3==3)})
ca<-rbind(c1,c2,c3)
ca<-t(ca)
ca<-as.data.frame(ca)

#Find fraction of each cluster for each measurement
ca$c1f <- ca$c1/(ca$c1+ca$c2+ca$c3)
ca$c2f <- ca$c2/(ca$c1+ca$c2+ca$c3)
ca$c3f <- ca$c3/(ca$c1+ca$c2+ca$c3)

ca<-subset(ca, select=c(c1f,c2f,c3f))

#Make sure the hour is correct
ca$hour<-seq.int(nrow(ca))
df1 <- melt(ca ,  id.vars = 'hour', variable.name = 'series')
breaksy=c(0, 0.1, 1)

#Plotting fraction of each cluster by hour
b<-ggplot(df1, aes(x=hour, y=value, group=series, colour=series) )+
  geom_line()+xlab("Hour")+ylab("fraction")+labs(colour="Cluster number")+scale_y_continuous(trans='log', labels=breaksy, breaks=breaksy)
b
```
#More cluster analysis

```{r}

out30$velocity<-out_merge$velocity
#out30$cluster5<-bins$cluster5
out30$backgroundcontribution<-ifelse(out30$pm25>0, (out30$pm25-out30$spm25)/out30$pm25, 0)

a1<-aggregate(cbind( n1, n10,  pm1, pm25, pm10, backgroundcontribution, velocity)~cluster3, data=out30, FUN=function(x){mean(x)})

a2<-aggregate(cbind( sn1,  sn10, n1, n10,  pm1, pm25, pm10, spm1, spm25, spm25_percentile, spm10, backgroundcontribution, velocity)~cluster3, data=out30, FUN=function(x){mean(x)})

a3<-aggregate(cbind( velocity)~cluster3, data=out30, FUN=function(x){sum(x==0)})

a4<-aggregate(cbind(date)~cluster3, data=out30, FUN=function(x){length(unique(x))})
table(out30$cluster3)

table(out30$cluster3)



a1
a2
a3
a4
```
#Calculating the most frequent cluster
```{r}
Mode<-function(x){
  ux<-unique(x)
  ux[which.max(tabulate(match(x,ux)))]
}
a4<-aggregate(cbind(cluster3)~ID, data=out30, FUN=function(x){Mode(x)})
summary(a4$cluster3)
print(unique(a4$cluster3))
write.csv(a4, file="most_freqcluster.csv")
#write.csv(out30, file="clusters.csv")
'''
