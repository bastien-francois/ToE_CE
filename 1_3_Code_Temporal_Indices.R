# # #### 3. Generate temporal indices
####### Periode 1850_2100 a decouper en periode de 30 ans

rm(list=ls())
set.seed(42)
library(timeDate)

#Dropping 29 of February for each year in the reference
DATES_regular_1850_2100 = atoms(timeSequence(from="1850-01-01",to="2100-12-31",by='day'))
head(DATES_regular_1850_2100)

Ind_29Feb_1850_2100 = which((DATES_regular_1850_2100[,2]==02) &
                                    (DATES_regular_1850_2100[,3]==29))
DATES_365_1850_2100 = DATES_regular_1850_2100[-Ind_29Feb_1850_2100,] #remove rows with 29 Feb.

#For Seasons (meteorological)
Ind_spring_1850_2100 =  which(DATES_365_1850_2100[,2]==03 |
                                  DATES_365_1850_2100[,2]==04 |
                                  DATES_365_1850_2100[,2]==05)
Ind_summer_1850_2100 = which(DATES_365_1850_2100[,2]==06|
                                 DATES_365_1850_2100[,2]==07|
                                 DATES_365_1850_2100[,2]==08)
Ind_automn_1850_2100 =  which(DATES_365_1850_2100[,2]==09|
                                  DATES_365_1850_2100[,2]==10 |
                                  DATES_365_1850_2100[,2]==11)
Ind_winter_1850_2100 =  which(DATES_365_1850_2100[,2]==12|
                                  DATES_365_1850_2100[,2]==01 |
                                  DATES_365_1850_2100[,2]==02)
#to del
Ind_latespringMJ_1850_2100 =  which(DATES_365_1850_2100[,2]==05 |
                                  DATES_365_1850_2100[,2]==06)
Ind_latespringAMJ_1850_2100 =  which(DATES_365_1850_2100[,2]==04 |
                                     DATES_365_1850_2100[,2]==05 |
                                     DATES_365_1850_2100[,2]==06)
Ind_latespringJ_1850_2100 =  which(DATES_365_1850_2100[,2]==06)

Ind_latesummer_1850_2100 =  which(DATES_365_1850_2100[,2]==07 |
                                  DATES_365_1850_2100[,2]==08)


# For months
Ind_months_1850_2100=list()
for(i in 1:12){
        Ind_months_1850_2100[[i]]=which(DATES_365_1850_2100[,2]==i)
}
### For sliding windows of 30 years
Ind_SlidingWindow30_winter_1850_2100=list()
Ind_SlidingWindow30_spring_1850_2100=list()
Ind_SlidingWindow30_summer_1850_2100=list()
Ind_SlidingWindow30_automn_1850_2100=list()
#to del
Ind_SlidingWindow30_latespringMJ_1850_2100=list()
Ind_SlidingWindow30_latespringAMJ_1850_2100=list()
Ind_SlidingWindow30_latespringJ_1850_2100=list()

Ind_SlidingWindow30_latesummer_1850_2100=list()


start=1850
end=1879

nb_SlidingWindow30=221
label_sliding_window=c()

nb_days_by_winter=length(Ind_winter_1850_2100)/251
nb_days_by_summer=length(Ind_summer_1850_2100)/251
nb_days_by_automn=length(Ind_automn_1850_2100)/251
nb_days_by_spring=length(Ind_spring_1850_2100)/251

nb_days_by_latespringMJ=length(Ind_latespringMJ_1850_2100)/251
nb_days_by_latespringAMJ=length(Ind_latespringAMJ_1850_2100)/251
nb_days_by_latespringJ=length(Ind_latespringJ_1850_2100)/251
nb_days_by_latesummer=length(Ind_latesummer_1850_2100)/251

### Between 1850 and 2100, 201 sliding window of size 50
for(i in 0: nb_SlidingWindow30){
  label_sliding_window=c(label_sliding_window,paste0(start+i,"_",end+i))
  print(label_sliding_window)
  Ind_SlidingWindow30_winter_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-Ind_winter_1850_2100[((i)*nb_days_by_winter+1):((i*nb_days_by_winter+30*nb_days_by_winter))]
  Ind_SlidingWindow30_automn_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-Ind_automn_1850_2100[((i)*nb_days_by_automn+1):((i*nb_days_by_automn+30*nb_days_by_automn))]
  Ind_SlidingWindow30_summer_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-Ind_summer_1850_2100[((i)*nb_days_by_summer+1):((i*nb_days_by_summer+30*nb_days_by_summer))]
  Ind_SlidingWindow30_spring_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-Ind_spring_1850_2100[((i)*nb_days_by_spring+1):((i*nb_days_by_spring+30*nb_days_by_spring))]

  Ind_SlidingWindow30_latesummer_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-
  Ind_latesummer_1850_2100[((i)*nb_days_by_latesummer+1):((i*nb_days_by_latesummer+30*nb_days_by_latesummer))]
  Ind_SlidingWindow30_latespringMJ_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-
    Ind_latespringMJ_1850_2100[((i)*nb_days_by_latespringMJ+1):((i*nb_days_by_latespringMJ+30*nb_days_by_latespringMJ))]
  Ind_SlidingWindow30_latespringAMJ_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-
    Ind_latespringAMJ_1850_2100[((i)*nb_days_by_latespringAMJ+1):((i*nb_days_by_latespringAMJ+30*nb_days_by_latespringAMJ))]
  Ind_SlidingWindow30_latespringJ_1850_2100[[label_sliding_window[length(label_sliding_window)]]]<-
    Ind_latespringJ_1850_2100[((i)*nb_days_by_latespringJ+1):((i*nb_days_by_latespringJ+30*nb_days_by_latespringJ))]

}


Label_SlidingWindow30_1850_2100=label_sliding_window

setwd(dir="/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data")
save(Ind_spring_1850_2100,Ind_summer_1850_2100,Ind_automn_1850_2100,Ind_winter_1850_2100,
     Ind_latespringMJ_1850_2100,
     Ind_latespringJ_1850_2100,
     Ind_latespringAMJ_1850_2100,
     Ind_latesummer_1850_2100, Ind_months_1850_2100,
     #### CVchrono
     Label_SlidingWindow30_1850_2100,
     Ind_SlidingWindow30_winter_1850_2100,
     Ind_SlidingWindow30_automn_1850_2100,
     Ind_SlidingWindow30_summer_1850_2100,
     Ind_SlidingWindow30_spring_1850_2100,
     Ind_SlidingWindow30_latesummer_1850_2100,
     Ind_SlidingWindow30_latespringMJ_1850_2100,
     Ind_SlidingWindow30_latespringAMJ_1850_2100,
     Ind_SlidingWindow30_latespringJ_1850_2100,
     file="Compound_Event_Temporal_indices_SlidingWindow30_1850_2100.RData")

##################################################################################################################################################
##################################################################################################################################################
####### Periode 1950_2020 a decouper

rm(list=ls())
set.seed(42)
library(timeDate)

#Dropping 29 of February for each year in the reference
DATES_regular_1950_2020 = atoms(timeSequence(from="1950-01-01",to="2020-12-31",by='day'))
head(DATES_regular_1950_2020)

Ind_29Feb_1950_2020 = which((DATES_regular_1950_2020[,2]==02) &
                                    (DATES_regular_1950_2020[,3]==29))
DATES_365_1950_2020 = DATES_regular_1950_2020[-Ind_29Feb_1950_2020,] #remove rows with 29 Feb.

#For Seasons (meteorological)
Ind_spring_1950_2020 =  which(DATES_365_1950_2020[,2]==03 |
                                  DATES_365_1950_2020[,2]==04 |
                                  DATES_365_1950_2020[,2]==05)
Ind_summer_1950_2020 = which(DATES_365_1950_2020[,2]==06|
                                 DATES_365_1950_2020[,2]==07|
                                 DATES_365_1950_2020[,2]==08)
Ind_automn_1950_2020 =  which(DATES_365_1950_2020[,2]==09|
                                  DATES_365_1950_2020[,2]==10 |
                                  DATES_365_1950_2020[,2]==11)
Ind_winter_1950_2020 =  which(DATES_365_1950_2020[,2]==12|
                                  DATES_365_1950_2020[,2]==01 |
                                  DATES_365_1950_2020[,2]==02)

#to del
Ind_latespringMJ_1950_2020 =  which(DATES_365_1950_2020[,2]==05 |
                                  DATES_365_1950_2020[,2]==06)
Ind_latespringAMJ_1950_2020 =  which(DATES_365_1950_2020[,2]==04 |
                                      DATES_365_1950_2020[,2]==05 |
                                  DATES_365_1950_2020[,2]==06)
Ind_latespringJ_1950_2020 =  which(DATES_365_1950_2020[,2]==06)
Ind_latesummer_1950_2020 =  which(DATES_365_1950_2020[,2]==07 |
                                  DATES_365_1950_2020[,2]==08)

# For months
Ind_months_1950_2020=list()
for(i in 1:12){
        Ind_months_1950_2020[[i]]=which(DATES_365_1950_2020[,2]==i)
}

### For sliding windows of 50 years
Ind_SlidingWindow30_winter_1950_2020=list()
Ind_SlidingWindow30_spring_1950_2020=list()
Ind_SlidingWindow30_summer_1950_2020=list()
Ind_SlidingWindow30_automn_1950_2020=list()
Ind_SlidingWindow30_latespringMJ_1950_2020=list()
Ind_SlidingWindow30_latespringAMJ_1950_2020=list()
Ind_SlidingWindow30_latespringJ_1950_2020=list()
Ind_SlidingWindow30_latesummer_1950_2020=list()



start=1950
end=1979

nb_SlidingWindow30=41
label_sliding_window=c()

nb_days_by_winter=length(Ind_winter_1950_2020)/71
nb_days_by_summer=length(Ind_summer_1950_2020)/71
nb_days_by_automn=length(Ind_automn_1950_2020)/71
nb_days_by_spring=length(Ind_spring_1950_2020)/71

nb_days_by_latespringMJ=length(Ind_latespringMJ_1950_2020)/71
nb_days_by_latespringAMJ=length(Ind_latespringAMJ_1950_2020)/71
nb_days_by_latespringJ=length(Ind_latespringJ_1950_2020)/71
nb_days_by_latesummer=length(Ind_latesummer_1950_2020)/71



### Between 1950 and 2020, 42 sliding window of size 50
for(i in 0: nb_SlidingWindow30){
  label_sliding_window=c(label_sliding_window,paste0(start+i,"_",end+i))
  print(label_sliding_window)
  Ind_SlidingWindow30_winter_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_winter_1950_2020[((i)*nb_days_by_winter+1):((i*nb_days_by_winter+30*nb_days_by_winter))]
  Ind_SlidingWindow30_automn_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_automn_1950_2020[((i)*nb_days_by_automn+1):((i*nb_days_by_automn+30*nb_days_by_automn))]
  Ind_SlidingWindow30_summer_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_summer_1950_2020[((i)*nb_days_by_summer+1):((i*nb_days_by_summer+30*nb_days_by_summer))]
  Ind_SlidingWindow30_spring_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_spring_1950_2020[((i)*nb_days_by_spring+1):((i*nb_days_by_spring+30*nb_days_by_spring))]

  Ind_SlidingWindow30_latesummer_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_latesummer_1950_2020[((i)*nb_days_by_latesummer+1):((i*nb_days_by_latesummer+30*nb_days_by_latesum
mer))]
  Ind_SlidingWindow30_latespringMJ_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_latespringMJ_1950_2020[((i)*nb_days_by_latespringMJ+1):((i*nb_days_by_latespringMJ+30*nb_days_by
_latespringMJ))]
  Ind_SlidingWindow30_latespringAMJ_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_latespringAMJ_1950_2020[((i)*nb_days_by_latespringAMJ+1):((i*nb_days_by_latespringAMJ+30*nb_day
s_by_latespringAMJ))]
  Ind_SlidingWindow30_latespringJ_1950_2020[[label_sliding_window[length(label_sliding_window)]]]<-Ind_latespringJ_1950_2020[((i)*nb_days_by_latespringJ+1):((i*nb_days_by_latespringJ+30*nb_days_by_lat
espringJ))]
}

Label_SlidingWindow30_1950_2020=label_sliding_window

setwd(dir="/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data")
save(Ind_spring_1950_2020,Ind_summer_1950_2020,Ind_automn_1950_2020,Ind_winter_1950_2020,
     Ind_latespringMJ_1950_2020,
     Ind_latespringAMJ_1950_2020,
     Ind_latespringJ_1950_2020,
     Ind_latesummer_1950_2020,
     Ind_months_1950_2020,
     #### CVchrono
     Label_SlidingWindow30_1950_2020,
     Ind_SlidingWindow30_winter_1950_2020,
     Ind_SlidingWindow30_automn_1950_2020,
     Ind_SlidingWindow30_summer_1950_2020,
     Ind_SlidingWindow30_spring_1950_2020,
     Ind_SlidingWindow30_latespringMJ_1950_2020,
     Ind_SlidingWindow30_latespringAMJ_1950_2020,
     Ind_SlidingWindow30_latespringJ_1950_2020,
     Ind_SlidingWindow30_latesummer_1950_2020,
     file="Compound_Event_Temporal_indices_SlidingWindow30_1950_2020.RData")

