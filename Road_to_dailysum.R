rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Code/function_TCE.R")
library(copula)
library(VineCopula)
library(viridis)
library(RColorBrewer)
library(scales)

load("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/Compound_Event_Temporal_indices_SlidingWindow30_1850_2100.RData")
load("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/Compound_Event_Temporal_indices_SlidingWindow30_1950_2020.RData")

name_var1="sfcWindmax"
name_var2="pr"
season="winter"
Mod="CNRMCM6" #IPSL ERA5
Region="Bretagne"

if(Mod %in% c("CNRMCM6", "IPSL")){
  period="1850_2100"
  nb_years_=251
  Label_SlidingWindow30= get(paste0("Label_SlidingWindow30_", period))
}
if(Mod=="ERA5"){
  period="1950_2020"
  nb_years_=71
  Label_SlidingWindow30= get(paste0("Label_SlidingWindow30_", period))
}
Ind_season=get(paste0("Ind_",season,"_",period))

load(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/RefMod_Data/",
            Mod, "/", name_var1, "_day_", Mod, "_", period, "_", Region,
            ".RData"))
load(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/RefMod_Data/",
            Mod, "/", name_var2, "_day_", Mod, "_", period, "_", Region,
            ".RData"))

var1=get(paste0(name_var1,"_day_", Mod, "_",period, "_", Region))[,,Ind_season]
var2=get(paste0(name_var2,"_day_", Mod, "_",period, "_", Region))[,,Ind_season]

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Plots/",
             Region, "/", name_var1, "_", name_var2,"/", Mod))


#### Daily mean of Windmax and daily sum of pr over the region
index1=compute_daily_spat_mean(var1)
index2=compute_daily_spat_sum(var2)

nb_index_by_year_=length(Ind_season)/nb_years_

reord_in_list<-function(index, label_sw, nb_index_by_year,
                        length_sw=30){
  res_list=list()
  k=0
  for(l in label_sw){
    res_list[[l]]=index[(k*nb_index_by_year+1):(k*nb_index_by_year+nb_index_by_year*length_sw)]
    k=k+1
  }
  return(res_list)
}


tmp_list30_index1<-reord_in_list(index1, Label_SlidingWindow30,
                                 nb_index_by_year_)
tmp_list30_index2<-reord_in_list(index2, Label_SlidingWindow30,  nb_index_by_year_)
uncond_list30_index1<- tmp_list30_index1
uncond_list30_index2<- tmp_list30_index2

quant90_index1=quantile(tmp_list30_index1[[1]],probs=0.9, na.rm=T)
quant90_index2=quantile(tmp_list30_index2[[1]],probs=0.9, na.rm=T)

### Plot
par(mfrow=c(2,5))
for(l in Label_SlidingWindow30[trunc(seq(1,length(Label_SlidingWindow30),
                                                 length.out=10))]){
  print(l)
  tmp_cor=cor(tmp_list30_index1[[l]],
                        tmp_list30_index2[[l]], method="spearman")
  plot(tmp_list30_index1[[l]],
                        tmp_list30_index2[[l]],
                        ylim=c(min(unlist(tmp_list30_index2),
                                   na.rm=T),max(unlist(tmp_list30_index2),
                                   na.rm=T)),
                        xlim=c(min(unlist(tmp_list30_index1),
                                          na.rm=T),max(unlist(tmp_list30_index1),
                                   na.rm=T)),
                        main=paste0(l,", cor: ",
                                                  round(tmp_cor,2)),
       ylab="Daily sum PR",
       xlab="Spat. Mean of max Wind")
  abline(h=quant90_index2, col="red")
  abline(v=quant90_index1, col="red")
}


### Conditionning to 90 percentiles (AND)
list30_index1=list()
list30_index2=list()

res_cormax=res_cormaxpval=c()
k=0
par(mfrow=c(2,5))
for(l in Label_SlidingWindow30){
  print(l)
  coord_phy=which(tmp_list30_index1[[l]]>=quant90_index1
                  & tmp_list30_index2[[l]]>=quant90_index2)
  list30_index1[[l]]=tmp_list30_index1[[l]][coord_phy]
  list30_index2[[l]]=tmp_list30_index2[[l]][coord_phy]
  print(length(coord_phy))

}


par(mfrow=c(2,5))
for(l in Label_SlidingWindow30[trunc(seq(1,length(Label_SlidingWindow30),
                                                 length.out=10))]){
  print(l)
  tmp_index1= list30_index1[[l]]
  tmp_index2= list30_index2[[l]]
  plot(tmp_index1, tmp_index2,
                        ylim=c(min(unlist(list30_index2),
                                   na.rm=T),max(unlist(list30_index2),
                                   na.rm=T)),
                        xlim=c(min(unlist(list30_index1),
                                          na.rm=T),max(unlist(list30_index1),
                                   na.rm=T)),
                        main=paste0(l,", cor: ",
                                                  round(cor(tmp_index1,
                                                            tmp_index2,
                                                            method="spearman"),2)),
        ylab="Daily sum PR",
       xlab="Spat. Mean of max Wind")
}




#### Begin A famille fixe
probs_to_eval=seq(0.05, 0.95,0.05)
probs_to_eval_index1= probs_to_eval_index2=probs_to_eval

Init_index1<-list30_index1[[1]]
Init_index2<-list30_index2[[1]]

###What is the quantile associated with probs_to_eval in 1850_1899?
quant_Init_index1=Fm1_GPD(Init_index1, quant90_index1, probs_to_eval_index1)
quant_Init_index2=Fm1_GPD(Init_index2, quant90_index2, probs_to_eval_index2)

probs_list30_index1=matrix(NaN, ncol=length(probs_to_eval_index1), nrow=
                              length(Label_SlidingWindow30)) ###matrix SlidingWindow x probs
probs_list30_index2=matrix(NaN, ncol=length(probs_to_eval_index2), nrow=
                                   length(Label_SlidingWindow30)) ###matrix SlidingWindow x probs


res_Goftest_best_Copula=c()
res_theta_best_Copula_68=list()
res_theta_best_Copula_95=list()

res_name_best_Copula=c()
array_prob_margdep_68_best_Copula <- array_prob_marg_68_best_Copula <-
  array_prob_dep_68_best_Copula <-array(NaN,dim=c(length(probs_to_eval_index1),
                                               length(probs_to_eval_index2),
                  length(Label_SlidingWindow30),3)) #drought x proba HW x Slid x CI
array_prob_margdep_95_best_Copula <- array_prob_marg_95_best_Copula <-
  array_prob_dep_95_best_Copula <-array(NaN,dim=c(length(probs_to_eval_index1),
                                               length(probs_to_eval_index2),
                  length(Label_SlidingWindow30),3)) #drought x proba HW x Slid x CI


family_to_try=c("Frank", "Joe", "Clayton", "Gumbel")
for(f in family_to_try){
  assign(paste0("res_Goftest_",f),c())
  assign(paste0("res_theta_",f, "_68"),list())
  assign(paste0("res_theta_",f, "_95"),list())
  assign(paste0("array_prob_margdep_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
  assign(paste0("array_prob_marg_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
  assign(paste0("array_prob_dep_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
  assign(paste0("array_prob_margdep_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
  assign(paste0("array_prob_marg_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
  assign(paste0("array_prob_dep_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindow30),3))) #drought x proba HW x Slid x CI
}

k=0
for(l in Label_SlidingWindow30){
  print(l)
  k=k+1
  ### faire la methodo super genial
  tmp_list30_index1<-list30_index1[[l]]
  tmp_list30_index2<-list30_index2[[l]]
  ### Evaluate the change of probs (seuil univarie pour marginal prop.)
  #### to which proba the quantile associated with probs 0.5 in 1850-1899 correspond to in
  #Sliding_Window?
  probs_list30_index1[k,]=F_GPD(tmp_list30_index1, quant90_index1, quant_Init_index1)
  probs_list30_index2[k,]=F_GPD(tmp_list30_index2, quant90_index2, quant_Init_index2)
  ### Evaluate the change of probs (seuil bivarie pour marginal+dep prop.)
  ### Fit Copula 
  pobs_list30_index1=pobs(tmp_list30_index1)#F_GPD(tmp_list30_index1, quant90_index1, tmp_list30_index1)
  pobs_list30_index2=pobs(tmp_list30_index2)#F_GPD(tmp_list30_index2, quant90_index2, tmp_list30_index2)

  U_pobs=matrix(NaN,ncol=2,nrow=length(pobs_list30_index1))
  U_pobs[,1]=pobs_list30_index1
  U_pobs[,2]=pobs_list30_index2
  ####1. Selection among Archimedean family (without rotations, i.e. survival)
  cop_best_Copula=selectedCopula(pobs_list30_index1, pobs_list30_index2,
                     pobs_to_calc=FALSE)
  res_name_best_Copula<-c(res_name_best_Copula, cop_best_Copula$familyname)

  ### Save cop of the first Sliding Window
  if(k==1){
    cop_Init_best_Copula<-cop_best_Copula
    estim_ci_theta_Init_best_Copula_68=determine_ci_theta(cop_Init_best_Copula, U_pobs, ci_level=0.68)
    estim_ci_theta_Init_best_Copula_95=determine_ci_theta(cop_Init_best_Copula, U_pobs, ci_level=0.95)
  }

  if(cop_best_Copula$family==3){name_family="Clayton"}
  if(cop_best_Copula$family==4){name_family="Gumbel"}
  if(cop_best_Copula$family==5){name_family="Frank"}
  if(cop_best_Copula$family==6){name_family="Joe"}
  ## 2. Goodness of fit: is the fit good?
  if(cop_best_Copula$family==0){
    print("Test Indep.")
    test=testIndepCopula(pobs_list30_index1, pobs_list30_index2,
    pobs_to_calc=FALSE)
  }else{
    test=testCopula(pobs_list30_index1, pobs_list30_index2,
                    cop_best_Copula$family, pobs_to_calc=FALSE)
  }
  #H0: the fit is good. p<0.05 we reject
  res_Goftest_best_Copula=c(res_Goftest_best_Copula,test$p.value)
  print(paste0(name_family, ", par: ", round(cop_best_Copula$par,2),", GoF test: ",
               round(test$p.value,3)))
  ## 3. Estimate theta_hat and the CI of the theta_hat
  if(name_family !="Independence"){
    print(paste0("Selected family: ", name_family))
    estim_ci_theta_list30_68=determine_ci_theta(cop_best_Copula, U_pobs, ci_level=0.68)
    estim_ci_theta_list30_95=determine_ci_theta(cop_best_Copula, U_pobs, ci_level=0.95)

    res_theta_best_Copula_68[[l]]<-estim_ci_theta_list30_68
    res_theta_best_Copula_95[[l]]<-estim_ci_theta_list30_95
    ##4. Margdep
    #Estimate the proba and the CI of the desired proba
    array_prob_margdep_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_best_Copula,
                                                                           estim_ci_theta_list30_68, probs_list30_index1[k,], probs_list30_index2[k,])
    array_prob_margdep_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_best_Copula,
                                                                           estim_ci_theta_list30_95, probs_list30_index1[k,], probs_list30_index2[k,])
    ## 5. Marginal
    array_prob_marg_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_68, probs_list30_index1[k,], probs_list30_index2[k,])
    array_prob_marg_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_95, probs_list30_index1[k,], probs_list30_index2[k,])
    ## 6. Dep
    array_prob_dep_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_best_Copula,
                                                                    estim_ci_theta_list30_68,
                                                                    probs_to_eval_index1,
                                                                    probs_to_eval_index2)
    array_prob_dep_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array(cop_best_Copula,
                                                                    estim_ci_theta_list30_95,
                                                                    probs_to_eval_index1,
                                                                    probs_to_eval_index2)

  }

  #### 7. Re-do the methodo for each f family fixed
  for(f in family_to_try){
   print(f)
   if(f=="Clayton"){f_number=3}
   if(f=="Gumbel"){f_number=4}
   if(f=="Frank"){f_number=5}
   if(f=="Joe"){f_number=6}
   ### 7.2 Goofness of fit test
   assign(paste0("res_Goftest_",f), c(get((paste0("res_Goftest_",f))),testCopula(pobs_list30_index1, pobs_list30_index2,
                  f_number, pobs_to_calc=FALSE)$p.value))
   ### 7.3 Estim. theta at fixed family
   est_cop_f=BiCopEst(pobs_list30_index1, pobs_list30_index2, f_number)
  ### Save cop of the first list30 Window
   if(k==1){
    assign(paste0("cop_Init_",f), est_cop_f)
    assign(paste0("estim_ci_theta_Init_",f, "_68"),determine_ci_theta(get(paste0("cop_Init_",f)),
                                                               U_pobs,
                                                               ci_level=0.68))
    assign(paste0("estim_ci_theta_Init_",f, "_95"),determine_ci_theta(get(paste0("cop_Init_",f)),
                                                               U_pobs,
                                                               ci_level=0.95))
   }
   ### Define cop object
   tmp_cop_f=BiCop(f_number, est_cop_f$par, par2 = 0)
   print(tmp_cop_f)
   estim_ci_theta_f_68=determine_ci_theta(tmp_cop_f, U_pobs, ci_level=0.68)
   estim_ci_theta_f_95=determine_ci_theta(tmp_cop_f, U_pobs, ci_level=0.95)
   eval(parse(text=paste0("res_theta_",f,"_68[['", l, "']] <- estim_ci_theta_f_68")))
   eval(parse(text=paste0("res_theta_",f,"_95[['", l, "']] <- estim_ci_theta_f_95")))
   ##4. Margdep
   #Estimate the proba and the CI of the desired proba
   tmp_prob_margdep_68<-determine_ci_prob_bivar_array(tmp_cop_f,
                                                      estim_ci_theta_f_68, probs_list30_index1[k,], probs_list30_index2[k,])
   tmp_prob_margdep_95<-determine_ci_prob_bivar_array(tmp_cop_f,
                                                      estim_ci_theta_f_95, probs_list30_index1[k,], probs_list30_index2[k,])
   eval(parse(text=paste0("array_prob_margdep_68_", f,"[,,k,]<-
                          tmp_prob_margdep_68")))
   eval(parse(text=paste0("array_prob_margdep_95_", f,"[,,k,]<-
                          tmp_prob_margdep_95")))

   tmp_prob_marg_68<-determine_ci_prob_bivar_array(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_68")), probs_list30_index1[k,], probs_list30_index2[k,])
   tmp_prob_marg_95<-determine_ci_prob_bivar_array(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_95")), probs_list30_index1[k,], probs_list30_index2[k,])
   eval(parse(text=paste0("array_prob_marg_68_", f,"[,,k,]<-
                          tmp_prob_marg_68")))

   eval(parse(text=paste0("array_prob_marg_95_", f,"[,,k,]<-
                          tmp_prob_marg_95")))

   tmp_prob_dep_68<-determine_ci_prob_bivar_array(tmp_cop_f, estim_ci_theta_f_68,
                                                probs_to_eval_index1,
                                                probs_to_eval_index2)
   tmp_prob_dep_95<-determine_ci_prob_bivar_array(tmp_cop_f,
                                                  estim_ci_theta_f_95,
                                                probs_to_eval_index1,
                                                probs_to_eval_index2)
   eval(parse(text=paste0("array_prob_dep_68_", f,"[,,k,]<-
                          tmp_prob_dep_68")))
   eval(parse(text=paste0("array_prob_dep_95_", f,"[,,k,]<-
                          tmp_prob_dep_95")))
  }
}
name_array=c()
for(f in c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton")){
  for(v in c("margdep", "marg", "dep")){
    for(ci in c("68", "95")){
    name_array=c(name_array,paste0("array_prob_",v,"_",ci,"_", f))
    }
  }
}

for(f in c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton")){
  name_array=c(name_array, paste0("res_Goftest_",f),
               paste0("res_theta_",f,"_68"), paste0("res_theta_",f,"_68"))
}


object_to_save=c(name_array,
                 "res_name_best_Copula",
                 "probs_to_eval_index1",
                 "probs_to_eval_index2",
                 "probs_list30_index1",
                 "probs_list30_index2",
                 "quant_Init_index1",
                 "quant_Init_index2",
                 "list30_index1",
                 "list30_index2",
                 "Label_SlidingWindow30",
                 "uncond_list30_index1",
                 "uncond_list30_index2")

name_to_save=c()
for(obj in object_to_save){
  name_to_save=c(name_to_save, paste0(obj, "_", Mod))
  assign(paste0(obj, "_", Mod), get(paste0(obj)))
}


save_wd=getwd()
setwd(paste0(save_wd,"/Save_Data"))

save(list=name_to_save,
     file=paste0("array_prob_list30_Family_fixed_", name_var1, "_", name_var2, "_", Mod,
                 "_", season, "_", period, "_", Region, ".RData"))


setwd(paste0(save_wd))

print(ok)
