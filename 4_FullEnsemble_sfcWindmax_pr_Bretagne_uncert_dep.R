#!/usr/bin/env Rscript
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)
ModRef=args[1]
Ref_SlidingWindow=args[2]
#ModRef="CNRMCM6"
#Ref_SlidingWindow="1871_1900"
print(Ref_SlidingWindow)
library(CDFt)
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_1dQQ.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
library(copula)
library(VineCopula)
library(viridis)
library(RColorBrewer)
library(scales)
load("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Compound_Event_Temporal_indices_SlidingWindows_1850_2100.RData")

name_var1="sfcWindmax"
name_var2="pr"
season="winter"
Region="Bretagne"
is_TardiveFrost=FALSE
Length_SlidingWindow=30
if(is_TardiveFrost==FALSE){zone_="topright"}
period="1850_2100"
nb_years_=251
Label_SlidingWindowX= get(paste0("Label_SlidingWindow", Length_SlidingWindow,"_", period))
coord_Ref_SlidingWindow=which(Label_SlidingWindowX==Ref_SlidingWindow)
nb_listX=length(Label_SlidingWindowX)
Ind_season=get(paste0("Ind_",season,"_",period))

### Load Model indices 
list_Mod=c("CNRMCM6", "CMCCESM2", "CanESM5", "CNRMCM6HR", "FGOALSG3", "IPSL",
           "ECEARTH3", "INMCM48", "INMCM50", "GFDLCM4", "MIROC6", "MPIESM1LR",
           "MRIESM2")
# c("CNRMCM6", "IPSL", "MIROC6", "INMCM48", "CanESM5",
            #  "GFDLCM4", "MPIESM1LR", "MRIESM2", "CNRMCM6HR", "CMCCESM2",
            #  "ECEARTH3", "FGOALSG3", "INMCM50")
for(Mod in list_Mod){
  print(Mod)
  load(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/",
            Region, "/", name_var1, "_", name_var2, "/SlidingWindow30_Ref",
            Ref_SlidingWindow, "/Res_array_prob_list30_Ref", Ref_SlidingWindow,
            "_", name_var1, "_", name_var2, "_", Mod, "_", season, "_",
            period, "_", Region, ".RData"))
}

### Define ModRef
tmp_ess1=tmp_ess2=c()
for(Mod in list_Mod){
  print(Mod)
  print(get(paste0("quant90_index1_", Mod)))
  tmp_ess1=c(tmp_ess1,get(paste0("quant90_index1_", Mod)))
}
for(Mod in list_Mod){
  print(Mod)
  print(get(paste0("quant90_index2_", Mod)))
  tmp_ess2=c(tmp_ess2,get(paste0("quant90_index2_", Mod)))
}

plot(tmp_ess2~tmp_ess1, col="red", ylab="q90 Rain acc.", xlab="q90 Wind speed")
text(tmp_ess2 ~tmp_ess1, labels=list_Mod, cex=0.75, font=1)

ModRef_vecX_index1=get(paste0("listX_index1_", ModRef))[[Ref_SlidingWindow]]
ModRef_vecX_index2=get(paste0("listX_index2_", ModRef))[[Ref_SlidingWindow]]

list_Mod_without_ModRef=list_Mod[which(ModRef != list_Mod)]
 
### Load indices for the other Mod and do CDFt
npas_=floor(length(ModRef_vecX_index1)*2/3)

### For seuillage CDFt
min_index1_ModRef=min(unlist(get(paste0("listX_index1_", ModRef))))
min_index2_ModRef=min(unlist(get(paste0("listX_index2_", ModRef))))
quant90_index1_ModRef=get(paste0("quant90_index1_", ModRef))
quant90_index2_ModRef=get(paste0("quant90_index2_", ModRef))

for(Mod in list_Mod_without_ModRef){
  print(Mod)
 
  tmp_Mod_listX_index1=get(paste0("listX_index1_", Mod))
  tmp_Mod_listX_index2=get(paste0("listX_index2_", Mod))

  tmp_CDFt_listX_index1=list()
  tmp_CDFt_listX_index2=list()

  #### Correction
  for(l in Label_SlidingWindowX){
    print(l)
#    if(l==Ref_SlidingWindow){#QQ
##      tmp_CDFt_listX_index1[[l]]=QQb_new(ModRef_vecX_index1,
##                                         tmp_Mod_listX_index1[[coord_Ref_SlidingWindow]],
##                                         tmp_Mod_listX_index1[[coord_Ref_SlidingWindow]])$Mch[,1] 
##      tmp_CDFt_listX_index2[[l]]=QQb_new(ModRef_vecX_index2,
##                                         tmp_Mod_listX_index2[[coord_Ref_SlidingWindow]],
##                                         tmp_Mod_listX_index2[[coord_Ref_SlidingWindow]])$Mch[,1] 
##      ### new version QQ parametric
#      tmp_q0_index1=get(paste0("quant90_index1_", Mod))
#      tmp_q0_index2=get(paste0("quant90_index2_", Mod))
#      ModRef_q0_index1=get(paste0("quant90_index1_", ModRef))
#      ModRef_q0_index2=get(paste0("quant90_index2_", ModRef))
#
#      tmp_CDFt_listX_index1[[l]]=QQb_new_GPD(ModRef_vecX_index1,
#                                         tmp_Mod_listX_index1[[coord_Ref_SlidingWindow]],
#                                         tmp_Mod_listX_index1[[coord_Ref_SlidingWindow]],
#                                         ModRef_q0_index1, tmp_q0_index1)$Mch[,1] 
#      tmp_CDFt_listX_index2[[l]]=QQb_new_GPD(ModRef_vecX_index2,
#                                         tmp_Mod_listX_index2[[coord_Ref_SlidingWindow]],
#                                         tmp_Mod_listX_index2[[coord_Ref_SlidingWindow]],
#                                         ModRef_q0_index2, tmp_q0_index2)$Mch[,1] 
#
#
#
#    }else{
      tmp_CDFt_listX_index1[[l]]=CDFt(ModRef_vecX_index1,
                                         tmp_Mod_listX_index1[[coord_Ref_SlidingWindow]],
                                         tmp_Mod_listX_index1[[l]], npas=npas_)$DS 
      tmp_CDFt_listX_index2[[l]]=CDFt(ModRef_vecX_index2,
                                         tmp_Mod_listX_index2[[coord_Ref_SlidingWindow]],
                                         tmp_Mod_listX_index2[[l]], npas=npas_)$DS 
#    }
#    ### Seuillage
#    #Idea: tous les points suite a CDFt sous quant90_index1, je fais un runif
#    #entre quant90_index1 et min(listX_index1_ModRef)
#    Windex1 = which(tmp_CDFt_listX_index1[[l]]< quant90_index1_ModRef)
#   #if(length(Windex1)>0){print(l)}
#    tmp_CDFt_listX_index1[[l]][Windex1] <- runif(length(Windex1),quant90_index1_ModRef, min_index1_ModRef)
#    Windex2 = which(tmp_CDFt_listX_index2[[l]]< quant90_index2_ModRef)
#    tmp_CDFt_listX_index2[[l]][Windex2] <- runif(length(Windex2),quant90_index2_ModRef, min_index2_ModRef)
    #Idea: tous les points suite a CDFt sous quant90_index1, on les vire
    Windex1 = which(tmp_CDFt_listX_index1[[l]]< quant90_index1_ModRef)
    Windex2 = which(tmp_CDFt_listX_index2[[l]]< quant90_index2_ModRef)
    if(length(c(Windex1, Windex2))>0){
      print(l)
      tmp_CDFt_listX_index1[[l]] <- tmp_CDFt_listX_index1[[l]][-c(Windex1,
                                                                Windex2)]
      tmp_CDFt_listX_index2[[l]] <- tmp_CDFt_listX_index2[[l]][-c(Windex1, Windex2)]
    }
  }
  assign(paste0("listX_index1_CDFt_", Mod, "_ModRef", ModRef), tmp_CDFt_listX_index1)
  assign(paste0("listX_index2_CDFt_", Mod, "_ModRef", ModRef), tmp_CDFt_listX_index2)
}  
assign(paste0("listX_index1_ModRef", ModRef), get(paste0("listX_index1_", ModRef)))
assign(paste0("listX_index2_ModRef", ModRef), get(paste0("listX_index2_", ModRef)))

print("ok")
### Pooling of the CDF-t correction
assign(paste0("listX_index1_CDFt_pooled_ModRef", ModRef), list())
assign(paste0("listX_index2_CDFt_pooled_ModRef", ModRef), list())

for(l in Label_SlidingWindowX){
  tmp_pooled_index1=tmp_pooled_index2=c()
  ### Mod_without_ModRef
  for(Mod in list_Mod_without_ModRef){
    tmp_pooled_index1=c(tmp_pooled_index1, unname(unlist(get(paste0("listX_index1_CDFt_", Mod,
                                               "_ModRef", ModRef))[[l]])))
    tmp_pooled_index2=c(tmp_pooled_index2, unname(unlist(get(paste0("listX_index2_CDFt_", Mod,
                                               "_ModRef", ModRef))[[l]])))
  }
  ### Put ModRef
  tmp_pooled_index1=c(tmp_pooled_index1, unname(unlist(get(paste0("listX_index1_", ModRef))[[l]])))
  tmp_pooled_index2=c(tmp_pooled_index2, unname(unlist(get(paste0("listX_index2_", ModRef))[[l]])))
  
  eval(parse(text=paste0("listX_index1_CDFt_pooled_ModRef", ModRef,"[[l]] <-
                         tmp_pooled_index1")))
  eval(parse(text=paste0("listX_index2_CDFt_pooled_ModRef", ModRef,"[[l]] <-
                         tmp_pooled_index2")))
} 

print("okok")
par(mfrow=c(2,5))
for(l in Label_SlidingWindowX[trunc(seq(22,31,
                                                 length.out=10))]){
  tmp_index1= get(paste0("listX_index1_CDFt_pooled_ModRef", ModRef))[[l]]
  tmp_index2= get(paste0("listX_index2_CDFt_pooled_ModRef", ModRef))[[l]]
  plot(tmp_index1, tmp_index2,
                        xlim=c(min(unlist(get(paste0("listX_index1_CDFt_pooled_ModRef",
                                                     ModRef))),
                                   na.rm=T),max(unlist(get(paste0("listX_index1_CDFt_pooled_ModRef",
                                                     ModRef))),
                                   na.rm=T)),
                        ylim=c(min(unlist(get(paste0("listX_index2_CDFt_pooled_ModRef",
                                                     ModRef))),
                                   na.rm=T),max(unlist(get(paste0("listX_index2_CDFt_pooled_ModRef",
                                                     ModRef))),
                                   na.rm=T)),
                        main=paste0(l,", cor: ", round(cor(tmp_index1, tmp_index2,
                                                            method="spearman"),2),
                                    ", pval: ", round(cor.test(tmp_index1, tmp_index2,
                                                            method="spearman")$p.value,2)),
        ylab="Daily sum PR",
       xlab="Spat. Mean of max Wind")
  
}





################## Copula proba with pooled_CDFt ################
print("copula")
listX_index1<-get(paste0("listX_index1_CDFt_pooled_ModRef", ModRef))
listX_index2<-get(paste0("listX_index2_CDFt_pooled_ModRef", ModRef))

quant90_index1<- quant90_index1_ModRef
quant90_index2<- quant90_index2_ModRef



### Exact copy from 1_4 
#### Results
probs_to_eval=seq(0.05, 0.95,0.05)
probs_to_eval_index1= probs_to_eval_index2=probs_to_eval

probs_listX_index1=matrix(NaN, ncol=length(probs_to_eval_index1), nrow=
                              length(Label_SlidingWindowX)) ###matrix SlidingWindow x probs
probs_listX_index2=matrix(NaN, ncol=length(probs_to_eval_index2), nrow=
                                   length(Label_SlidingWindowX)) ###matrix SlidingWindow x probs
res_name_best_Copula=rep(NaN, nb_listX)

family_to_try=c("Frank", "Joe", "Clayton", "Gumbel")
for(f in c(family_to_try, "best_Copula")){
  assign(paste0("res_Goftest_",f),rep(NaN, nb_listX))
  assign(paste0("res_theta_",f, "_68"),list())
  assign(paste0("res_theta_",f, "_95"),list())
  assign(paste0("array_prob_margdep_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3))) 
  assign(paste0("array_prob_marg_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3))) 
  assign(paste0("array_prob_dep_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3))) 
  assign(paste0("array_prob_margdep_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3)))
  assign(paste0("array_prob_marg_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3)))
  assign(paste0("array_prob_dep_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3)))
  ### Infl_marg1: comme marg mais en fixant index2
  assign(paste0("array_prob_Inflmarg1_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3))) 
  assign(paste0("array_prob_Inflmarg1_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3)))
  ### Infl_marg2: comme marg mais en fixant index1
  assign(paste0("array_prob_Inflmarg2_68_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3))) 
  assign(paste0("array_prob_Inflmarg2_95_",f), array(NaN,dim=c(length(probs_to_eval_index1), length(probs_to_eval_index1),
                  length(Label_SlidingWindowX),3)))
} 


#### Compute Copulas
print("###################### Init ####################")
print("###################### Best Copula ####################")

Init_index1<-listX_index1[[coord_Ref_SlidingWindow]]
Init_index2<-listX_index2[[coord_Ref_SlidingWindow]]
###What is the quantile associated with probs_to_eval in Ref_SlidingWindow?
quant_Init_index1=Fm1_GPD(Init_index1, quant90_index1, probs_to_eval_index1)
quant_Init_index2=Fm1_GPD(Init_index2, quant90_index2, probs_to_eval_index2)

probs_listX_index1[coord_Ref_SlidingWindow,]=F_GPD(Init_index1, quant90_index1,
                             quant_Init_index1)
probs_listX_index2[coord_Ref_SlidingWindow,]=F_GPD(Init_index2, quant90_index2,
                             quant_Init_index2)
pobs_Init_index1=pobs(Init_index1)#F_GPD(tmp_listX_index1,quant90_index1, tmp_listX_index1)
pobs_Init_index2=pobs(Init_index2)#F_GPD(tmp_listX_index2,quant90_index2, tmp_listX_index2)
U_pobs_Init=matrix(NaN,ncol=2,nrow=length(pobs_Init_index1))
U_pobs_Init[,1]=pobs_Init_index1
U_pobs_Init[,2]=pobs_Init_index2

#### Selection among Archimedean family (without rotations, i.e. survival)
cop_Init_best_Copula=selectedCopula3(U_pobs_Init, pobs_to_calc=FALSE)
estim_ci_theta_Init_best_Copula_68=determine_ci_theta2(cop_Init_best_Copula,
                                                      U_pobs_Init, ci_level=0.68)
estim_ci_theta_Init_best_Copula_95=determine_ci_theta2(cop_Init_best_Copula,
                                                      U_pobs_Init, ci_level=0.95)
print(paste0("Best Copula: ", cop_Init_best_Copula$family, ", par: ",
             round(cop_Init_best_Copula$par,4))) 
 
print("###################### Family fixed ####################")
for(f in family_to_try){
   print(f)
   ### 7.3 Estim. theta at fixed family
   est_cop_f=determine_theta3(U_pobs_Init, f, pobs_to_calc=FALSE)
   ### Save cop_Init
   assign(paste0("cop_Init_",f), est_cop_f)
   assign(paste0("estim_ci_theta_Init_",f, "_68"),determine_ci_theta2(get(paste0("cop_Init_",f)),
                                                               U_pobs_Init,
                                                               ci_level=0.68))
   assign(paste0("estim_ci_theta_Init_",f, "_95"),determine_ci_theta2(get(paste0("cop_Init_",f)),
                                                               U_pobs_Init,ci_level=0.95))
   print(paste0("Copula Family fixed: ", f, ", ", est_cop_f$family, ", par: ",
             round(est_cop_f$par,4))) 
}

print("###################### Sliding Window ####################")
time_loop=Sys.time() 
k=0
for(l in Label_SlidingWindowX){
  print(Sys.time()-time_loop) 
  time_loop=Sys.time() 
  print(paste0("##################################### ", l, " #################################"))
  print("###################### Best Copula ####################")

  k=k+1
  ### Methodology
  tmp_listX_index1<-listX_index1[[l]]
  tmp_listX_index2<-listX_index2[[l]]
  #### to which proba the quantile associated with probs 0.5 in Ref_SlidingWindow correspond to in
  #Sliding_Window?
  probs_listX_index1[k,]=F_GPD(tmp_listX_index1, quant90_index1, quant_Init_index1)
  probs_listX_index2[k,]=F_GPD(tmp_listX_index2, quant90_index2, quant_Init_index2)
  ### Pseudo obs.
  pobs_listX_index1=pobs(tmp_listX_index1)#F_GPD(tmp_listX_index1, quant90_index1, tmp_listX_index1)
  pobs_listX_index2=pobs(tmp_listX_index2)#F_GPD(tmp_listX_index2, quant90_index2, tmp_listX_index2)
  U_pobs=matrix(NaN,ncol=2,nrow=length(pobs_listX_index1))
  U_pobs[,1]=pobs_listX_index1
  U_pobs[,2]=pobs_listX_index2
  #### Selection among Archimedean family (without rotations, i.e. survival)
  cop_best_Copula=selectedCopula3(U_pobs,pobs_to_calc=FALSE)
  res_name_best_Copula[k]<- cop_best_Copula$family
  best_Copula_name_family=cop_best_Copula$family
  ### Goodness of fit test
  test=testCopula3(U_pobs,
                    best_Copula_name_family, par_=cop_best_Copula$par,pobs_to_calc=FALSE)
  #H0: the fit is good. p<0.05 we reject
  res_Goftest_best_Copula[k] <-test$p.value
  print(paste0("Best Copula: ", best_Copula_name_family, ", par: ",
               round(cop_best_Copula$par,4),", GoF test: ",
               round(test$p.value,4))) 
  if(best_Copula_name_family !="Independence"){
    ## Estimate CI of theta_hat 
    estim_ci_theta_listX_68=determine_ci_theta2(cop_best_Copula, U_pobs, ci_level=0.68)
    estim_ci_theta_listX_95=determine_ci_theta2(cop_best_Copula, U_pobs, ci_level=0.95)
    res_theta_best_Copula_68[[l]]<-estim_ci_theta_listX_68
    res_theta_best_Copula_95[[l]]<-estim_ci_theta_listX_95
    ## Proba. Margdep
    array_prob_margdep_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_best_Copula,
                                                                           estim_ci_theta_listX_68,
                                                                           probs_listX_index1[k,],
                                                                           probs_listX_index2[k,],
                                                                           zone_bivar=zone_)
    array_prob_margdep_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_best_Copula,
                                                                           estim_ci_theta_listX_95,
                                                                           probs_listX_index1[k,],
                                                                           probs_listX_index2[k,],
                                                                           zone_bivar=zone_)
    ## Proba. Marg
    array_prob_marg_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_68,
                                                                        probs_listX_index1[k,],
                                                                        probs_listX_index2[k,],
                                                                        zone_bivar=zone_)
    array_prob_marg_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_95,
                                                                        probs_listX_index1[k,],
                                                                        probs_listX_index2[k,],
                                                                        zone_bivar=zone_)
    ## Proba. Dep
    array_prob_dep_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_best_Copula,
                                                                    estim_ci_theta_listX_68,
                                                                    probs_to_eval_index1,
                                                                    probs_to_eval_index2,
                                                                    zone_bivar=zone_)
    array_prob_dep_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_best_Copula,
                                                                    estim_ci_theta_listX_95,
                                                                    probs_to_eval_index1,
                                                                    probs_to_eval_index2,
                                                                    zone_bivar=zone_)
    ## Proba. Inflmarg1
    array_prob_Inflmarg1_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_68,
                                                                        probs_listX_index1[k,],
                                                                        probs_to_eval_index2,
                                                                        zone_bivar=zone_)
    array_prob_Inflmarg1_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_95,
                                                                        probs_listX_index1[k,],
                                                                        probs_to_eval_index2,
                                                                        zone_bivar=zone_)
    ## Proba. Inflmarg2
    array_prob_Inflmarg2_68_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_68,
                                                                        probs_to_eval_index1,
                                                                        probs_listX_index2[k,],
                                                                        zone_bivar=zone_)
    array_prob_Inflmarg2_95_best_Copula[,,k,]<-determine_ci_prob_bivar_array2(cop_Init_best_Copula,
                                                                        estim_ci_theta_Init_best_Copula_95,
                                                                        probs_to_eval_index1,
                                                                        probs_listX_index2[k,],
                                                                        zone_bivar=zone_)
   
  }
  print("###################### Family fixed ####################")

  #### Re-do the whole methodology for each family
  for(f in family_to_try){
   print(paste0("######### ",f, " ########"))
   if(f=="Clayton"){f_number=3}
   if(f=="Gumbel"){f_number=4}
   if(f=="Frank"){f_number=5}
   if(f=="Joe"){f_number=6}
   ### Estim. theta at fixed family
   est_cop_f=determine_theta3(U_pobs, f, pobs_to_calc=FALSE)
   ### Goofness of fit test
   eval(parse(text=paste0("res_Goftest_",f,"[k] <- testCopula3(U_pobs,f,
                          par_=est_cop_f$par, pobs_to_calc=FALSE)$p.value")))
   tmp_cop_f=est_cop_f
   eval(parse(text=paste0("tmp_pval_f <-res_Goftest_",f,"[k]")))
   print(paste0("Best Copula: ", tmp_cop_f$family, ", par: ",
                round(tmp_cop_f$par,4),", GoF test: ",
               round(tmp_pval_f,4))) 
   ## Estimate CI of theta_hat 
   estim_ci_theta_f_68=determine_ci_theta2(tmp_cop_f, U_pobs, ci_level=0.68)
   estim_ci_theta_f_95=determine_ci_theta2(tmp_cop_f, U_pobs, ci_level=0.95)
   eval(parse(text=paste0("res_theta_",f,"_68[['", l, "']] <- estim_ci_theta_f_68")))
   eval(parse(text=paste0("res_theta_",f,"_95[['", l, "']] <- estim_ci_theta_f_95")))
   ##Proba. Margdep
   tmp_prob_margdep_68<-determine_ci_prob_bivar_array2(tmp_cop_f,
                                                      estim_ci_theta_f_68,
                                                      probs_listX_index1[k,],
                                                      probs_listX_index2[k,],
                                                      zone_bivar=zone_)
   tmp_prob_margdep_95<-determine_ci_prob_bivar_array2(tmp_cop_f,
                                                      estim_ci_theta_f_95,
                                                      probs_listX_index1[k,],
                                                      probs_listX_index2[k,],
                                                      zone_bivar=zone_)
   eval(parse(text=paste0("array_prob_margdep_68_", f,"[,,k,]<-
                          tmp_prob_margdep_68")))
   eval(parse(text=paste0("array_prob_margdep_95_", f,"[,,k,]<-
                          tmp_prob_margdep_95")))
   ##Proba. Marg
   tmp_prob_marg_68<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_68")),
                                                probs_listX_index1[k,],
                                                probs_listX_index2[k,],
                                                zone_bivar=zone_)
   tmp_prob_marg_95<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_95")),
                                                probs_listX_index1[k,],
                                                probs_listX_index2[k,],
                                                zone_bivar=zone_)
   eval(parse(text=paste0("array_prob_marg_68_", f,"[,,k,]<-
                          tmp_prob_marg_68")))
   eval(parse(text=paste0("array_prob_marg_95_", f,"[,,k,]<-
                          tmp_prob_marg_95")))
   ##Proba. Dep
   tmp_prob_dep_68<-determine_ci_prob_bivar_array2(tmp_cop_f, estim_ci_theta_f_68, 
                                                probs_to_eval_index1,
                                                probs_to_eval_index2,
                                                zone_bivar=zone_)
   tmp_prob_dep_95<-determine_ci_prob_bivar_array2(tmp_cop_f,
                                                  estim_ci_theta_f_95, 
                                                probs_to_eval_index1,
                                                probs_to_eval_index2,
                                                zone_bivar=zone_)
   eval(parse(text=paste0("array_prob_dep_68_", f,"[,,k,]<-
                          tmp_prob_dep_68")))
   eval(parse(text=paste0("array_prob_dep_95_", f,"[,,k,]<-
                          tmp_prob_dep_95")))
   ##Proba. Inflmarg1
   tmp_prob_Inflmarg1_68<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_68")),
                                                probs_listX_index1[k,],
                                                probs_to_eval_index2,
                                                zone_bivar=zone_)
   tmp_prob_Inflmarg1_95<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_95")),
                                                probs_listX_index1[k,],
                                                probs_to_eval_index2,
                                                zone_bivar=zone_)
   eval(parse(text=paste0("array_prob_Inflmarg1_68_", f,"[,,k,]<-
                          tmp_prob_Inflmarg1_68")))
   eval(parse(text=paste0("array_prob_Inflmarg1_95_", f,"[,,k,]<-
                          tmp_prob_Inflmarg1_95")))

   ##Proba. Inflmarg2
   tmp_prob_Inflmarg2_68<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_68")),
                                                probs_to_eval_index1,
                                                probs_listX_index2[k,],
                                                zone_bivar=zone_)
   tmp_prob_Inflmarg2_95<-determine_ci_prob_bivar_array2(get(paste0("cop_Init_",f)),
                                                get(paste0("estim_ci_theta_Init_",f,"_95")),
                                                probs_to_eval_index1,
                                                probs_listX_index2[k,],
                                                zone_bivar=zone_)
   eval(parse(text=paste0("array_prob_Inflmarg2_68_", f,"[,,k,]<-
                          tmp_prob_Inflmarg2_68")))
   eval(parse(text=paste0("array_prob_Inflmarg2_95_", f,"[,,k,]<-
                          tmp_prob_Inflmarg2_95")))

 }
}


##### best_Copula + Union for CIs: best_Copula_Union 
array_prob_margdep_68_best_Copula_Union <- array_prob_marg_68_best_Copula_Union <-
  array_prob_dep_68_best_Copula_Union <-array(NaN,dim=c(length(probs_to_eval_index1),
                                               length(probs_to_eval_index2),
                  length(Label_SlidingWindowX),3)) 
array_prob_margdep_95_best_Copula_Union <- array_prob_marg_95_best_Copula_Union <-
  array_prob_dep_95_best_Copula_Union <-array(NaN,dim=c(length(probs_to_eval_index1),
                                               length(probs_to_eval_index2),
                  length(Label_SlidingWindowX),3))
version_=c("margdep", "marg", "dep")

for(v in version_){
  tmp_proba_68=get(paste0("array_prob_",v,"_68_best_Copula"))
  tmp_proba_95=get(paste0("array_prob_",v,"_95_best_Copula"))
  eval(parse(text=paste0("array_prob_", v, "_68_best_Copula_Union[,,,1]<-
                          tmp_proba_68[,,,1]")))
  eval(parse(text=paste0("array_prob_", v, "_95_best_Copula_Union[,,,1]<-
                          tmp_proba_95[,,,1]")))
}
#
nb_listX=length(Label_SlidingWindowX)
p_of_index2=1
p_of_index1=1
for(p_of_index2 in seq(1,19,by=1)){
    for(p_of_index1 in seq(1,19, by=1)){
      for(v in version_){
        print(v)
        mat_min_CI_68=matrix(NaN,ncol=nb_listX,nrow=4)
        mat_max_CI_68=matrix(NaN,ncol=nb_listX,nrow=4)
        mat_min_CI_95=matrix(NaN,ncol=nb_listX,nrow=4)
        mat_max_CI_95=matrix(NaN,ncol=nb_listX,nrow=4)
        k=0
        for(f in c("Clayton", "Frank", "Gumbel", "Joe")){
              k=k+1
              print(f)
              tmp_prob_68=get(paste0("array_prob_",v,"_68_",f))
              tmp_prob_95=get(paste0("array_prob_",v,"_95_",f))
              mat_min_CI_68[k,]<-tmp_prob_68[p_of_index2, p_of_index1,,2]
              mat_max_CI_68[k,]<-tmp_prob_68[p_of_index2, p_of_index1,,3]
              mat_min_CI_95[k,]<-tmp_prob_95[p_of_index2, p_of_index1,,2]
              mat_max_CI_95[k,]<-tmp_prob_95[p_of_index2, p_of_index1,,3]
          }
        Union_min_CI_68<-apply(mat_min_CI_68,2,min)
        Union_max_CI_68<-apply(mat_max_CI_68,2,max)
        Union_min_CI_95<-apply(mat_min_CI_95,2,min)
        Union_max_CI_95<-apply(mat_max_CI_95,2,max)
        
        eval(parse(text=paste0("array_prob_", v, "_68_best_Copula_Union[p_of_index2, p_of_index1,,2]<- Union_min_CI_68")))
        eval(parse(text=paste0("array_prob_", v, "_68_best_Copula_Union[p_of_index2,p_of_index1,,3]<- Union_max_CI_68")))
        eval(parse(text=paste0("array_prob_", v, "_95_best_Copula_Union[p_of_index2, p_of_index1,,2]<- Union_min_CI_95")))
        eval(parse(text=paste0("array_prob_", v, "_95_best_Copula_Union[p_of_index2,p_of_index1,,3]<- Union_max_CI_95")))
      }
   }
}

#### end best_Copula Union for CIs: 



name_array=c()
for(f in c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton",
           "best_Copula_Union")){
  for(v in c("margdep", "marg", "dep")){
    for(ci in c("68", "95")){
    name_array=c(name_array,paste0("array_prob_",v,"_",ci,"_", f))
    }
  }
}

### For Inflmarg1 and 2
for(f in c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton")){
  for(v in c("Inflmarg1", "Inflmarg2")){
    for(ci in c("68", "95")){
      name_array=c(name_array,paste0("array_prob_",v,"_",ci,"_", f))
    }
  }
}

for(f in c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton")){
  name_array=c(name_array, paste0("res_Goftest_",f),
               paste0("res_theta_",f,"_68"), paste0("res_theta_",f,"_95"))
}


### End exact copy

object_to_save=c(name_array, 
                 "res_name_best_Copula", 
                 "probs_to_eval_index1",
                 "probs_to_eval_index2", 
                 "probs_listX_index1",
                 "probs_listX_index2", 
                 "quant_Init_index1",
                 "quant_Init_index2", 
                 "quant90_index1",
                 "quant90_index2", 
                 "listX_index1", 
                 "listX_index2",
                 "Label_SlidingWindowX"
                 #"uncond_listX_index1",
                 #"uncond_listX_index2"
                 )

name_to_save=c()
for(obj in object_to_save){
  name_to_save=c(name_to_save, paste0(obj, "_CDFt_pooled_ModRef", ModRef))
  assign(paste0(obj, "_CDFt_pooled_ModRef", ModRef), get(paste0(obj)))
}


save_wd=paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/",
               Region, "/", name_var1, "_", name_var2, "/SlidingWindow30_Ref",
               Label_SlidingWindowX[[coord_Ref_SlidingWindow]],"/Pooled_Res/")
setwd(paste0(save_wd))

save(list=name_to_save,
     file=paste0("Res_array_prob_list", Length_SlidingWindow, "_Ref",
     Label_SlidingWindowX[[coord_Ref_SlidingWindow]], "_", name_var1, "_",
     name_var2, "_CDFt_pooled_ModRef", ModRef,
                 "_", season, "_", period, "_", Region, ".RData"))

setwd(paste0(save_wd))
q()










