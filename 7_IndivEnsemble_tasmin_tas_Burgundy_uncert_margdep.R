#!/usr/bin/env Rscript
######### Version with uncertainty of marginals
rm(list=ls())
args <- commandArgs(trailingOnly=TRUE)
Mod=args[1]
name_month_Frost=args[2]
Ref_SlidingWindow=args[3]
#Mod="CNRMCM6"
#name_month_Frost="FrostApril"
#Ref_SlidingWindow="1871_1900"
print(Mod)
print(name_month_Frost)
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
library(copula)
library(VineCopula)
library(viridis)
library(RColorBrewer)
library(scales)
load("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Compound_Event_Temporal_indices_SlidingWindows_1850_2100.RData")

name_var1="tasmin" ### for Frost
name_var2="tas" ### for GDD
Region="Burgundy"
is_TardiveFrost=TRUE
Length_SlidingWindow=30
if(is_TardiveFrost==TRUE){zone_="topleft"}
lab_index1="TasminApr"
lab_index2="Spat. mean GDD"


setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/",
             Region, "/", name_var1, "_", name_var2, "/SlidingWindow",
             Length_SlidingWindow, "_", name_month_Frost, "_Ref",
             Ref_SlidingWindow))


load(paste0("Res_array_prob_list30_", name_month_Frost,"_", name_var1, "_",
            name_var2, "_", Mod, "_1850_2100_", Region, ".RData"))




if(Mod %in% c("CNRMCM6", "IPSL", "MIROC6", "INMCM48", "CanESM5",
              "GFDLCM4", "MPIESM1LR", "MRIESM2", "CNRMCM6HR", "CMCCESM2",
              "ECEARTH3", "FGOALSG3", "INMCM50")){
  period="1850_2100"
  nb_years_=251
  Label_SlidingWindowX= get(paste0("Label_SlidingWindow", Length_SlidingWindow,"_", period))
  coord_Ref_SlidingWindow=which(Label_SlidingWindowX==Ref_SlidingWindow)
}

nb_listX=length(Label_SlidingWindowX)

coord_sub_ts=c(coord_Ref_SlidingWindow:222)

### Determine the best family
best_family_final=list()
for(Mod in Mod){
  seq_best_copula<-get(paste0("res_name_best_Copula_", Mod))[coord_sub_ts]

  best_family_Mod=names(which.max(table(seq_best_copula)))
  res_Gof=round(mean(get(paste0("res_Goftest_", best_family_Mod, "_", 
                                Mod))[coord_sub_ts]<0.05),3)
  if(res_Gof>=0.05){
    print("!!!! check for other family !!!")
    tmp_res_Gof=res_Gof
    ordered_family=sort(table(seq_best_copula),decreasing=TRUE)
    k=1
    while(k<5 && tmp_res_Gof>=0.05){
      tmp_best_family_Mod= names(ordered_family[k])
      tmp_res_Gof=round(mean(get(paste0("res_Goftest_", tmp_best_family_Mod,
                                        "_",
                                                               Mod))[coord_sub_ts]<0.05),3)
      k=k+1
    }
    if(k==5){
      best_family_Mod="No family"
      res_Gof=NaN
    }else{
      best_family_Mod=tmp_best_family_Mod
      res_Gof=tmp_res_Gof
    }
  }
  print(paste0("Best family for ", Mod, ": ", best_family_Mod, ", Gof test prop.: ", res_Gof)) 
  best_family_final[[Mod]]=best_family_Mod
}




print("###################### Sliding Window ####################")
time_loop=Sys.time() 
k=0
### For Bootstrap
B=100
b_probs_68_listX_index1=b_probs_68_listX_index2=array(NaN,
                                                dim=c(length(Label_SlidingWindowX),19,B))
b_probs_95_listX_index1=b_probs_95_listX_index2=array(NaN,
                                                dim=c(length(Label_SlidingWindowX),19,B))

#quant90_index1<-get(paste0("quant90_index1_", Mod))
#quant90_index2<-get(paste0("quant90_index2_", Mod))
quant_Init_index1<-get(paste0("quant_Init_index1_", Mod))
quant_Init_index2<-get(paste0("quant_Init_index2_", Mod))
probs_listX_index1<-get(paste0("probs_listX_index1_", Mod))
probs_listX_index2<-get(paste0("probs_listX_index2_", Mod))
 
for(v in c("margdep", "marg", "dep")){
  for(ci in c(68,95)){
    for(f in best_family_final[[Mod]]){
      assign(paste0("tmp_prob_", v, "_v2bis_", ci),  array(NaN,dim=c(19,19,
                                                                         length(Label_SlidingWindowX),3)))
   }
  }
}

l="1850_1879"
### Init Bootstrap for coord_Ref_SlidingWindow
tmp_listX_index1<-get(paste0("listX_index1_", Mod))[[coord_Ref_SlidingWindow]]
tmp_listX_index2<-get(paste0("listX_index2_", Mod))[[coord_Ref_SlidingWindow]]

b_probs_68_listX_index1[coord_Ref_SlidingWindow,,]=boot_param_GEV_for_min(tmp_listX_index1, quant_Init_index1,B, level_=0.68)
b_probs_68_listX_index2[coord_Ref_SlidingWindow,,]=boot_param_Gaussian(tmp_listX_index2, quant_Init_index2, B, level_=0.68)
b_probs_95_listX_index1[coord_Ref_SlidingWindow,,]=boot_param_GEV_for_min(tmp_listX_index1,
                                                                          quant_Init_index1,B,
                                                                          level_=0.95)
b_probs_95_listX_index2[coord_Ref_SlidingWindow,,]=boot_param_Gaussian(tmp_listX_index2, quant_Init_index2, B, level_=0.95)



for(l in Label_SlidingWindowX){
  print(Sys.time()-time_loop) 
  time_loop=Sys.time() 
  print(paste0("##################################### ", l, " #################################"))
  print("###################### Best Copula ####################")

  k=k+1
  ### Methodology
  tmp_listX_index1<-get(paste0("listX_index1_", Mod))[[l]]
  tmp_listX_index2<-get(paste0("listX_index2_", Mod))[[l]]

  #### Re-do the whole methodology for each family
  for(f in best_family_final[[Mod]]){#family_to_try){
   print(paste0("######### ",f, " ########"))
   if(f=="Clayton"){f_number=3}
   if(f=="Gumbel"){f_number=4}
   if(f=="Frank"){f_number=5}
   if(f=="Joe"){f_number=6}
 
   b_prob_margdep_68<-array(NaN,dim=c(19,19,3,B))
   b_prob_marg_68<-array(NaN,dim=c(19,19,3,B))
   b_prob_dep_68<-array(NaN,dim=c(19,19,3,B))
   ess_to_del<-array(NaN,dim=c(19,19,3,B))
   b_prob_margdep_95<-array(NaN,dim=c(19,19,3,B))
   b_prob_marg_95<-array(NaN,dim=c(19,19,3,B))
   b_prob_dep_95<-array(NaN,dim=c(19,19,3,B))

   if(k!=coord_Ref_SlidingWindow){
    b_probs_68_listX_index1[k,,]=boot_param_GEV_for_min(tmp_listX_index1, quant_Init_index1,B, level_=0.68)
    b_probs_68_listX_index2[k,,]=boot_param_Gaussian(tmp_listX_index2, quant_Init_index2, B, level_=0.68)
    b_probs_95_listX_index1[k,,]=boot_param_GEV_for_min(tmp_listX_index1, quant_Init_index1,B, level_=0.95)
    b_probs_95_listX_index2[k,,]=boot_param_Gaussian(tmp_listX_index2, quant_Init_index2, B, level_=0.95)
   }

   tmp_cop_f=list(family=f)
   tmp_theta_Init_f_68=get(paste0("res_theta_", f, "_68_", Mod))[[coord_Ref_SlidingWindow]]
   tmp_theta_Init_f_95=get(paste0("res_theta_", f, "_95_", Mod))[[coord_Ref_SlidingWindow]]
   tmp_theta_l_f_68=get(paste0("res_theta_", f, "_68_", Mod))[[l]]
   tmp_theta_l_f_95=get(paste0("res_theta_", f, "_95_", Mod))[[l]]
   for(b in 1:B){
    print(b)
    b_prob_margdep_68[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_l_f_68[1],
                                                  b_probs_68_listX_index1[k,,b],
                                                  b_probs_68_listX_index2[k,,b],
                                                  zone_bivar=zone_)
    b_prob_marg_68[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_Init_f_68[1],
                                                  b_probs_68_listX_index1[k,,b],
                                                  b_probs_68_listX_index2[k,,b],
                                                  zone_bivar=zone_)
    b_prob_dep_68[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_l_f_68[1],
                                                  b_probs_68_listX_index1[coord_Ref_SlidingWindow,,b],
                                                  b_probs_68_listX_index2[coord_Ref_SlidingWindow,,b],
                                                  zone_bivar=zone_)
    b_prob_margdep_95[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_l_f_95[1],
                                                  b_probs_95_listX_index1[k,,b],
                                                  b_probs_95_listX_index2[k,,b],
                                                  zone_bivar=zone_)
    b_prob_marg_95[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_Init_f_95[1],
                                                  b_probs_95_listX_index1[k,,b],
                                                  b_probs_95_listX_index2[k,,b],
                                                  zone_bivar=zone_)
    b_prob_dep_95[,,1,b]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                  tmp_theta_l_f_95[1],
                                                  b_probs_95_listX_index1[coord_Ref_SlidingWindow,,b],
                                                  b_probs_95_listX_index2[coord_Ref_SlidingWindow,,b],
                                                  zone_bivar=zone_)
  }
      
   ### For ci 68
   #for each proba, Which bootstrap is responsible of the 84% percentile? 16% percentile?
   if(zone_=="topright"){
    coord_b_sup_margdep_68=apply(b_prob_margdep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_inf_margdep_68=apply(b_prob_margdep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
    coord_b_sup_marg_68=apply(b_prob_marg_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_inf_marg_68=apply(b_prob_marg_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
    coord_b_sup_dep_68=apply(b_prob_dep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_inf_dep_68=apply(b_prob_dep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
  }
   if(zone_=="topleft"){
    coord_b_inf_margdep_68=apply(b_prob_margdep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_sup_margdep_68=apply(b_prob_margdep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
    coord_b_inf_marg_68=apply(b_prob_marg_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_sup_marg_68=apply(b_prob_marg_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
    coord_b_inf_dep_68=apply(b_prob_dep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.84)])[1])})
    coord_b_sup_dep_68=apply(b_prob_dep_68[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.16)])[1])})
  }


   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
     tmp_prob_margdep_v2bis_68[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f, tmp_theta_l_f_68[3],
                                prob_vector1=b_probs_68_listX_index1[k,p_i1,unlist(coord_b_sup_margdep_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[k,p_i2,
                                                                  unlist(coord_b_sup_margdep_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_margdep_v2bis_68[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_68[2],
                                prob_vector1=b_probs_68_listX_index1[k,p_i1,unlist(coord_b_inf_margdep_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[k,p_i2,
                                                                  unlist(coord_b_inf_margdep_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_margdep_v2bis_68[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_68[1],
                                prob_vector1=probs_listX_index1[k,p_i1],
                                prob_vector2=probs_listX_index2[k,p_i2],
                                zone_bivar=zone_)
     }
   }
   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
     tmp_prob_marg_v2bis_68[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                   tmp_theta_Init_f_68[3],
                                prob_vector1=b_probs_68_listX_index1[k,p_i1,unlist(coord_b_sup_marg_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[k,p_i2,
                                                                  unlist(coord_b_sup_marg_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_marg_v2bis_68[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_Init_f_68[2],
                                prob_vector1=b_probs_68_listX_index1[k,p_i1,unlist(coord_b_inf_marg_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[k,p_i2,
                                                                  unlist(coord_b_inf_marg_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_marg_v2bis_68[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_Init_f_68[1],
                                prob_vector1=probs_listX_index1[k,p_i1],
                                prob_vector2=probs_listX_index2[k,p_i2],
                                zone_bivar=zone_)
     }
   }
   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
      tmp_prob_dep_v2bis_68[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f, tmp_theta_l_f_68[3],
                                prob_vector1=b_probs_68_listX_index1[coord_Ref_SlidingWindow,p_i1,unlist(coord_b_sup_dep_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[coord_Ref_SlidingWindow,p_i2,
                                                                  unlist(coord_b_sup_dep_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_dep_v2bis_68[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_68[2],
                                prob_vector1=b_probs_68_listX_index1[coord_Ref_SlidingWindow,p_i1,unlist(coord_b_inf_dep_68[p_i2,p_i1])],
                                prob_vector2=b_probs_68_listX_index2[coord_Ref_SlidingWindow,p_i2,
                                                                  unlist(coord_b_inf_dep_68[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_dep_v2bis_68[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_68[1],
                                prob_vector1=probs_listX_index1[coord_Ref_SlidingWindow,p_i1],
                                prob_vector2=probs_listX_index2[coord_Ref_SlidingWindow,p_i2],
                                zone_bivar=zone_)
     }
   }
   ### For ci 95: 
   #for each proba, Which bootstrap is responsible of the 84% percentile? 16% percentile?
   if(zone_=="topright"){
    coord_b_sup_margdep_95=apply(b_prob_margdep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_inf_margdep_95=apply(b_prob_margdep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
    coord_b_sup_marg_95=apply(b_prob_marg_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_inf_marg_95=apply(b_prob_marg_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
    coord_b_sup_dep_95=apply(b_prob_dep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_inf_dep_95=apply(b_prob_dep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
   } 
   if(zone_=="topleft"){
    coord_b_inf_margdep_95=apply(b_prob_margdep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_sup_margdep_95=apply(b_prob_margdep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
    coord_b_inf_marg_95=apply(b_prob_marg_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_sup_marg_95=apply(b_prob_marg_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
    coord_b_inf_dep_95=apply(b_prob_dep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.97)])[1])})
    coord_b_sup_dep_95=apply(b_prob_dep_95[,,1,],c(1,2),
                      function(x){return(which(x==sort(x)[trunc(B*0.03)])[1])})
   }

   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
     tmp_prob_margdep_v2bis_95[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f, tmp_theta_l_f_95[3],
                                prob_vector1=b_probs_95_listX_index1[k,p_i1,unlist(coord_b_sup_margdep_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[k,p_i2,
                                                                  unlist(coord_b_sup_margdep_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_margdep_v2bis_95[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_95[2],
                                prob_vector1=b_probs_95_listX_index1[k,p_i1,unlist(coord_b_inf_margdep_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[k,p_i2,
                                                                  unlist(coord_b_inf_margdep_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_margdep_v2bis_95[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_95[1],
                                prob_vector1=probs_listX_index1[k,p_i1],
                                prob_vector2=probs_listX_index2[k,p_i2],
                                zone_bivar=zone_)
     }
   }
   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
     tmp_prob_marg_v2bis_95[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                   tmp_theta_Init_f_95[3],
                                prob_vector1=b_probs_95_listX_index1[k,p_i1,unlist(coord_b_sup_marg_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[k,p_i2,
                                                                  unlist(coord_b_sup_marg_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_marg_v2bis_95[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_Init_f_95[2],
                                prob_vector1=b_probs_95_listX_index1[k,p_i1,unlist(coord_b_inf_marg_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[k,p_i2,
                                                                  unlist(coord_b_inf_marg_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_marg_v2bis_95[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_Init_f_95[1],
                                prob_vector1=probs_listX_index1[k,p_i1],
                                prob_vector2=probs_listX_index2[k,p_i2],
                                zone_bivar=zone_)
     }
   }
   for(p_i2 in 1:19){
     for(p_i1 in 1:19){
      tmp_prob_dep_v2bis_95[p_i2,p_i1,k,3]<-determine_bootstrap_prob_bivar(tmp_cop_f, tmp_theta_l_f_95[3],
                                prob_vector1=b_probs_95_listX_index1[coord_Ref_SlidingWindow,p_i1,unlist(coord_b_sup_dep_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[coord_Ref_SlidingWindow,p_i2,
                                                                  unlist(coord_b_sup_dep_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_dep_v2bis_95[p_i2,p_i1,k,2]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_95[2],
                                prob_vector1=b_probs_95_listX_index1[coord_Ref_SlidingWindow,p_i1,unlist(coord_b_inf_dep_95[p_i2,p_i1])],
                                prob_vector2=b_probs_95_listX_index2[coord_Ref_SlidingWindow,p_i2,
                                                                  unlist(coord_b_inf_dep_95[p_i2,p_i1])],
                                zone_bivar=zone_)
      tmp_prob_dep_v2bis_95[p_i2,p_i1,k,1]<-determine_bootstrap_prob_bivar(tmp_cop_f,
                                                                                 tmp_theta_l_f_95[1],
                                prob_vector1=probs_listX_index1[coord_Ref_SlidingWindow,p_i1],
                                prob_vector2=probs_listX_index2[coord_Ref_SlidingWindow,p_i2],
                                zone_bivar=zone_)
     }
   }
  }
}

#### On inverse low_ci and high_ci for zone_bivar="topleft"
if(zone_=="topleft"){
  for(v in c("margdep", "marg", "dep")){
    for(ci in c("68", "95")){
      tmp_save_low=get(paste0("tmp_prob_",v,"_v2bis_", ci))[,,, 2] 
      tmp_save_high=get(paste0("tmp_prob_",v,"_v2bis_", ci))[,,, 3] 
      eval(parse(text=(paste0("tmp_prob_", v, "_v2bis_", ci,"[,,,2]<-tmp_save_high"))))
      eval(parse(text=(paste0("tmp_prob_", v, "_v2bis_", ci,"[,,,3]<-tmp_save_low"))))
    }
  }
}


 


name_array=c()
for(f in c(best_family_final[[Mod]])){#c("best_Copula", "Frank", "Gumbel", "Joe", "Clayton",
#           "best_Copula_Union")){
  for(v in c("margdep", "marg", "dep")){
    for(ci in c("68", "95")){#, "95")){
      eval(parse(text=(paste0("array_prob_",v,"_v2bis_",ci,"_", f,"<-tmp_prob_", v,
                    "_v2bis_", ci))))
      name_array=c(name_array,paste0("array_prob_",v,"_v2bis_",ci,"_", f))
    }
  }
}

object_to_save=c(name_array)

name_to_save=c()
for(obj in object_to_save){
  name_to_save=c(name_to_save, paste0(obj, "_marg_uncert_", Mod))
  assign(paste0(obj, "_marg_uncert_", Mod), get(paste0(obj)))
}
save_wd=getwd()
setwd(paste0(save_wd,"/marg_uncert/"))
save(list=name_to_save,
     file=paste0("Res_array_v2bis_marg_uncert_B",B,"_prob_list", Length_SlidingWindow, "_Ref",
     Label_SlidingWindowX[[coord_Ref_SlidingWindow]], "_", name_var1, "_", name_var2, "_", Mod,
                 "_", period, "_", Region, ".RData"))
q()
#
#   
