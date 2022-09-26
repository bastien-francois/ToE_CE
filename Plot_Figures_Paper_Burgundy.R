rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

rdylgnPal<-colorRampPalette(brewer.pal(11, "RdYlGn"))(500)
brbgPal <- rev(colorRampPalette(brewer.pal(11, "BrBG"))(500))
rbPal<-rev(colorRampPalette(brewer.pal(11, "RdBu"))(500))
SpecPaldisc13<-c(rbPal[1], rev(brewer.pal(n=11, name="Spectral")), rbPal[500])

Ref_SlidingWindow="1871_1900"
Region="Burgundy"
x_axis_lab="Minimal temperature [°C]"
y_axis_lab="GDD [°C.d]"

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Burgundy/tasmin_tas/SlidingWindow30_FrostApril_Ref",Ref_SlidingWindow))
load(paste0("Res_array_prob_list30_FrostApril_tasmin_tas_CNRMCM6_1850_2100_Burgundy.RData"))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Burgundy/tasmin_tas/SlidingWindow30_FrostApril_Ref",Ref_SlidingWindow,"/marg_uncert/"))
load(paste0("Res_array_v2bis_marg_uncert_B100_prob_list30_Ref",Ref_SlidingWindow,"_tasmin_tas_CNRMCM6_1850_2100_Burgundy.RData"))

coord_sub_ts=c(which(Ref_SlidingWindow==Label_SlidingWindowX_CNRMCM6):222)
coord_sub_ts
col_main=c("dodgerblue","darkorange", "gray", "indianred1") # "chartreuse3"

# setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Burgundy/tasmin_tas/fastplot_v2bis/SlidingWindow30_Ref", Ref_SlidingWindow))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Plot_Figures"))

coord_baseline=22
print(coord_baseline)
period_baseline=paste0(substrLeft(Label_SlidingWindowX_CNRMCM6[coord_baseline[1]],4), "_", substrRight(Label_SlidingWindowX_CNRMCM6[coord_baseline[length(coord_baseline)]],4))
nb_listX=201
xlab_name=seq.int(1886,2086,20)
xlab_pos=seq(1,nb_listX,20)
p_i2=5
p_i1=19

#### Figure 1
pdf(paste0("fastplot_", Region, "_Figure1_",period_baseline ,".pdf"), width=15, height=10)#,width=1333, height=666)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
# ## Proba
fastplot_time_series(array_prob_margdep_v2bis_68_Clayton_marg_uncert_CNRMCM6, array_prob_margdep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_margdep_68_Clayton_CNRMCM6, array_prob_margdep_95_Clayton_CNRMCM6,ylim_=c(-0.03, 0.35),
                     main="Marg.-dep.",coord_baseline,coord_sub_ts, subplot_name="(a)")
fastplot_time_series(array_prob_marg_v2bis_68_Clayton_marg_uncert_CNRMCM6, array_prob_marg_v2bis_95_Clayton_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_marg_68_Clayton_CNRMCM6,array_prob_marg_95_Clayton_CNRMCM6,ylim_=c(-0.03, 0.35),main="Marg.",coord_baseline,coord_sub_ts, subplot_name="(b)")
fastplot_time_series(array_prob_dep_v2bis_68_Clayton_marg_uncert_CNRMCM6, array_prob_dep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_dep_68_Clayton_CNRMCM6,array_prob_dep_95_Clayton_CNRMCM6,ylim_=c(-0.03, 0.35),main="Dep.",coord_baseline,coord_sub_ts, subplot_name="(c)")
plot.new()
## Contrib
res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
                                              array_prob_marg_v2bis_68_Clayton_marg_uncert_CNRMCM6,
                                              array_prob_dep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
                                              coord_baseline)
fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=TRUE, plot_matrix=FALSE, plot_barplot=TRUE,coord_sub_ts)
dev.off()



source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
# #### Figure 2
# pdf(paste0("fastplot_", Region, "_Figure2_68_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# ## Proba
# fastplot_toe_mat(array_prob_margdep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Marg.-dep.",coord_baseline,coord_sub_ts)
# fastplot_toe_mat(array_prob_marg_v2bis_68_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Marg.",coord_baseline,coord_sub_ts)
# fastplot_toe_mat(array_prob_dep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Dep.",coord_baseline,coord_sub_ts)
# ## Contrib
# res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
#                                               array_prob_marg_v2bis_68_Clayton_marg_uncert_CNRMCM6, 
#                                               array_prob_dep_v2bis_68_Clayton_marg_uncert_CNRMCM6,
#                                               coord_baseline)
# fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=FALSE, plot_matrix=TRUE, plot_barplot=FALSE,coord_sub_ts)
# dev.off()
# 
# pdf(paste0("fastplot_", Region, "_Figure2_95_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# ## Proba
# fastplot_toe_mat(array_prob_margdep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Marg.-dep.",coord_baseline,coord_sub_ts)
# fastplot_toe_mat(array_prob_marg_v2bis_95_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Marg.",coord_baseline,coord_sub_ts)
# fastplot_toe_mat(array_prob_dep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
#                  p_i2, p_i1, main="Dep.",coord_baseline,coord_sub_ts)
# ## Contrib
# res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
#                                               array_prob_marg_v2bis_95_Clayton_marg_uncert_CNRMCM6, 
#                                               array_prob_dep_v2bis_95_Clayton_marg_uncert_CNRMCM6,
#                                               coord_baseline)
# fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=FALSE, plot_matrix=TRUE, plot_barplot=FALSE,coord_sub_ts)
# dev.off()


### Figure 3
rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

Ref_SlidingWindow="1871_1900"
Region="Burgundy"
x_axis_lab="Minimal temperature [°C]"
y_axis_lab="GDD [°C.d]"

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Burgundy/tasmin_tas/SlidingWindow30_FrostApril_Ref",Ref_SlidingWindow))
list_Mod=c("CNRMCM6", "CNRMCM6HR", "IPSL", "MIROC6", "ECEARTH3", "INMCM48", "INMCM50", "CanESM5",
           "GFDLCM4", "MPIESM1LR", "MRIESM2", "CMCCESM2", "FGOALSG3")
for(Mod in list_Mod){
  print(Mod)
  load(paste0("Res_array_prob_list30_FrostApril_tasmin_tas_", Mod, "_1850_2100_Burgundy.RData"))
}

coord_sub_ts=c(which(Ref_SlidingWindow==Label_SlidingWindowX_CNRMCM6):222)

### Determine the best family
best_family_final=list()
for(Mod in list_Mod){
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
      tmp_res_Gof=round(mean(get(paste0("res_Goftest_", tmp_best_family_Mod, "_", Mod))[coord_sub_ts]<0.05),3)
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


plot_v2bis=TRUE
rm(list=setdiff(ls(), c("best_family_final","Label_SlidingWindowX_CNRMCM6",
                        "quant_Init_index1_CNRMCM6","quant_Init_index2_CNRMCM6",
                        "probs_to_eval_index1_CNRMCM6","probs_to_eval_index2_CNRMCM6",
                        "plot_v2bis","coord_sub_ts","Ref_SlidingWindow", "Region", "x_axis_lab", "y_axis_lab")))
#### Indiv
list_array_Indiv_margdep_68=list()
list_array_Indiv_marg_68=list()
list_array_Indiv_dep_68=list()
list_array_Indiv_margdep_95=list()
list_array_Indiv_marg_95=list()
list_array_Indiv_dep_95=list()
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Burgundy/tasmin_tas/SlidingWindow30_FrostApril_Ref",Ref_SlidingWindow,"/marg_uncert/"))
list_Mod=c("CNRMCM6", "CanESM5","ECEARTH3", "INMCM48", "INMCM50", "CNRMCM6HR",
           "GFDLCM4", "MPIESM1LR", "MRIESM2", "FGOALSG3", "IPSL", "MIROC6")
for(Mod in list_Mod){
  print(Mod)
  load(paste0("Res_array_v2bis_marg_uncert_B100_prob_list30_Ref", Ref_SlidingWindow,"_tasmin_tas_", Mod, "_1850_2100_Burgundy.RData"))
}
for(Mod in list_Mod){
  print(Mod)
  print(best_family_final[[Mod]])
  list_array_Indiv_margdep_68[[Mod]]<-get(paste0("array_prob_margdep_v2bis_68_",best_family_final[[Mod]],"_marg_uncert_", Mod))
  list_array_Indiv_marg_68[[Mod]]<-get(paste0("array_prob_marg_v2bis_68_",best_family_final[[Mod]],"_marg_uncert_", Mod))
  list_array_Indiv_dep_68[[Mod]]<-get(paste0("array_prob_dep_v2bis_68_",best_family_final[[Mod]],"_marg_uncert_", Mod))
  list_array_Indiv_margdep_95[[Mod]]<-get(paste0("array_prob_margdep_v2bis_95_",best_family_final[[Mod]],"_marg_uncert_", Mod))
  list_array_Indiv_marg_95[[Mod]]<-get(paste0("array_prob_marg_v2bis_95_",best_family_final[[Mod]],"_marg_uncert_", Mod))
  list_array_Indiv_dep_95[[Mod]]<-get(paste0("array_prob_dep_v2bis_95_",best_family_final[[Mod]],"_marg_uncert_", Mod))
}
# setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Burgundy/tasmin_tas/fastplot_v2bis/SlidingWindow30_Ref", Ref_SlidingWindow))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Plot_Figures"))

#### Full
tmp_dir=paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Burgundy/tasmin_tas/SlidingWindow30_FrostApril_Ref", Ref_SlidingWindow, "/")
load(paste0(tmp_dir,"Pooled_Res/Res_array_prob_list30_FrostApril_tasmin_tas_CDFt_pooled_ModRefNone_1850_2100_Burgundy.RData"))

### Determine the best family
best_family_final_Full=list()
for(Mod in "CDFt_pooled_ModRefNone"){
  seq_best_copula<-get(paste0("res_name_best_Copula_", Mod))
  
  best_family_Mod=names(which.max(table(seq_best_copula)))
  res_Gof=round(mean(get(paste0("res_Goftest_", best_family_Mod, "_",
                                Mod))<0.05),3)
  if(res_Gof>=0.05){
    print("!!!! check for other family !!!")
    tmp_res_Gof=res_Gof
    ordered_family=sort(table(seq_best_copula),decreasing=TRUE)
    k=1
    while(k<5 && tmp_res_Gof>=0.05){
      tmp_best_family_Mod= names(ordered_family[k])
      tmp_res_Gof=round(mean(get(paste0("res_Goftest_", tmp_best_family_Mod,
                                        "_",
                                        Mod))<0.05),3)
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
  best_family_final_Full[[Mod]]=best_family_Mod
}


load(paste0(tmp_dir,"Pooled_Res/marg_uncert/Res_array_v2bis_marg_uncert_B100_prob_list30_Ref1871_1900_tasmin_tas_CDFt_pooled_ModRefNone_1850_2100_Burgundy.RData"))
list_array_Full_margdep_68=list()
list_array_Full_marg_68=list()
list_array_Full_dep_68=list()
list_array_Full_margdep_95=list()
list_array_Full_marg_95=list()
list_array_Full_dep_95=list()
for(Mod in "CDFt_pooled_ModRefNone"){
  print(Mod)
  print(best_family_final_Full[[Mod]])
  list_array_Full_margdep_68[[Mod]]<-get(paste0("array_prob_margdep_v2bis_68_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
  list_array_Full_marg_68[[Mod]]<-get(paste0("array_prob_marg_v2bis_68_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
  list_array_Full_dep_68[[Mod]]<-get(paste0("array_prob_dep_v2bis_68_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
  list_array_Full_margdep_95[[Mod]]<-get(paste0("array_prob_margdep_v2bis_95_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
  list_array_Full_marg_95[[Mod]]<-get(paste0("array_prob_marg_v2bis_95_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
  list_array_Full_dep_95[[Mod]]<-get(paste0("array_prob_dep_v2bis_95_",best_family_final_Full[[Mod]],"_marg_uncert_", Mod))
}




rdylgnPal<-colorRampPalette(brewer.pal(11, "RdYlGn"))(500)
brbgPal <- rev(colorRampPalette(brewer.pal(11, "BrBG"))(500))
rbPal<-rev(colorRampPalette(brewer.pal(11, "RdBu"))(500))
SpecPaldisc13<-c(rbPal[1], rev(brewer.pal(n=11, name="Spectral")), rbPal[500])

coord_baseline=22
p_i2=5
p_i1=19
nb_listX=201
xlab_name=seq.int(1886,2086,20)
xlab_pos=seq(1,nb_listX,20)
period_baseline=paste0(substrLeft(Label_SlidingWindowX_CNRMCM6[coord_baseline[1]],4), "_", substrRight(Label_SlidingWindowX_CNRMCM6[coord_baseline[length(coord_baseline)]],4))

pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_68[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5),  main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(f)", subplot_pos="topright")
dev.off()

pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_95[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5), main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_95[[1]], p_i2, p_i1,  NaN,  ylim_=c(0, 0.5),main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_95[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(f)", subplot_pos="topright")
dev.off()

p_i2=1
p_i1=19
pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_68[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5),  main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(f)", subplot_pos="topright")
dev.off()

pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_95[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5), main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_95[[1]], p_i2, p_i1,  NaN,  ylim_=c(0, 0.5),main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_95[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(f)", subplot_pos="topright")
dev.off()

p_i2=9
p_i1=19
pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_68[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5),  main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_68[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="68",subplot_name="(f)", subplot_pos="topright")
dev.off()

pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure3_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(a)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(0, 0.5), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(b)", subplot_pos="topright")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted,  ylim_=c(0, 0.5),main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95",is_TardiveFrost=TRUE,subplot_name="(c)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_margdep_95[[1]], p_i2, p_i1, NaN,  ylim_=c(0, 0.5), main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(d)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_marg_95[[1]], p_i2, p_i1,  NaN,  ylim_=c(0, 0.5),main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(e)", subplot_pos="topright")
fastplot_Full_time_series(list_array_Full_dep_95[[1]], p_i2, p_i1,  NaN, ylim_=c(0, 0.5), main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="95",subplot_name="(f)", subplot_pos="topright")
dev.off()


#### Contribution Indiv Full
p_i2=5
p_i1=19
col_main=c("dodgerblue","darkorange", "gray", "indianred1") # "chartreuse3"
res_Indiv_Contrib=list()
for(Mod in list_Mod){
  print(Mod)
  print(best_family_final[[Mod]])
  res_Indiv_Contrib[[Mod]]<-compute_Contrib_matrix_from_array(list_array_Indiv_margdep_68[[Mod]],
                                                              list_array_Indiv_marg_68[[Mod]], 
                                                              list_array_Indiv_dep_68[[Mod]],
                                                              coord_baseline)
}
res_Full_Contrib=list()
res_Full_Contrib[["CDFt_pooled_ModRefNone"]]<-compute_Contrib_matrix_from_array(list_array_Full_margdep_68[["CDFt_pooled_ModRefNone"]],
                                                                                   list_array_Full_marg_68[["CDFt_pooled_ModRefNone"]], 
                                                                                   list_array_Full_dep_68[["CDFt_pooled_ModRefNone"]],
                                                                                   coord_baseline)
sorted_list_Mod=sort(list_Mod)
mat_count=matrix(NaN,ncol=length(sorted_list_Mod)+2, nrow=3)
rownames(mat_count)<-c("Marg.", "Dep.", "Int.")
colnames(mat_count)<-c(sorted_list_Mod,"Indiv-Ens.", "Full-Ens.")

array_Med_C_marg_Indiv<-array(NaN,dim=c(19,19,length(sorted_list_Mod)))
array_Med_C_dep_Indiv<-array(NaN,dim=c(19,19,length(sorted_list_Mod)))
array_Med_C_int_Indiv<-array(NaN,dim=c(19,19,length(sorted_list_Mod)))

k=0
for(Mod in sorted_list_Mod){
  k=k+1
  tmp_Med_C_marg_Indiv<-apply(res_Indiv_Contrib[[Mod]]$Contrib_marg,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  tmp_Med_C_dep_Indiv<-apply(res_Indiv_Contrib[[Mod]]$Contrib_dep,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  tmp_Med_C_int_Indiv<-apply(res_Indiv_Contrib[[Mod]]$Contrib_int,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
  mat_count[,k]=c(tmp_Med_C_marg_Indiv[p_i2,p_i1], tmp_Med_C_dep_Indiv[p_i2,p_i1],
                  tmp_Med_C_int_Indiv[p_i2,p_i1])
  ### Feed the arrays
  array_Med_C_marg_Indiv[,,k]<-tmp_Med_C_marg_Indiv
  array_Med_C_dep_Indiv[,,k]<-tmp_Med_C_dep_Indiv
  array_Med_C_int_Indiv[,,k]<-tmp_Med_C_int_Indiv
}

k=k+1
mat_count[,k]=apply(mat_count,1,median,na.rm=TRUE)
k=k+1
Med_C_marg_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefNone"]]$Contrib_marg,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
Med_C_dep_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefNone"]]$Contrib_dep,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
Med_C_int_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefNone"]]$Contrib_int,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
mat_count[,k]=c(Med_C_marg_Full[p_i2,p_i1], Med_C_dep_Full[p_i2,p_i1],
                Med_C_int_Full[p_i2,p_i1])


pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure4_Contrib_Barplot_68_",period_baseline ,".pdf"),width=15, height=8)
par(mfrow=c(1,2),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Boxplot_ToE(list_array_Indiv_margdep_68, 
                     list_array_Indiv_marg_68,
                     list_array_Indiv_dep_68,
                     list_array_Full_margdep_68,
                     list_array_Full_marg_68,
                     list_array_Full_dep_68,
                     p_i2, p_i1, array_dotted, ylim_=c(1900,2086),
                     main_="NaN", coord_baseline,coord_sub_ts,
                     label_68_95="68",subplot_name="")
barplot(mat_count, col=c(col_main[1], col_main[2], col_main[3]), beside=TRUE,
        ylab="Median contribution (%)", ylim=c(-5,120), las=2)
legend("topleft",
       legend = c("Marg.", "Dep.", "Int."),
       fill = c(col_main[1], col_main[2], col_main[3]))
legend("topright",
       legend = "",bty="n", cex=2)
abline(h=0)
abline(h=100,lty=2)
abline(h=50,lty=2)
abline(v=16*3+0.5, lty=2)
abline(v=17*3+1.5, lty=2)
dev.off()

pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure4_Contrib_Barplot_95_",period_baseline ,".pdf"),width=15, height=8)
par(mfrow=c(1,2),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Boxplot_ToE(list_array_Indiv_margdep_95, 
                     list_array_Indiv_marg_95,
                     list_array_Indiv_dep_95,
                     list_array_Full_margdep_95,
                     list_array_Full_marg_95,
                     list_array_Full_dep_95,
                     p_i2, p_i1, array_dotted, ylim_=c(1900,2086),
                     main_="NaN", coord_baseline,coord_sub_ts,
                     label_68_95="95",subplot_name="")
barplot(mat_count, col=c(col_main[1], col_main[2], col_main[3]), beside=TRUE,
        ylab="Median contribution (%)", ylim=c(-5,120), las=2)
legend("topleft",
       legend = c("Marg.", "Dep.", "Int."),
       fill = c(col_main[1], col_main[2], col_main[3]))
legend("topright",
       legend = "",bty="n", cex=2)
abline(h=0)
abline(h=100,lty=2)
abline(h=50,lty=2)
abline(v=16*3+0.5, lty=2)
abline(v=17*3+1.5, lty=2)
dev.off()



# pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure5_mat_TOE_68_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, main_="Marg.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_margdep_68, p_i2, p_i1, NaN,  main_="Marg.-dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_marg_68, p_i2, p_i1,  NaN, main_="Marg.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_dep_68, p_i2, p_i1,  NaN, main_="Dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# dev.off()
# 
# pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure5_mat_TOE_95_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, main_="Marg.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_margdep_95, p_i2, p_i1, NaN,  main_="Marg.-dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_marg_95, p_i2, p_i1,  NaN, main_="Marg.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Full_dep_95, p_i2, p_i1,  NaN, main_="Dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# dev.off()

# 
# ### Mat Contribution Indv Full
# pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure6_Contrib_mat_IndivFull_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# fastplot_Contrib_matrix(apply(array_Med_C_marg_Indiv,c(1,2),median),
#                         apply(array_Med_C_dep_Indiv,c(1,2),median),
#                         apply(array_Med_C_int_Indiv,c(1,2),median),
#                         plot_matrix=TRUE,is_TardiveFrost=TRUE)
# fastplot_Contrib_matrix(Med_C_marg_Full,
#                         Med_C_dep_Full,
#                         Med_C_int_Full,
#                         plot_matrix=TRUE,is_TardiveFrost=TRUE)
# dev.off()

p_i2=5
p_i1=19
pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure7_Contrib_ts_IndivFull_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Contrib_Indiv(res_Indiv_Contrib,p_i2,p_i1,coord_sub_ts)
fastplot_Contrib(res_Full_Contrib[["CDFt_pooled_ModRefNone"]],p_i2,p_i1,plot_ts=TRUE, plot_matrix=FALSE, plot_barplot=FALSE,coord_sub_ts)
dev.off()


# source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
# p_i2=5
# p_i1=19
# pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure8_68_Sum_IQR_Indiv_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, 
#                main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, 
#                main_="Marg.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, 
#                main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, 
#                main_="Marg.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68",is_TardiveFrost=TRUE)
# dev.off()
# 
# pdf(paste0("fastplot_", Region, "_", quant_Init_index2_CNRMCM6[p_i2],"_", quant_Init_index1_CNRMCM6[p_i1], "_Figure8_95_Sum_IQR_Indiv_",period_baseline ,".pdf"),width=1000, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
# fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, 
#                main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, 
#                main_="Marg.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, 
#                main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, 
#                main_="Marg.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, main_="Dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95",is_TardiveFrost=TRUE)
# dev.off()








#####################################################################################################

## a la retraite
# fastplot_Indiv<-function(list_array_prob,p_i2, p_i1, array_dotted=NaN, 
#                          main_=NaN, coord_baseline=NaN,coord_sub_ts=NaN, plot_proba=FALSE, plot_toe_mat=FALSE,plot_sum=FALSE, label_68_95="NaN"){
#   toe_red_68=col=rgb(1,0,0, alpha=0.2)
#   toe_red_95=col=rgb(1,0,0, alpha=0.05)
#   if(plot_proba==TRUE){
#     plot(list_array_prob[[1]][p_i2,p_i1,coord_sub_ts,1],type='l',col="white",
#          ylim=c(0,0.2), ylab="Proba.",xlab="Sliding windows", xaxt="n", main=paste0(main_))
#     axis(1,xlab_pos,xlab_name)
#     count_ToE=0
#     indiv_ToE=c()
#     for(i_Mod in names(list_array_prob)){
#       tmp_proba=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
#       lines(tmp_proba, type='l', col=toe_red_68)
#       tmp_ci_low_68=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,2]
#       tmp_ci_high_68=list_array_prob[[i_Mod]][p_i2,p_i1,coord_sub_ts,3]
#       polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(tmp_ci_low_68,rev(tmp_ci_high_68)),
#               col= toe_red_95, border = FALSE)
#       tmp_toe=compute_ToE_matrix_from_array(list_array_prob[[i_Mod]],30, "1850_2100", coord_baseline)
#       if(!is.na(tmp_toe$coord_TOE[p_i2,p_i1])){
#         indiv_ToE=c(indiv_ToE,tmp_toe$TOE[p_i2,p_i1])
#         count_ToE=count_ToE+1
#       }
#       abline(v=tmp_toe$coord_TOE[p_i2,p_i1]-coord_sub_ts[1]+1, col=toe_red_68,lty=1)
#     }
#     polygon(c(coord_baseline,rev(coord_baseline)),c(rep(-1000,length(coord_baseline)),rev(rep(6000,length(coord_baseline)))),
#             col=rgb(0,0,0, alpha=0.1), border = FALSE)
#     toe1=paste0(trunc(median(indiv_ToE))," (", trunc(median(indiv_ToE))-15,"-", trunc(median(indiv_ToE))+14, ")")
#     leg.txt <- c(paste0("Median ToE ", label_68_95,"%: ",  toe1),paste0("Nb. models w. ToE ", label_68_95,"%: ",  count_ToE))
#     legend("topleft",leg.txt,col=c("red",NaN),lty=c(1,NaN))
#     # legend("topleft",leg.txt,col=c(toe_red_68),lty=c(1), lwd=c(2))
#     abline(v=median(indiv_ToE)-1864-coord_sub_ts[1]+1, col="red",lty=1)
#   }
#   if(plot_toe_mat==TRUE){
#     ### ToE-Indiv
#     list_toe=array(NaN, dim=c(19,19,length(names(list_array_prob))))
#     k=0
#     for(i_Mod in names(list_array_prob)){
#       k=k+1
#       list_toe[,,k]=compute_ToE_matrix_from_array(list_array_prob[[i_Mod]],30, "1850_2100", coord_baseline)$TOE
#     }
#     res_median=apply(list_toe,c(1,2),median,na.rm=TRUE)
#     image.plot(1:19, 1:19,postproc_image.plot(res_median),zlim=c(1980,2086), col=rdylgnPal, main=paste0(main_), 
#                ylab="PR", xlab="Wind", xaxt="n", yaxt="n")
#     axis(1,seq(2,18,by=2),round(probs_to_eval_index1_CNRMCM6[seq(2,18,by=2)],2))
#     axis(2,seq(2,18,by=2),round(probs_to_eval_index2_CNRMCM6[seq(2,18,by=2)],2))
#     bool_infTOE=which(postproc_image.plot(res_median)<1980,arr.ind=TRUE)
#     points(bool_infTOE, pch=20, col="red")
#   }
#   if(plot_sum==TRUE){
#     lab.breaks=c(0,1,3,5,7,9,11, 13)
#     breaks = c(0,0.9 , 2.9, 4.9 ,6.9, 8.9 ,10.9, 12.9)
#     res_sum=apply(!is.na(list_toe),c(1,2),sum,na.rm=TRUE)
#     image.plot(1:19, 1:19,postproc_image.plot(res_sum),
#                lab.breaks=lab.breaks,
#                breaks= breaks,
#                # axis.args=list(cex.axis =1, at=breaks, labels= c(0,1,3,5,7,9,11, 13)),
#                col=SpecPaldisc13[7:13], zlim=c(0,14),
#                ylab="PR", xlab="Wind", xaxt="n", yaxt="n", main=paste0("nb models with ToE ", main_))
#     axis(1,seq(2,18,by=2),round(probs_to_eval_index1_CNRMCM6[seq(2,18,by=2)],2))
#     axis(2,seq(2,18,by=2),round(probs_to_eval_index2_CNRMCM6[seq(2,18,by=2)],2))
#   }
# }
# ####### Chantier Ruban
# fastplot_Contrib_Indiv<-function(list_Contrib, p_i2,p_i1, coord_sub_ts=NaN){
#   col_main=c("black","dodgerblue","darkorange", "chartreuse3", "indianred1")
#   col_k=0
#   eval(parse(text=paste0("plot(list_Contrib[['CNRMCM6']]$FAR_margdep[p_i2, p_i1,coord_sub_ts], xaxt='n', ylab='Bivar. FAR', xlab='Sliding windows',type='l',ylim=c(-1,1),col='white')")))
#   axis(1,xlab_pos,xlab_name)
#   abline(h=0)
#   
#   for(v in c("margdep", "marg", "dep","int")){
#     col_k=col_k+1
#     mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
#     k=0
#     for(i_Mod in names(list_Contrib)){
#       k=k+1
#       eval(parse(text=paste0("mat_far[,k]<-list_Contrib[[i_Mod]]$FAR_", v,"[p_i2, p_i1,coord_sub_ts]")))
#     }
#     lines(apply(mat_far,1,median),col=col_main[col_k])
#     if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
#     if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
#     if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
#     if(v=="margdep"){col_second=rgb(0, 0,0, alpha=0.1)}
#     mat_far[is.infinite(mat_far) & mat_far<0]<-(-1000)
#     mat_far[is.infinite(mat_far) & mat_far>0]<-(1000)
#     # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
#     #         col= col_second, border = FALSE)
#   }
#   legend("bottomright", c("Marg.+dep.+int.", "Marg.", "Dep.", "Int."),
#          lty=c(1,1,1, 1), col=c("black", col_main[2], col_main[3], col_main[4]))
#   col_k=0
#   eval(parse(text=paste0("plot(list_Contrib[['CNRMCM6']]$Ecart_relat_margdep[p_i2, p_i1,coord_sub_ts], xaxt='n', ylab='Evol. relat.', xlab='Sliding windows',type='l',ylim=c(-0.5,2),col='white')")))
#   axis(1,xlab_pos,xlab_name)
#   abline(h=0)
#   for(v in c("margdep", "marg", "dep", "int")){
#     col_k=col_k+1
#     mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
#     k=0
#     for(i_Mod in names(list_Contrib)){
#       k=k+1
#       eval(parse(text=paste0("mat_far[,k]<-list_Contrib[[i_Mod]]$Ecart_relat_", v,"[p_i2, p_i1,coord_sub_ts]")))
#     }
#     lines(apply(mat_far,1,median),col=col_main[col_k])
#     if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
#     if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
#     if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
#     if(v=="margdep"){col_second=rgb(0, 0,0, alpha=0.1)}
#     mat_far[is.infinite(mat_far) & mat_far<0]<-(-1000)
#     mat_far[is.infinite(mat_far) & mat_far>0]<-(1000)
#     # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
#     #         col= col_second, border = FALSE)
#   }
#   legend("topleft", c("Evol. marg.+dep.+int.", "Evol. marg.", "Evol. dep.", "Evol. int."),
#          lty=c(1,1,1, 1), col=c("black",col_main[2], col_main[3], col_main[4]))
#   tmp_Contrib<- eval(parse(text=paste0("list_Contrib[['CNRMCM6']]$Contrib_marg[p_i2, p_i1,coord_sub_ts]")))
#   eval(parse(text=paste0("plot(tmp_Contrib, xaxt='n', ylab='Contrib.', xlab='Sliding windows',type='l',ylim=c(-150,150),col='white')")))
#   axis(1,xlab_pos,xlab_name)
#   col_k=1
#   for(v in c("marg", "dep","int")){
#     col_k=col_k+1
#     tmp_Contrib<- eval(parse(text=paste0("list_Contrib[['CNRMCM6']]$Contrib_", v,"[p_i2, p_i1,coord_sub_ts]")))
#     coord_Contrib_high<- which(tmp_Contrib>150)
#     coord_Contrib_low<- which(tmp_Contrib<(-150))
#     tmp_Contrib[(tmp_Contrib>150 | tmp_Contrib<(-150))]<-NaN
#     mat_far=matrix(NaN,ncol=length(names(list_Contrib)), nrow=length(coord_sub_ts))
#     k=0
#     for(i_Mod in names(list_Contrib)){
#       k=k+1
#       tmp_Contrib_Mod=eval(parse(text=paste0("list_Contrib[['",i_Mod, "']]$Contrib_", v,"[p_i2, p_i1,coord_sub_ts]")))
#       # coord_Contrib_high<- which(tmp_Contrib_Mod>150)
#       # coord_Contrib_low<- which(tmp_Contrib_Mod<(-150))
#       # tmp_Contrib_Mod[(tmp_Contrib_Mod>150 | tmp_Contrib_Mod<(-150))]<-NaN
#       # eval(parse(text=paste0("lines(tmp_Contrib_Mod, ,type='l',col=col_main[",col_k,"])")))
#       eval(parse(text=paste0("mat_far[,k]<-tmp_Contrib_Mod")))
#     }
#     lines(apply(mat_far,1,median,na.rm=TRUE), col=col_main[col_k])
#     abline(h=median(apply(mat_far,1,median,na.rm=TRUE),na.rm=TRUE), col=col_main[col_k],lty=2,lwd=2)
#     if(v=="marg"){col_second=rgb(0, 0,1, alpha=0.1)}
#     if(v=="dep"){col_second=rgb(1, 0.3,0, alpha=0.1)}
#     if(v=="int"){col_second=rgb(0.1, 0.8,0.3, alpha=0.1)}
#     # polygon(c(1:length(coord_sub_ts),rev(1:length(coord_sub_ts))),c(apply(mat_far,1,min,na.rm=TRUE),rev(apply(mat_far,1,max,na.rm=TRUE))),
#     #         col= col_second, border = FALSE)
#   }
#   abline(h=-100,lty=2)
#   abline(h=100,lty=2)
#   abline(h=0,lty=2)
#   legend("bottomright", c("Contrib. marg.", "Contrib. dep.", "Contrib. int."),
#          lty=c(1,1, 1), col=c(col_main[2], col_main[3], col_main[4]))
# }
