rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

rdylgnPal<-colorRampPalette(brewer.pal(11, "RdYlGn"))(500)
brbgPal <- rev(colorRampPalette(brewer.pal(11, "BrBG"))(500))
rbPal<-rev(colorRampPalette(brewer.pal(11, "RdBu"))(500))
SpecPaldisc13<-c(rbPal[1], rev(brewer.pal(n=11, name="Spectral")), rbPal[500])
rdylbuPal<-rev(c(rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(250)),colorRampPalette(brewer.pal(9, "YlOrRd"))(250)))[1:250]
col_main=c("dodgerblue","darkorange", "gray", "indianred1") # "chartreuse3"
magma13<-rev(magma(14))[1:13]

Ref_SlidingWindow="1871_1900"
Region="Bretagne"
x_axis_lab="Wind speed [m/s]"
y_axis_lab="Precipitation [mm/d]"

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref",Ref_SlidingWindow))
load(paste0("Res_array_prob_list30_Ref",Ref_SlidingWindow,"_sfcWindmax_pr_CNRMCM6_winter_1850_2100_Bretagne.RData"))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref",Ref_SlidingWindow,"/marg_uncert/"))
load(paste0("Res_array_v2bis_marg_uncert_B100_prob_list30_Ref",Ref_SlidingWindow,"_sfcWindmax_pr_CNRMCM6_winter_1850_2100_Bretagne.RData"))

coord_sub_ts=c(which(Ref_SlidingWindow==Label_SlidingWindowX_CNRMCM6):222)
coord_sub_ts

# setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Bretagne/sfcWindmax_pr/fastplot_v2bis/SlidingWindow30_Ref", Ref_SlidingWindow))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Plot_Figures"))

coord_baseline=22
print(coord_baseline)
period_baseline=paste0(substrLeft(Label_SlidingWindowX_CNRMCM6[coord_baseline[1]],4), "_", substrRight(Label_SlidingWindowX_CNRMCM6[coord_baseline[length(coord_baseline)]],4))
nb_listX=201
xlab_name=seq.int(1886,2086,20)
xlab_pos=seq(1,nb_listX,20)
p_i2=16
p_i1=16

#### Figure 1
pdf(paste0("fastplot_", Region, "_Figure1_",period_baseline ,".pdf"), width=15, height=10)#,width=1333, height=666)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
## Proba
fastplot_time_series(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_margdep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_margdep_68_Joe_CNRMCM6, array_prob_margdep_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Marg.-dep.",coord_baseline,coord_sub_ts, subplot_name="(a)")
fastplot_time_series(array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_marg_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_marg_68_Joe_CNRMCM6,array_prob_marg_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Marg.",coord_baseline,coord_sub_ts, subplot_name="(b)")
fastplot_time_series(array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_dep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_dep_68_Joe_CNRMCM6,array_prob_dep_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Dep.",coord_baseline,coord_sub_ts, subplot_name="(c)")
# plot.new()
## Contrib
res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                                              array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6, 
                                              array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                                              coord_baseline)
fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=TRUE, plot_matrix=FALSE, plot_barplot=FALSE,coord_sub_ts)
dev.off()



source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
#### Figure 2
pdf(paste0("fastplot_", Region, "_Figure2_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
## Proba
fastplot_toe_mat(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Marg.-dep.",coord_baseline,coord_sub_ts,subplot_name="(a)",zlim_=c(2020,2086))
fastplot_toe_mat(array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Marg.",coord_baseline,coord_sub_ts,subplot_name="(b)",zlim_=c(2020,2086))
fastplot_toe_mat(array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Dep.",coord_baseline,coord_sub_ts,subplot_name="(c)",zlim_=c(2020,2086))
## Contrib
res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                                              array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6, 
                                              array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                                              coord_baseline)
fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=FALSE, plot_matrix=TRUE, plot_barplot=FALSE,coord_sub_ts)
dev.off()

pdf(paste0("fastplot_", Region, "_Figure2_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
## Proba
fastplot_toe_mat(array_prob_margdep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Marg.-dep.",coord_baseline,coord_sub_ts,subplot_name="(a)",zlim_=c(2020,2086))
fastplot_toe_mat(array_prob_marg_v2bis_95_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Marg.",coord_baseline,coord_sub_ts,subplot_name="(b)",zlim_=c(2020,2086))
fastplot_toe_mat(array_prob_dep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Dep.",coord_baseline,coord_sub_ts,subplot_name="(c)",zlim_=c(2020,2086))
## Contrib
res_CNRMCM6=compute_Contrib_matrix_from_array(array_prob_margdep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                                              array_prob_marg_v2bis_95_Joe_marg_uncert_CNRMCM6, 
                                              array_prob_dep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                                              coord_baseline)
fastplot_Contrib(res_CNRMCM6,p_i2,p_i1,plot_ts=FALSE, plot_matrix=TRUE, plot_barplot=FALSE,coord_sub_ts)
dev.off()


#### DiffFigure 2
pdf(paste0("fastplot_", Region, "_DiffFigure2_68_",period_baseline ,".pdf"),width=15, height=5)
par(mfrow=c(1,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)#oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
## Proba
fastplot_difftoe_mat(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6, 
                     array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6, 
                     array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6,
                 p_i2, p_i1, main="Marg.-(Marg.-dep.)",coord_baseline,coord_sub_ts)
dev.off()



### Figure 3
rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

Ref_SlidingWindow="1871_1900"
Region="Bretagne"
x_axis_lab="Wind speed [m/s]"
y_axis_lab="Precipitation [mm/d]"

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref",Ref_SlidingWindow))
list_Mod=c("CNRMCM6", "CNRMCM6HR", "IPSL", "MIROC6", "ECEARTH3", "INMCM48", "CanESM5",
           "GFDLCM4", "MPIESM1LR", "MRIESM2", "CMCCESM2", "FGOALSG3")
for(Mod in list_Mod){
  print(Mod)
  load(paste0("Res_array_prob_list30_Ref",Ref_SlidingWindow,"_sfcWindmax_pr_", Mod, "_winter_1850_2100_Bretagne.RData"))
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

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref", Ref_SlidingWindow, "/marg_uncert/"))
list_Mod=c("CNRMCM6", "CanESM5","ECEARTH3", "INMCM48", "CNRMCM6HR",
           "GFDLCM4", "MPIESM1LR", "MRIESM2", "CMCCESM2", "FGOALSG3", "IPSL", "MIROC6")
for(Mod in list_Mod){
  print(Mod)
  load(paste0("Res_array_v2bis_marg_uncert_B100_prob_list30_Ref", Ref_SlidingWindow,"_sfcWindmax_pr_", Mod, "_winter_1850_2100_Bretagne.RData"))
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
# setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Bretagne/sfcWindmax_pr/fastplot_v2bis/SlidingWindow30_Ref", Ref_SlidingWindow))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Plot_Figures"))

#### Full
tmp_dir=paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref", Ref_SlidingWindow, "/")
load(paste0(tmp_dir,"Pooled_Res/Res_array_prob_list30_Ref", Ref_SlidingWindow, "_sfcWindmax_pr_CDFt_pooled_ModRefCNRMCM6_winter_1850_2100_Bretagne.RData"))

### Determine the best family
best_family_final_Full=list()
for(Mod in "CDFt_pooled_ModRefCNRMCM6"){
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


load(paste0(tmp_dir,"Pooled_Res/marg_uncert/Res_array_v2bis_marg_uncert_B100_prob_list30_Ref", Ref_SlidingWindow, "_sfcWindmax_pr_CDFt_pooled_ModRefCNRMCM6_winter_1850_2100_Bretagne.RData"))
list_array_Full_margdep_68=list()
list_array_Full_marg_68=list()
list_array_Full_dep_68=list()
list_array_Full_margdep_95=list()
list_array_Full_marg_95=list()
list_array_Full_dep_95=list()
for(Mod in "CDFt_pooled_ModRefCNRMCM6"){
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
rdylbuPal<-rev(c(rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(250)),colorRampPalette(brewer.pal(9, "YlOrRd"))(250)))
col_main=c("dodgerblue","darkorange", "gray", "indianred1") # "chartreuse3"
magma13<-rev(magma(14))[1:13]

coord_baseline=22
p_i2=16
p_i1=16
nb_listX=201
xlab_name=seq.int(1886,2086,20)
xlab_pos=seq(1,nb_listX,20)
period_baseline=paste0(substrLeft(Label_SlidingWindowX_CNRMCM6[coord_baseline[1]],4), "_", substrRight(Label_SlidingWindowX_CNRMCM6[coord_baseline[length(coord_baseline)]],4))

source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
pdf(paste0("fastplot_", Region, "_Figure3_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2), main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68", subplot_name="(a)")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2), main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68", subplot_name="(b)")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="68", subplot_name="(c)")
fastplot_Full_time_series(list_array_Full_margdep_68[[1]], p_i2, p_i1, NaN, ylim_=c(-0.01, 0.2), main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="68", subplot_name="(d)")
fastplot_Full_time_series(list_array_Full_marg_68[[1]], p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="68", subplot_name="(e)")
fastplot_Full_time_series(list_array_Full_dep_68[[1]], p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="68", subplot_name="(f)")
dev.off()

pdf(paste0("fastplot_", Region, "_Figure3_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95", subplot_name="(a)")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95", subplot_name="(b)")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Dep.", coord_baseline,coord_sub_ts,plot_proba=TRUE,plot_sum=FALSE,label_68_95="95", subplot_name="(c)")
fastplot_Full_time_series(list_array_Full_margdep_95[[1]], p_i2, p_i1, NaN,  ylim_=c(-0.01, 0.2),main_="Marg.-dep.",coord_baseline,coord_sub_ts,label_68_95="95", subplot_name="(d)")
fastplot_Full_time_series(list_array_Full_marg_95[[1]], p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Marg.",coord_baseline,coord_sub_ts,label_68_95="95", subplot_name="(e)")
fastplot_Full_time_series(list_array_Full_dep_95[[1]], p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Dep.",coord_baseline,coord_sub_ts,label_68_95="95", subplot_name="(f)")
dev.off()





#### Contribution Indiv Full
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
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
res_Full_Contrib[["CDFt_pooled_ModRefCNRMCM6"]]<-compute_Contrib_matrix_from_array(list_array_Full_margdep_68[["CDFt_pooled_ModRefCNRMCM6"]],
                                                              list_array_Full_marg_68[["CDFt_pooled_ModRefCNRMCM6"]], 
                                                              list_array_Full_dep_68[["CDFt_pooled_ModRefCNRMCM6"]],
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
Med_C_marg_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefCNRMCM6"]]$Contrib_marg,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
Med_C_dep_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefCNRMCM6"]]$Contrib_dep,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
Med_C_int_Full<-apply(res_Full_Contrib[["CDFt_pooled_ModRefCNRMCM6"]]$Contrib_int,c(1,2), function(x){median(x[coord_sub_ts],na.rm=TRUE)})
mat_count[,k]=c(Med_C_marg_Full[p_i2,p_i1], Med_C_dep_Full[p_i2,p_i1],
                Med_C_int_Full[p_i2,p_i1])

pdf(paste0("fastplot_", Region, "_Figure4_Contrib_Barplot_68_",period_baseline ,".pdf"),width=15, height=8)
par(mfrow=c(1,2),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Boxplot_ToE(list_array_Indiv_margdep_68, 
                     list_array_Indiv_marg_68,
                     list_array_Indiv_dep_68,
                     list_array_Full_margdep_68,
                     list_array_Full_marg_68,
                     list_array_Full_dep_68,
                     p_i2, p_i1, array_dotted,  ylim_=c(1980,2086),
                     main_="NaN", coord_baseline,coord_sub_ts,
                     label_68_95="68",subplot_name="")
barplot(mat_count, col=c(col_main[1], col_main[2], col_main[3]), beside=TRUE,
        ylab="Median contribution (%)", ylim=c(-5,120), las=2)
legend("topleft",
       legend = c("Marg.", "Dep.", "Int."),
       fill = c(col_main[1], col_main[2], col_main[3]))
#legend("topright",
#       legend = "(b)",bty="n", cex=2)
abline(h=0)
abline(h=100,lty=2)
abline(h=50,lty=2)
abline(v=16*3+0.5, lty=2)
abline(v=17*3+1.5, lty=2)
dev.off()

pdf(paste0("fastplot_", Region, "_Figure4_Contrib_Barplot_95_",period_baseline ,".pdf"),width=15, height=8)
par(mfrow=c(1,2),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Boxplot_ToE(list_array_Indiv_margdep_95, 
                     list_array_Indiv_marg_95,
                     list_array_Indiv_dep_95,
                     list_array_Full_margdep_95,
                     list_array_Full_marg_95,
                     list_array_Full_dep_95,
                     p_i2, p_i1, array_dotted, ylim_=c(1980,2086),
                     main_="NaN", coord_baseline,coord_sub_ts,
                     label_68_95="95",subplot_name="")
barplot(mat_count, col=c(col_main[1], col_main[2], col_main[3]), beside=TRUE,
        ylab="Median contribution (%)", ylim=c(-5,120), las=2)
legend("topleft",
       legend = c("Marg.", "Dep.", "Int."),
       fill = c(col_main[1], col_main[2], col_main[3]))
#legend("topright",
#       legend = "(b)",bty="n", cex=2)
abline(h=0)
abline(h=100,lty=2)
abline(h=50,lty=2)
abline(v=16*3+0.5, lty=2)
abline(v=17*3+1.5, lty=2)
dev.off()


quant_Init_index1_CNRMCM6=probs_to_eval_index1_CNRMCM6
quant_Init_index2_CNRMCM6=probs_to_eval_index2_CNRMCM6
x_axis_lab="Conditional proba. (Wind speed)"
y_axis_lab="Conditional proba. (Precipitation)"
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

pdf(paste0("fastplot_", Region, "_Figure5_mat_TOE_68_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(a)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(b)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(c)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_margdep_68, p_i2, p_i1, NaN,  ylim_=c(-0.01, 0.2),main_="Marg.-dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(d)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_marg_68, p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Marg.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(e)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_dep_68, p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2), main_="Dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="68",subplot_name="(f)",zlim_=c(1980,2086))
dev.off()

pdf(paste0("fastplot_", Region, "_Figure5_mat_TOE_95_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(a)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Marg.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(b)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),main_="Dep.", coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(c)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_margdep_95, p_i2, p_i1, NaN, ylim_=c(-0.01, 0.2), main_="Marg.-dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(d)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_marg_95, p_i2, p_i1,  NaN, ylim_=c(-0.01, 0.2),main_="Marg.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(e)",zlim_=c(1980,2086))
fastplot_Indiv(list_array_Full_dep_95, p_i2, p_i1,  NaN,ylim_=c(-0.01, 0.2),main_="Dep.",coord_baseline,coord_sub_ts,plot_toe_mat=TRUE,label_68_95="95",subplot_name="(f)",zlim_=c(1980,2086))
dev.off()


### Mat Contribution Indv Full
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

pdf(paste0("fastplot_", Region, "_Figure6_Contrib_mat_IndivFull_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Contrib_matrix(apply(array_Med_C_marg_Indiv,c(1,2),median),
                        apply(array_Med_C_dep_Indiv,c(1,2),median),
                        apply(array_Med_C_int_Indiv,c(1,2),median),
                        plot_matrix=TRUE, subplot_name=c("(a)", "(b)", "(c)"))
fastplot_Contrib_matrix(Med_C_marg_Full,
                        Med_C_dep_Full,
                        Med_C_int_Full,
                        plot_matrix=TRUE,
                        subplot_name=c("(d)", "(e)", "(f)"))
dev.off()

source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
p_i2=16
p_i1=16
pdf(paste0("fastplot_", Region, "_Figure7_Contrib_ts_IndivFull_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(2,3,3,3.5),mar=c(5,4.5,4,1) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Contrib_Indiv(res_Indiv_Contrib,p_i2,p_i1,coord_sub_ts)
fastplot_Contrib(res_Full_Contrib[["CDFt_pooled_ModRefCNRMCM6"]],p_i2,p_i1,plot_ts=TRUE, plot_matrix=FALSE, plot_barplot=FALSE,coord_sub_ts)
dev.off()


source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")
p_i2=16
p_i1=16
pdf(paste0("fastplot_", Region, "_Figure8_68_Sum_IQR_Indiv_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted,  ylim_=c(-0.01, 0.2),
               main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68", subplot_name="(a)")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68", subplot_name="(b)")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="68", subplot_name="(c)")
fastplot_Indiv(list_array_Indiv_margdep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68", subplot_name="(d)")
fastplot_Indiv(list_array_Indiv_marg_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68", subplot_name="(e)")
fastplot_Indiv(list_array_Indiv_dep_68,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2), 
               main_="Dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="68", subplot_name="(f)")
dev.off()

pdf(paste0("fastplot_", Region, "_Figure8_95_Sum_IQR_Indiv_",period_baseline ,".pdf"),width=15, height=10)
par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,7,4,2) + 0.1)
#par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95", subplot_name="(a)")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95", subplot_name="(b)")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2), 
               main_="Dep.", coord_baseline,coord_sub_ts,plot_sum=TRUE,label_68_95="95", subplot_name="(c)")
fastplot_Indiv(list_array_Indiv_margdep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.-dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95", subplot_name="(d)")
fastplot_Indiv(list_array_Indiv_marg_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Marg.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95", subplot_name="(e)")
fastplot_Indiv(list_array_Indiv_dep_95,p_i2, p_i1, array_dotted, ylim_=c(-0.01, 0.2),
               main_="Dep.", coord_baseline,coord_sub_ts,plot_iqr=TRUE,label_68_95="95", subplot_name="(f)")
dev.off()






### Scatterplot ToE/Periode de Ref
rm(list=ls())
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_TCE.R")
source("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Code/function_Plot_Figures_Paper.R")

rdylgnPal<-colorRampPalette(brewer.pal(11, "RdYlGn"))(500)
brbgPal <- rev(colorRampPalette(brewer.pal(11, "BrBG"))(500))
rbPal<-rev(colorRampPalette(brewer.pal(11, "RdBu"))(500))
SpecPaldisc13<-c(rbPal[1], rev(brewer.pal(n=11, name="Spectral")), rbPal[500])
rdylbuPal<-rev(c(rev(colorRampPalette(brewer.pal(9, "YlGnBu"))(250)),colorRampPalette(brewer.pal(9, "YlOrRd"))(250)))
col_main=c("dodgerblue","darkorange", "gray", "indianred1") # "chartreuse3"
magma13<-rev(magma(14))[1:13]

Ref_SlidingWindow="1871_1900"
Region="Bretagne"
x_axis_lab="Wind speed [m/s]"
y_axis_lab="Precipitation [mm/d]"

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref",Ref_SlidingWindow))
load(paste0("Res_array_prob_list30_Ref",Ref_SlidingWindow,"_sfcWindmax_pr_CNRMCM6_winter_1850_2100_Bretagne.RData"))
setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Data/Results_Data/Bretagne/sfcWindmax_pr/SlidingWindow30_Ref",Ref_SlidingWindow,"/marg_uncert/"))
load(paste0("Res_array_v2bis_marg_uncert_B100_prob_list30_Ref",Ref_SlidingWindow,"_sfcWindmax_pr_CNRMCM6_winter_1850_2100_Bretagne.RData"))

coord_sub_ts=c(which(Ref_SlidingWindow==Label_SlidingWindowX_CNRMCM6):222)
coord_sub_ts

setwd(paste0("/home/starmip/bfran/LSCE_These/Compound_Event/TCE6/Plots/Plot_Figures"))
coord_baseline=22
nb_listX=201
xlab_name=seq.int(1886,2086,20)
xlab_pos=seq(1,nb_listX,20)
p_i2=16
p_i1=16

tmp_toe_scatt_margdep_68=tmp_toe_scatt_marg_68=tmp_toe_scatt_dep_68=c()
tmp_toe_scatt_margdep_95=tmp_toe_scatt_marg_95=tmp_toe_scatt_dep_95=c()

for(coord_baseline in 22:222){
  print(coord_baseline)
  tmp_toe_scatt_margdep_68=c(tmp_toe_scatt_margdep_68,compute_ToE_matrix_from_array(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  tmp_toe_scatt_marg_68=c(tmp_toe_scatt_marg_68,compute_ToE_matrix_from_array(array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  tmp_toe_scatt_dep_68=c(tmp_toe_scatt_dep_68,compute_ToE_matrix_from_array(array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  
  tmp_toe_scatt_margdep_95=c(tmp_toe_scatt_margdep_95,compute_ToE_matrix_from_array(array_prob_margdep_v2bis_95_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  tmp_toe_scatt_marg_95=c(tmp_toe_scatt_marg_95,compute_ToE_matrix_from_array(array_prob_marg_v2bis_95_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  tmp_toe_scatt_dep_95=c(tmp_toe_scatt_dep_95,compute_ToE_matrix_from_array(array_prob_dep_v2bis_95_Joe_marg_uncert_CNRMCM6,30, "1850_2100", coord_baseline)$TOE[p_i2,p_i1])
  
}

pdf(paste0("fastplot_", Region, "_Figure9_68_Scatterplot_Toe_BaselinePeriod.pdf"),width=5, height=5)
plot(tmp_toe_scatt_margdep_68, xaxt="n",ylab="ToE",type='l',xlab="Baseline period")
lines(tmp_toe_scatt_marg_68, xaxt="n",ylab="ToE",col="dodgerblue")
lines(tmp_toe_scatt_dep_68, xaxt="n",ylab="ToE",col="darkorange")
axis(1,xlab_pos,xlab_name)
legend("bottomright", c("Marg.-Dep.","Marg.", "Dep."), lty=c(1,1,1), col=c("black","dodgerblue","darkorange"))
dev.off()

pdf(paste0("fastplot_", Region, "_Figure9_95_Scatterplot_Toe_BaselinePeriod.pdf"),width=5, height=5)
plot(tmp_toe_scatt_margdep_95, xaxt="n",ylab="ToE",type='l',xlab="Baseline period")
lines(tmp_toe_scatt_marg_95, xaxt="n",ylab="ToE",col="dodgerblue")
lines(tmp_toe_scatt_dep_95, xaxt="n",ylab="ToE",col="darkorange")
axis(1,xlab_pos,xlab_name)
legend("bottomright", c("Marg.-Dep.","Marg.", "Dep."), lty=c(1,1,1), col=c("black","dodgerblue","darkorange"))
dev.off()


#### Figure 1
# pdf(paste0("fastplot_", Region, "_Figure1_",period_baseline ,".pdf"), width=15, height=10)#,width=1333, height=666)
# par(mfrow=c(2,3),oma=c(3,3,3,3),mar=c(5,4,4,5) + 0.1)
## Proba
fastplot_time_series(array_prob_margdep_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_margdep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_margdep_68_Joe_CNRMCM6, array_prob_margdep_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Marg.-dep.",coord_baseline,coord_sub_ts, subplot_name="(a)")
fastplot_time_series(array_prob_marg_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_marg_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_marg_68_Joe_CNRMCM6,array_prob_marg_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Marg.",coord_baseline,coord_sub_ts)
fastplot_time_series(array_prob_dep_v2bis_68_Joe_marg_uncert_CNRMCM6, array_prob_dep_v2bis_95_Joe_marg_uncert_CNRMCM6,
                     p_i2, p_i1, array_prob_dep_68_Joe_CNRMCM6,array_prob_dep_95_Joe_CNRMCM6, ylim_=c(-0.01, 0.2),main="Dep.",coord_baseline,coord_sub_ts)
dev.off()




# fastplot_Boxplot_ToE<-function(list_array_prob_Indiv_margdep, 
#                                list_array_prob_Indiv_marg,
#                                list_array_prob_Indiv_dep,
#                                list_array_prob_Full_margdep,
#                                list_array_prob_Full_marg,
#                                list_array_prob_Full_dep,
#                                p_i2, p_i1, array_dotted=NaN,
#                                main_=NaN, coord_baseline=NaN,coord_sub_ts=NaN,label_68_95="NaN",
#                                is_TardiveFrost=FALSE){
#   nb_Mod_final_IndivEnsemble=length(names(get(paste0("list_array_prob_Indiv_margdep"))))
#   for(v in c("margdep", "marg", "dep")){
#     count_ToE=0
#     indiv_ToE=c()
#     for(i_Mod in names(get(paste0("list_array_prob_Indiv_", v)))){
#       print(v)
#       tmp_proba=get(paste0("list_array_prob_Indiv_", v))[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
#       tmp_toe=compute_ToE_matrix_from_array(get(paste0("list_array_prob_Indiv_", v))[[i_Mod]],30, "1850_2100", coord_baseline)
#       indiv_ToE=c(indiv_ToE,tmp_toe$TOE[p_i2,p_i1])
#       if(!is.na(tmp_toe$coord_TOE[p_i2,p_i1])){
#         count_ToE=count_ToE+1
#       }
#     }
#     print(count_ToE)
#     assign(paste0("tmp_TOE_Indiv_", v), indiv_ToE)
#   }
#   df <- data.frame( id = c(rep("Indiv",nb_Mod_final_IndivEnsemble))
#                     ,
#                     #        rep("Full", nb_Mod_final_FullEnsemble)),
#                     Margdep = c(tmp_TOE_Indiv_margdep),
#                     Marg = c(tmp_TOE_Indiv_marg),
#                     Dep = c(tmp_TOE_Indiv_dep))
#   proportion <- c(sum(!is.na(tmp_TOE_Indiv_margdep)),
#                   sum(!is.na(tmp_TOE_Indiv_marg)),
#                   sum(!is.na(tmp_TOE_Indiv_dep)))/(13)
#   boxplot(df[,-1],boxwex=proportion, boxfill="coral",  names=c("Marg. dep.","Marg.", "Dep."), ylim=c(1980,2086)) 
#   ### Full
#   for(v in c("margdep", "marg", "dep")){
#     full_ToE=c()
#     for(i_Mod in names(get(paste0("list_array_prob_Full_", v)))){
#       print(v)
#       tmp_proba=get(paste0("list_array_prob_Full_", v))[[i_Mod]][p_i2,p_i1,coord_sub_ts,1]
#       tmp_toe=compute_ToE_matrix_from_array(get(paste0("list_array_prob_Full_", v))[[i_Mod]],30, "1850_2100", coord_baseline)
#       full_ToE=c(full_ToE,tmp_toe$TOE[p_i2,p_i1])
#     }
#     assign(paste0("tmp_TOE_Full_", v), full_ToE)
#     print(full_ToE)
#   }
#   points(1.25,tmp_TOE_Full_margdep, pch=8, col="royalblue",
#          lwd=3)
#   points(2.25,tmp_TOE_Full_marg, pch=8, col="royalblue", lwd=3)
#   points(3.25,tmp_TOE_Full_dep, pch=8, col="royalblue", lwd=3)
#   legend("bottomright",
#          legend = c("Indiv-Ensemble", "Full-Ensemble"),
#          col = c("coral", "royalblue"),
#          pch=c(15, 8))
# }








