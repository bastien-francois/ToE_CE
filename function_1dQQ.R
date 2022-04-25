
#### ne pas toucher #### 
# QQb<-function(Rc,Mc,Mp){
#   N=ncol(Mc)
#   I_Mc=nrow(Mc)
#   I_Mp=nrow(Mp)
#   if(is.null(N)&is.null(I_Mc)){
#     N=1
#     I_Mc=length(Mc)
#     I_Mp=length(Mp)
#     Rc=matrix(Rc,ncol=1)
#     Mc=matrix(Mc,ncol=1)
#     Mp=matrix(Mp,ncol=1)
#   }
#   Mch=matrix(rep(NA,I_Mc*N),ncol=N)
#   Mph=matrix(rep(NA,I_Mp*N),ncol=N)
#   for(k in 1:N){ #for each column (variable)
#     FMc=ecdf(Mc[,k])
#     Mch[,k]=quantile(Rc[,k],probs=FMc(Mc[,k]),type=7)
#     # Mph[,k]=quantile(Rc[,k],probs=FMc(Mp[,k]),type=7) not good for proj...
#     ### approx probs for Mph to avoid changing rank structure for projection period.
#     FMC=FMc(Mc[,k])
#     probs_FMc_Mp=approx(Mc[,k], FMC, Mp[,k], yleft = min(FMC), yright = max(FMC), 
#                         ties = "mean")$y
#     Mph[,k]=quantile(Rc[,k],probs=probs_FMc_Mp,type=7)
#   }
#   return(list(Mch=Mch,Mph=Mph))
# }

#### end ne pas toucher #### 

#Quantile-Quantile version that allow to correct Mp in CC context (out-of-sample values) (Gudmundsson et al., 2012)
# QQb_new<-function(Rc,Mc,Mp,p=0){ #the parameter p is necessary to implement QQ for MRec (see Pegram and Bardossy, 2012): p=0 classical QQ
#   N=ncol(Mc)
#   I_Mc=nrow(Mc)
#   I_Mp=nrow(Mp)
#   if(is.null(N)&is.null(I_Mc)){ #for 1d vectors
#     N=1
#     I_Mc=length(Mc)
#     I_Mp=length(Mp)
#     Rc=matrix(Rc,ncol=1)
#     Mc=matrix(Mc,ncol=1)
#     Mp=matrix(Mp,ncol=1)
#   }
#   Mch=matrix(rep(NA,I_Mc*N),ncol=N)
#   Mph=matrix(rep(NA,I_Mp*N),ncol=N)
#   for(k in 1:N){ #for each column (variable)
#     #Classic quantile-quantile for Mc
#     FMc=ecdf(Mc[,k])
#     Mch[,k]=quantile(Rc[,k],probs=FMc(Mc[,k])*(1-p)+p,type=4)
#     #Save the correction done for highest and lowest quantiles (will be used later to correct Mp in a context of climate change)
#     correc_high_quntl=max(Mc[,k])-max(Mch[,k])
#     correc_low_quntl=min(Mc[,k])-min(Mch[,k])
#     
#     #Quantile-quantile for Mp
#     # which value in Mp are within [min(Mc),max(Mc)]?
#     in_range=((Mp[,k]<=max(Mc[,k]))&(Mp[,k]>=min(Mc[,k])))
#     #for these values, classic quantile quantile
#     Mph[in_range,k]=quantile(Rc[,k],probs=FMc(Mp[in_range,k])*(1-p)+p,type=4)
#     #for out-of-sample values of Mp, same correction than for Mc
#     out_range_low=(Mp[,k]<min(Mc[,k]))
#     out_range_high=(Mp[,k]>max(Mc[,k]))
#     Mph[out_range_low,k]<-Mp[out_range_low,k]-correc_low_quntl
#     Mph[out_range_high,k]<-Mp[out_range_high,k]-correc_high_quntl
#   }
#   return(list(Mch=Mch,Mph=Mph))
# }

#Quantile-Quantile version that allow to correct Mp in CC context (out-of-sample values) (Gudmundsson et al., 2012)
#as in Francois2021
QQb_new<-function(Rc,Mc,Mp){ 
  N=ncol(Mc)
  I_Mc=nrow(Mc)
  I_Mp=nrow(Mp)
  if(is.null(N)&is.null(I_Mc)){ #for 1d vectors
    N=1
    I_Mc=length(Mc)
    I_Mp=length(Mp)
    Rc=matrix(Rc,ncol=1)
    Mc=matrix(Mc,ncol=1)
    Mp=matrix(Mp,ncol=1)
  }
  Mch=matrix(rep(NA,I_Mc*N),ncol=N)
  Mph=matrix(rep(NA,I_Mp*N),ncol=N)
  for(k in 1:N){ #for each column (variable)
    #Classic quantile-quantile for Mc
    FMc=ecdf(Mc[,k])
    FMC=FMc(Mc[,k])
    FRc=ecdf(Rc[,k])
    FRC=FRc(Rc[,k])
    # Mch[,k]=quantile(Rc[,k],probs=FMc(Mc[,k]),type=7)
    Mch[,k]=approx(FRC,Rc[,k], FMC, yleft=min(Rc[,k]), yright= max(Rc[,k]))$y
    #Save the correction done for highest and lowest quantiles (will be used later to correct Mp in a context of climate change)
    correc_high_quntl=max(Mc[,k])-max(Mch[,k])
    correc_low_quntl=min(Mc[,k])-min(Mch[,k])
    
    #Quantile-quantile for Mp
    # which value in Mp are within [min(Mc),max(Mc)]?
    in_range=((Mp[,k]<=max(Mc[,k]))&(Mp[,k]>=min(Mc[,k])))
    #for these values, classic quantile quantile with approximation for interpolation
    # Mph[in_range,k]=quantile(Rc[,k],probs=FMc(Mp[in_range,k])*(1-p)+p,type=4)
    
    probs_FMc_Mp=approx(Mc[,k], FMC, Mp[in_range,k], yleft = min(FMC), yright = max(FMC), 
                        ties = "mean")$y
    # Mph[in_range,k]=quantile(Rc[,k],probs=probs_FMc_Mp,type=7)
    Mph[in_range,k]=approx(FRC,Rc[,k],probs_FMc_Mp,  yleft=min(Rc[,k]), yright= max(Rc[,k]))$y
    
    #for out-of-sample values of Mp, same correction than for Mc (Boe, Deque 2007)
    out_range_low=(Mp[,k]<min(Mc[,k]))
    out_range_high=(Mp[,k]>max(Mc[,k]))
    Mph[out_range_low,k]<-Mp[out_range_low,k]-correc_low_quntl
    Mph[out_range_high,k]<-Mp[out_range_high,k]-correc_high_quntl
  }
  return(list(Mch=Mch,Mph=Mph))
}
QQb_new_for_MRec<-function(Rc,Mc,Mp,p=0){ #the parameter p is necessary to implement QQ for MRec (see Pegram and Bardossy, 2012): p=0 classical QQ
  N=ncol(Mc)
  I_Mc=nrow(Mc)
  I_Mp=nrow(Mp)
  if(is.null(N)&is.null(I_Mc)){ #for 1d vectors
    N=1
    I_Mc=length(Mc)
    I_Mp=length(Mp)
    Rc=matrix(Rc,ncol=1)
    Mc=matrix(Mc,ncol=1)
    Mp=matrix(Mp,ncol=1)
  }
  Mch=matrix(rep(NA,I_Mc*N),ncol=N)
  Mph=matrix(rep(NA,I_Mp*N),ncol=N)
  for(k in 1:N){ #for each column (variable)
    #Classic quantile-quantile for Mc
    FMc=ecdf(Mc[,k])
    FMC=FMc(Mc[,k])
    FRc=ecdf(Rc[,k])
    FRC=FRc(Rc[,k])
    # Mch[,k]=quantile(Rc[,k],probs=FMc(Mc[,k])*(1-p)+p,type=3)
    Mch[,k]=approx(FRC,Rc[,k], FMC*(1-p)+p, yleft=min(Rc[,k]), yright= max(Rc[,k]))$y
    #Save the correction done for highest and lowest quantiles (will be used later to correct Mp in a context of climate change)
    correc_high_quntl=max(Mc[,k])-max(Mch[,k])
    correc_low_quntl=min(Mc[,k])-min(Mch[,k])
    
    #Quantile-quantile for Mp
    # which value in Mp are within [min(Mc),max(Mc)]?
    in_range=((Mp[,k]<=max(Mc[,k]))&(Mp[,k]>=min(Mc[,k])))
    #for these values, classic quantile quantile
    probs_FMc_Mp=approx(Mc[,k], FMC, Mp[in_range,k], yleft = min(FMC), yright = max(FMC), 
                        ties = "mean")$y
    Mph[in_range,k]=approx(FRC,Rc[,k],probs_FMc_Mp*(1-p)+p,  yleft=min(Rc[,k]), yright= max(Rc[,k]))$y
    # Mph[in_range,k]=quantile(Rc[,k],probs=FMc(Mp[in_range,k])*(1-p)+p,type=3)
    #for out-of-sample values of Mp, same correction than for Mc
    out_range_low=(Mp[,k]<min(Mc[,k]))
    out_range_high=(Mp[,k]>max(Mc[,k]))
    Mph[out_range_low,k]<-Mp[out_range_low,k]-correc_low_quntl
    Mph[out_range_high,k]<-Mp[out_range_high,k]-correc_high_quntl
  }
  return(list(Mch=Mch,Mph=Mph))
} #used in Francois2020

### reorder data_marg with struct. of data_struct data of dim. nb_point x nb_var
reorder_data<-function(data_marg,data_struct){ 
  nb_dim=ncol(data_struct)
  nb_point=nrow(data_struct)
  res = matrix( nrow = nb_point , ncol = nb_dim)
  for(k in 1:nb_dim){
    sorted_data_marg=sort(data_marg[,k])
    idx=rank(data_struct[,k],ties.method="first")
    res[,k] = sorted_data_marg[idx]
  }
  return(res)
}


#### For IPSLbis with MRec
gaussianization_MRec_for_IPSLbis<-function(mat_Mc, is.ratio){
  nb_dim=ncol(mat_Mc)
  nb_point_Mc=nrow(mat_Mc)
  gauss_mat_Mc = matrix(NaN, nrow = nb_point_Mc , ncol = nb_dim)
  for(j in 1:nb_dim){
    if(is.ratio[j]==TRUE){ #if is precip
      #Pj0 is the proportion of dry event for precip in Mc
      Pj0=sum(mat_Mc[,j]==0)/nb_point_Mc
      if(Pj0>0){ #(TRUE for most of the cases)
        mat_Mc_Positif_bool=mat_Mc[,j]!=0
        #For 0 values, projection to the inverse of the normal distrib with proba Pj0/2 (see Bardossy and Pegram, 2012)
        gauss_mat_Mc[!mat_Mc_Positif_bool,j]<-rep(qnorm(Pj0/2), sum(!mat_Mc_Positif_bool))
        #For values >0, QQ applied with parameter p=Pj0 (see Bardossy and Pegram, 2012)
        FMc_Positif=ecdf(mat_Mc[mat_Mc_Positif_bool,j])
        FMC_Positif = FMc_Positif(mat_Mc[mat_Mc_Positif_bool,j])
        ### For probs=1, on a qnorm=inf. On approche avec eps
        epsilon=1/length(FMC_Positif) #valeur d'écart entre deux proba
        FMC_Positif[which(FMC_Positif==1)]<-1-epsilon/2
        # qq_correc=QQb_new(sample_gaussienQQ,mat_Mc[mat_Mc_Positif_bool,j],mat_Mp[mat_Mp_Positif_bool,j],p=Pj0)
        gauss_mat_Mc[mat_Mc_Positif_bool,j]<-qnorm(FMC_Positif*(1-Pj0)+Pj0)
      }else{ #rare case to treat problem when Pj0=0, i.e. no dry event at all in Mc 
        #(je ne sais plus bien pourquoi j avais fait ce cas, la correction me semble maintenant equivalente a celle ci-dessus)
        #je laisse pour te laisser mon code original mais c'est peut etre a enlever 
        FMc=ecdf(mat_Mc[,j])
        FMC = FMc(mat_Mc[,j])
        ### For probs=1, on a qnorm=inf. On approche avec eps
        epsilon=1/length(FMC) #valeur d'écart entre deux proba
        FMC[which(FMC==1)]<-1-epsilon/2
        gauss_mat_Mc[,j]<-qnorm(FMC)
      }
    }
    if(is.ratio[j]==FALSE){ #if is not precip
      FMc=ecdf(mat_Mc[,j])
      FMC = FMc(mat_Mc[,j])
      epsilon=1/length(FMC) #valeur d'écart entre deux proba
      FMC[which(FMC==1)]<-1-epsilon/2
      gauss_mat_Mc[,j]<-qnorm(FMC)
    }
  }
  return(list(gauss_mat_Mc=gauss_mat_Mc))
}

MRec_for_IPSLbis<-function(Rc,Rp,Mc,Mp,Rp_from_CDFt,ratio.seq){ #Mrec 
  #Initialization of the results
  nb_dim=ncol(Rc)
  nb_point_Mc=nrow(Mc)
  nb_point_Mp=nrow(Mp)
  nb_point_Rc=nrow(Rc)
  nb_point_Rp=nrow(Rp)
  Rch= matrix( nrow = nb_point_Rc , ncol = nb_dim)
  Rph= matrix( nrow = nb_point_Rp , ncol = nb_dim)
  
  #Save RC for calib.
  initial_Rc=Rc
  
  #Following the steps from Pegram and Bardossy, 2012
  #Gaussianization of Rc
  GRc=gaussianization_MRec_for_IPSLbis(Rc,is.ratio=ratio.seq)$gauss_mat_Mc
  GRp=gaussianization_MRec_for_IPSLbis(Rp,is.ratio=ratio.seq)$gauss_mat_Mc
  
  #Spearman Cor. and SVD for W
  C0=cor(GRc,method="pearson")
  C1=cor(GRp,method="pearson")
  
  #Gaussianization of the model Mc and Mp
  GMc=gaussianization_MRec_for_IPSLbis(Mc,is.ratio=ratio.seq)$gauss_mat_Mc
  GMp=gaussianization_MRec_for_IPSLbis(Mp,is.ratio=ratio.seq)$gauss_mat_Mc
  
  #For Mc first
  R0=cor(GMc,method="pearson")
  R1=cor(GMp,method="pearson")
  
  #Formula from Pegram, 2012
  nom=(1+R1)/(1+R0)*(1+C0)-(1-R1)/(1-R0)*(1-C0)
  denom=(1+R1)/(1+R0)*(1+C0)+(1-R1)/(1-R0)*(1-C0)
  C1_new=nom/denom
  diag(C1_new)=1
  
  print(max(C1_new))
  ### Learn to decorrelate GRp
  SVD_C1=svd(C1)
  A2=SVD_C1$u
  B2=SVD_C1$v
  D2=diag(nb_dim)
  diag(D2)<-SVD_C1$d
  sqrtD2=sqrtm(D2)
  sqrtD2m1=solve(sqrtD2)
  T_=A2%*%sqrtD2m1%*%t(A2)
  
  #### Learn to recorrelate to C1_new
  SVD_C=svd(C1_new)
  A1=SVD_C$u
  B1=SVD_C$v
  D1=diag(nb_dim)
  diag(D1)<-SVD_C$d
  S_=A1%*%sqrtm(D1)%*%t(A1)
  
  # Decorr GRp and recorr to C1_new
  F_=T_%*%S_
  Recorr_GRp=GRp%*%F_
  
  #QQ to correct marginals of V1
  tmpRph= matrix( nrow = nb_point_Rp , ncol = nb_dim)
  tmpRch = initial_Rc
  #### put marginal from CDFt output (IPSLMRbili_SAFRANdetbili)
  tmpRph=reorder_data(Rp_from_CDFt,Recorr_GRp) 
  
  #### Feeding the final results matrices with corrected variables with MRec
  Rch<-tmpRch
  Rph<-tmpRph
  
  return(list(Rch=Rch,Rph=Rph)) #GRc = GRc, GRp=GRp, Recorr_GRp = Recorr_GRp, GMc = GMc, GMp=GMp, GRecorr_GRp=GRecorr_GRp
}

### end For IPSLbis with MRec



#### CDFt_SSR as in Francois2020
CDFt_SSR <- function(ObsRp, DataGp, DataGf, th_O=NaN, th_M=NaN, npas = 1000, dev = 2){ #Same as Vrac et al. with npas=1000 by default
  if(is.na(th_O) | is.na(th_M)){
    if(length(ObsRp[which(ObsRp>0)])==0){
      th_O=max(DataGp,DataGf)
    }
    else{
      th_O = min(ObsRp[which(ObsRp>0)],na.rm=TRUE)
    }
    th_M = min(DataGp[which(DataGp>0)], DataGf[which(DataGf>0)], na.rm=TRUE)
    #cat("th_O=",th_O,"th_M=",th_M,"\n")
  }
  else{
    th_O = th_O
    th_M = th_M
  }
  # print(th_O)
  # print(th_M)
  ### st for stoch simulations : from 0 to Unif [0,th]
  ObsRp_st = ObsRp
  DataGp_st = DataGp
  DataGf_st = DataGf
  WObs = which(ObsRp<=th_O)
  ObsRp_st[WObs] = runif(length(WObs),0,th_O)
  WGp = which(DataGp<=th_M)
  DataGp_st[WGp] = runif(length(WGp),0,th_M)
  WGf = which(DataGf<=th_M)
  DataGf_st[WGf] = runif(length(WGf),0,th_M)
  ### 
  
  
  ### Normalization based on the 90th quantile
  Q90O = quantile(ObsRp_st, probs=0.9, na.rm=TRUE)
  Q90Gp = quantile(DataGp_st, probs=0.9, na.rm=TRUE)
  DataGp2 = DataGp_st * (Q90O/Q90Gp)
  DataGf2 = DataGf_st * (Q90O/Q90Gp)
  
  FRp=ecdf(ObsRp_st)
  FGp=ecdf(DataGp2)
  FGf=ecdf(DataGf2)
  
  a=abs(max(DataGf_st, na.rm=TRUE)-max(DataGp_st, na.rm=TRUE))
  m=0
  M=max(ObsRp_st, DataGp_st, DataGf_st, na.rm=TRUE)+dev*a
  x=seq(m,M,length.out=npas)
  
  FGF=FGf(x)
  FGP=FGp(x)
  FRP=FRp(x)
  
  FGPm1.FGF=quantile(DataGp2,probs=FGF, na.rm=TRUE)
  
  FRF=FRp(FGPm1.FGF)
  
  ######################################
  # FRf=FRp with shift for x<min(DataGf)
  
  if(min(ObsRp_st, na.rm=TRUE)<min(DataGf2, na.rm=TRUE)){
    # cat("FRf=FRp with shift for x < max(DataGf)\n")
    i=1
    while(x[i]<=quantile(ObsRp_st,probs=FRF[1],na.rm=TRUE)){
      i=i+1
    }
    
    j=1
    while(x[j]<min(DataGf2, na.rm=TRUE)){
      j=j+1
    }
    
    k=i
    while(j>0 && k>0){
      FRF[j]=FRP[k]
      j=j-1
      k=k-1
    }
    
    ##########    
    
    if(j>0){
      for(k in j:1){
        FRF[k]=0
      }
    }
    
  }
  
  ######################################
  # FRf=FRp with shift for x>max(DataGf)
  
  
  if(FRF[length(x)]<1){
    #cat("FRf=FRp with shift for x > max(DataGf)\n")
    i=length(x)
    QQ=quantile(ObsRp_st,probs=FRF[length(x)], na.rm=TRUE)
    while(x[i]>=QQ){
      i=i-1
    }
    i=i+1
    
    j=length(x)-1
    while(j>0 && FRF[j]==FRF[length(x)]){
      j=j-1
    }
    
    if(j==0){
      stop("In CDFt, dev must be higher\n")
    }
    
    dif=min((length(x)-j),(length(x)-i))
    FRF[j:(j+dif)]=FRP[i:(i+dif)]
    k=j+dif
    
    if(k<length(x)){
      FRF[k:(length(x))]=1
    }
    
  }
  
  
  ######################################################################################
  ### Quantile-matching based on the new large-scale CDF and downscaled local-scale CDF.
  
  ############
  NaNs.indices = which(is.na(DataGf2))
  No.NaNs.indices = which(!is.na(DataGf2))
  
  qntl = array(NaN, dim=length(DataGf2))
  qntl[No.NaNs.indices] = FGf(DataGf2[No.NaNs.indices])
  
  xx = array(NaN, dim=length(DataGf2))
  xx = approx(FRF,x,qntl,yleft=x[1],yright=x[length(x)],ties='mean')
  
  ##############################################
  
  
  
  #################
  #################
  
  # CALCULER LE FRf(0)
  # POUR TOUS i TQ FRf(xx$y[i])<=FRf(0), THEN xx$y[i]=0
  # COMME C'EST DU QQ (i.e., FRf(xx$y[i])==FGf(DataGf2[i])), POUR TOUS i TQ FGf(DataGf2[i]))<=FRf(0), THEN xx$y[i]=0
  #                                                 <=>                              qntl[i]<=FRf(0), THEN xx$y[i]=0
  
  
  EmpGf2_th = (ecdf(DataGf2))(th_O)
  EmpGp2m1.Gf2_th = quantile(DataGp2, probs=EmpGf2_th, na.rm=TRUE)
  FRf_th = EmpFp.Gp2m1.Gf2_th = (ecdf(ObsRp_st))(EmpGp2m1.Gf2_th) ##### METTRE ObsRp_st pour avoir la proba qui evolue
  #cat("FRf_th =",FRf_th,"\n")
  xx$y[which(qntl<=(FRf_th))] = 0 ### A METTRE ABSOLUMENT
  xx$y[which(xx$y<=th_O)] = 0 ### NE CHANGE RIEN MAIS A METTRE PAR SECURITE
  
  #cat("length(which(xx$y<=th)) =",length(which(xx$y<=th)),"\n")
  #################
  #################
  
  
  
  #######################################################################################
  FGp=ecdf(DataGp)
  FGf=ecdf(DataGf)
  FGP=FGp(x)
  FGF=FGf(x)
  
  return(list(x=x,FRp=FRP,FGp=FGP,FGf=FGF,FRf=FRF,FRf_0=FRf_th,DS=xx$y))
  
}

CDFt_SSRbis <- function(ObsRp, DataGp, DataGf, th_O=NaN, th_M=NaN, npas = 1000, dev = 2){ #Same as Vrac et al. with npas=1000 by default 
  #No normalization to avoid cases when (Q90O/Q90Gp) is really high see below
  if(is.na(th_O) | is.na(th_M)){
    if(length(ObsRp[which(ObsRp>0)])==0){
      th_O=max(DataGp,DataGf)
    }
    else{
      th_O = min(ObsRp[which(ObsRp>0)],na.rm=TRUE)
    }
    th_M = min(DataGp[which(DataGp>0)], DataGf[which(DataGf>0)], na.rm=TRUE)
    cat("th_O=",th_O,"th_M=",th_M,"\n")
  }
  else{
    th_O = th_O
    th_M = th_M
  }
  # print(th_O)
  # print(th_M)
  ### st for stoch simulations : from 0 to Unif [0,th]
  ObsRp_st = ObsRp
  DataGp_st = DataGp
  DataGf_st = DataGf
  WObs = which(ObsRp<=th_O)
  ObsRp_st[WObs] = runif(length(WObs),0,th_O)
  WGp = which(DataGp<=th_M)
  DataGp_st[WGp] = runif(length(WGp),0,th_M)
  WGf = which(DataGf<=th_M)
  DataGf_st[WGf] = runif(length(WGf),0,th_M)
  ### 
  
  
  ### Normalization based on the 90th quantile
  # Q90O = quantile(ObsRp_st, probs=0.9, na.rm=TRUE)
  # Q90Gp = quantile(DataGp_st, probs=0.9, na.rm=TRUE)
  # DataGp2 = DataGp_st * (Q90O/Q90Gp)
  # DataGf2 = DataGf_st * (Q90O/Q90Gp)
  
  # No normalization to avoid cases when (Q90O/Q90Gp) is really high
  DataGp2=DataGp_st
  DataGf2=DataGf_st
  
  FRp=ecdf(ObsRp_st)
  FGp=ecdf(DataGp2)
  FGf=ecdf(DataGf2)
  
  a=abs(max(DataGf_st, na.rm=TRUE)-max(DataGp_st, na.rm=TRUE))
  m=0
  M=max(ObsRp_st, DataGp_st, DataGf_st, na.rm=TRUE)+dev*a
  x=seq(m,M,length.out=npas)
  
  FGF=FGf(x)
  FGP=FGp(x)
  FRP=FRp(x)
  
  FGPm1.FGF=quantile(DataGp2,probs=FGF, na.rm=TRUE)
  
  FRF=FRp(FGPm1.FGF)
  
  ######################################
  # FRf=FRp with shift for x<min(DataGf)
  
  if(min(ObsRp_st, na.rm=TRUE)<min(DataGf2, na.rm=TRUE)){
    # cat("FRf=FRp with shift for x < max(DataGf)\n")
    i=1
    while(x[i]<=quantile(ObsRp_st,probs=FRF[1],na.rm=TRUE)){
      i=i+1
    }
    
    j=1
    while(x[j]<min(DataGf2, na.rm=TRUE)){
      j=j+1
    }
    
    k=i
    while(j>0 && k>0){
      FRF[j]=FRP[k]
      j=j-1
      k=k-1
    }
    
    ##########    
    
    if(j>0){
      for(k in j:1){
        FRF[k]=0
      }
    }
    
  }
  
  ######################################
  # FRf=FRp with shift for x>max(DataGf)
  
  
  if(FRF[length(x)]<1){
    #cat("FRf=FRp with shift for x > max(DataGf)\n")
    i=length(x)
    QQ=quantile(ObsRp_st,probs=FRF[length(x)], na.rm=TRUE)
    while(x[i]>=QQ){
      i=i-1
    }
    i=i+1
    
    j=length(x)-1
    while(j>0 && FRF[j]==FRF[length(x)]){
      j=j-1
    }
    
    if(j==0){
      stop("In CDFt, dev must be higher\n")
    }
    
    dif=min((length(x)-j),(length(x)-i))
    FRF[j:(j+dif)]=FRP[i:(i+dif)]
    k=j+dif
    
    if(k<length(x)){
      FRF[k:(length(x))]=1
    }
    
  }
  
  
  ######################################################################################
  ### Quantile-matching based on the new large-scale CDF and downscaled local-scale CDF.
  
  ############
  NaNs.indices = which(is.na(DataGf2))
  No.NaNs.indices = which(!is.na(DataGf2))
  
  qntl = array(NaN, dim=length(DataGf2))
  qntl[No.NaNs.indices] = FGf(DataGf2[No.NaNs.indices])
  
  xx = array(NaN, dim=length(DataGf2))
  xx = approx(FRF,x,qntl,yleft=x[1],yright=x[length(x)],ties='mean')
  
  ##############################################
  
  
  
  #################
  #################
  
  # CALCULER LE FRf(0)
  # POUR TOUS i TQ FRf(xx$y[i])<=FRf(0), THEN xx$y[i]=0
  # COMME C'EST DU QQ (i.e., FRf(xx$y[i])==FGf(DataGf2[i])), POUR TOUS i TQ FGf(DataGf2[i]))<=FRf(0), THEN xx$y[i]=0
  #                                                 <=>                              qntl[i]<=FRf(0), THEN xx$y[i]=0
  
  
  EmpGf2_th = (ecdf(DataGf2))(th_O)
  EmpGp2m1.Gf2_th = quantile(DataGp2, probs=EmpGf2_th, na.rm=TRUE)
  FRf_th = EmpFp.Gp2m1.Gf2_th = (ecdf(ObsRp_st))(EmpGp2m1.Gf2_th) ##### METTRE ObsRp_st pour avoir la proba qui evolue
  #cat("FRf_th =",FRf_th,"\n")
  xx$y[which(qntl<=(FRf_th))] = 0 ### A METTRE ABSOLUMENT
  xx$y[which(xx$y<=th_O)] = 0 ### NE CHANGE RIEN MAIS A METTRE PAR SECURITE
  
  #cat("length(which(xx$y<=th)) =",length(which(xx$y<=th)),"\n")
  #################
  #################
  
  
  
  #######################################################################################
  FGp=ecdf(DataGp)
  FGf=ecdf(DataGf)
  FGP=FGp(x)
  FGF=FGf(x)
  
  return(list(x=x,FRp=FRP,FGp=FGP,FGf=FGF,FRf=FRF,FRf_0=FRf_th,DS=xx$y))
  
}

CDFt_for_PR<-function(Rc,Mc,Mp,npas,th_O=NaN, th_M=NaN){
  res_PR=rep(NaN,length=length(Rc))
  test_bug=tryCatch({
    CDFt_SSR(Rc,Mc,Mp,npas=npas,th_O=th_O, th_M = th_M)$DS
  }, error=function(e){})
  if(is.null(test_bug)){
    res_PR<-CDFt_SSRbis(Rc,Mc,Mp,npas=npas,th_O=th_O, th_M = th_M)$DS
  }
  else{
    res_PR<-test_bug
  }
  return(list(DS=res_PR))
}
#### end CDFt_SSR as in Francois2020
