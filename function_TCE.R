library(fields)
library(igraph)
library(ismev)
library(energy)
library(MASS)
library(evd)
library(VineCopula)
library(goftest)
library(copula)
library(viridis)
library(RColorBrewer)
library(scales)

flatten_array <- function(data_array){ #data_array of type LON x LAT x TIME; output: TIME x VARPHY
  res = matrix(NaN, ncol= dim(data_array)[1] * dim(data_array)[2], nrow = dim(data_array)[3])
  k=0
  for(i in 1:dim(data_array)[1]){
    for(j in 1:dim(data_array)[2]){
      k=k+1
      res[,k]<-data_array[i,j,]
    }
  }
  return(res)
}

convert_matrix_to_array <- function(data_mat,nb_lon = 2, nb_lat = 2){ #data_mat of type TIME x VARPHY; output: LON x LAT x TIME
  res = array(NaN, dim = c(nb_lon, nb_lat, nrow(data_mat)))
  k=0
  for(i in 1:nb_lat){
    for(j in 1:nb_lon){
      k=k+1
      res[i,j,]<-data_mat[,k]
    }
  }
  return(res)
}

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


# CvM mail Mathieu 25.05.2021 
########## Cette fonction est appelee par la fonction "cvm" et ne doit pas etre appelee par l'utilisateur.
########## Elle est a sourcer
CramerVonMisesTwoSamples = function (S1, S2){
  xS1 = sort(S1)
  M = length(xS1)
  xS2 = sort(S2)
  N = length(xS2)
  a = data.frame(val = xS1, rang = seq(M), ens = rep(1, M))
  b = data.frame(val = xS2, rang = seq(N), ens = rep(2, N))
  c = rbind(a, b)
  c = c[order(c$val), ]
  c = data.frame(c, rangTot = seq(M + N))
  dtfM = c[which(c$ens == 1), ]
  dtfN = c[which(c$ens == 2), ]
  somN = sum((dtfN$rang - dtfN$rangTot)^2)
  somM = sum((dtfM$rang - dtfM$rangTot)^2)
  U = N * somN + M * somM
  CvM = ((U/(N * M))/(N + M)) - ((4 * M * N - 1)/(6 * (M + N)))
  return(CvM)
}

################
cvm = function(TS1, TS2){
  # compute 2 things:
  # - value of the CvM "distance" T
  # and
  # - p-value of the CvM test performed between time series T1 and T2
  # pval must be < alpha to reject "HO: equality of distribution of the 2 datasets" with a confidence (1-alpha)*100 % (if alpha= 0.05, confidence = 95%)
  # pval > alpha => cannot reject equality => distributions are the same
  # pval < alpha => reject equality        => ditributions are significantly different
  
  L = length(TS1) # je fais l'hypothese que col1 et col2 ont la meme longueur. Si ce n'est pas la meme longueur mais proche ca ne change rien.
  T = CramerVonMisesTwoSamples(TS1,TS2)  # from CDFt
  pval = pCvM(T, n = L, lower.tail = FALSE) # from goftest
  return(list(pval=pval,T=T))
}

## Ici je genere des donnees aleatoirement. En pratique, tu mettras tes donnees.
#Obs = rnorm(10000)
#Model = rnorm(10000)
#
## calcul de cvm
#CVM = cvm(Obs, Model)
#CVM$T # donne la valeur de la "distance" entre pdfs
#CVM$pval # p-value of the CvM test
#
#alpha = 0.05 # Fixe confidence at 95%
##### conclusion
#if(CVM$pval>=alpha){
#  cat("pval > alpha => cannot reject equality => distributions are the same")
#}
#if(CVM$pval<alpha){
#  cat("pval < alpha => reject equality        => ditributions are significantly different (at the chosen confidence level)")
#}
#
 
compute_daily_spat_sum<-function(daily_data,Ind_season){
  nb_lat=dim(daily_data)[1]
  nb_lon=dim(daily_data)[2]
  tmp_daily_data <- t(matrix(daily_data, nrow = nrow(daily_data)
                                  * ncol(daily_data)))
  spavg_daily_data <- apply(tmp_daily_data, 1, sum, na.rm=T) #Spatial average for the season
  
  return(spavg_daily_data)
}

compute_daily_spat_max<-function(daily_data,Ind_season){
  nb_lat=dim(daily_data)[1]
  nb_lon=dim(daily_data)[2]
  tmp_daily_data <- t(matrix(daily_data, nrow = nrow(daily_data)
                                  * ncol(daily_data)))
  spavg_daily_data <- apply(tmp_daily_data, 1, max, na.rm=T) #Spatial average for the season
  
  return(spavg_daily_data)
}

compute_daily_spat_mean<-function(daily_data,Ind_season){
  nb_lat=dim(daily_data)[1]
  nb_lon=dim(daily_data)[2]
  tmp_daily_data <- t(matrix(daily_data, nrow = nrow(daily_data)
                                  * ncol(daily_data)))
  spavg_daily_data <- apply(tmp_daily_data, 1, mean, na.rm=T) #Spatial average for the season
  
  return(spavg_daily_data)
}

## More or less like Vautard2021: for each year, minimum of tas at each grid cell + mean
compute_monthly_mean_of_mins<-function(data_month_array, nb_years){
  nb_day_by_month=dim(data_month_array)[3]/nb_years
  res=c()
  for(k in 1:nb_years){
    ### min in april for each year
    tmp_min_by_gridcells<-
      apply(data_month_array[,,(1+(k-1)*nb_day_by_month):(nb_day_by_month+(k-1)*nb_day_by_month)],c(1,2),min)
    ### spatial mean of min by year
    res=c(res, mean(tmp_min_by_gridcells, na.rm=TRUE))
  }
  return(res)
}

## Like Vautard: Solstice d'hiver for starting date
compute_spat_mean_of_GDD<-function(data_array, nb_years, month_studied_){
  nb_day_by_month=dim(data_array)[3]/nb_years
  res=c()
  for(k in 1:nb_years){
    #tmp_days=c(355:365,1:90)+365*(k-1)
    if(month_studied==3){
      print("GDD from jan. to end of feb.")
      tmp_days=c(1:59)+365*(k-1)}
    if(month_studied==4){
      print("GDD from jan. to end of march")
     tmp_days=c(1:90)+365*(k-1)}
    #tmp_days=c(1:90)+365*(k-1)
    vec_to_sum=data_array[,,tmp_days]-5
    vec_to_sum[vec_to_sum<0]<-0
    gdd=apply(vec_to_sum,c(1,2),cumsum)
    ## Get mean of GDD at the end of March
    #res_spatmean=mean(gdd[90,,],na.rm=T)
    res_spatmean=mean(gdd[length(tmp_days),,],na.rm=T)
    res=c(res, res_spatmean)
  }
  return(res)
}


#### Fit kernel and estimate its inverse
F_Kernel<-function(data_x, points_to_estimate){
  pdf <- density(data_x)
  # Interpolate the density
  f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
  # Get the cdf by numeric integration
  cdf <- function(x){
    integrate(f, -Inf, x)$value
  }
  p_of_x=c()
  for(i in points_to_estimate){
    p_of_x=c(p_of_x,cdf(i))
  }
  # # Use a root finding function to invert the cdf
  # invcdf <- function(p){
  #   uniroot(function(data_x){cdf(data_x) - p}, range(data_x))$root
  # }
  # q_of_x=c()
  # for(j in p_of_x){
  #   q_of_x=c(q_of_x, invcdf(j))
  # }
  # print(paste0("Max of diff.: ", max(data_x-q_of_x)))
  # plot(q_of_x,data_x)
  # abline(a=0,b=1)
  return(p_of_x)
}

Fm1_Kernel<-function(data_x, probs_to_estimate){
  pdf <- density(data_x)
  # Interpolate the density
  f <- approxfun(pdf$x, pdf$y, yleft=0, yright=0)
  # Get the cdf by numeric integration
  cdf <- function(x){
    integrate(f, -Inf, x)$value
  }
  p_of_x=c()
  for(i in data_x){
    p_of_x=c(p_of_x,cdf(i))
  }
  # Use a root finding function to invert the cdf
  invcdf <- function(p){
    uniroot(function(data_x){cdf(data_x) - p}, range(data_x))$root
  }
  q_of_x=c()
  for(j in probs_to_estimate){
    q_of_x=c(q_of_x, invcdf(j))
  }
  return(q_of_x)
}

#### Fit GEV and estimate its inverse
F_GEV<-function(data_x, points_to_estimate){
  ### Fit F_hat with data_x
  (fitted <- fgev(data_x))
  print(AIC(fitted))
  ### Estimate F_hat(points_to_estimate)
  p_of_x=evd::pgev(points_to_estimate,loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  return(p_of_x)
}

Fm1_GEV<-function(data_x, probs_to_estimate){ 
  ### Fit F_hat with data_x
  (fitted <- fgev(data_x))
  print(AIC(fitted))
  ### Estimate Fm1_hat
  q_of_x=evd::qgev(probs_to_estimate, loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  return(q_of_x)
}

#### FIT GPD over thresh and its inverse
F_GPD<-function(data_x, thresh_GPD, points_to_estimate){
  ### Fit F_hat with data_x
  # can have pb as: "observed information matrix is singular; use std.err = FALSE"
  test_fit=tryCatch({
    (fitted<- evd::fpot(data_x, thresh_GPD))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- evd::fpot(data_x, thresh_GPD, std.err=FALSE))
  }else{
    fitted=test_fit
  }
#  print(AIC(fitted))
  print(paste0("loc: ", thresh_GPD, ", scale: ", fitted$estimate[1], ", shape:", fitted$estimate[2]))
  ### Estimate F_hat(points_to_estimate)
  p_of_x=evd::pgpd(points_to_estimate,loc =thresh_GPD,
              scale=fitted$estimate[1], shape=fitted$estimate[2])
  return(p_of_x)
}

Fm1_GPD<-function(data_x, thresh_GPD, probs_to_estimate){
  ### Fit F_hat with data_x
  test_fit=tryCatch({
    (fitted<- evd::fpot(data_x, thresh_GPD))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- evd::fpot(data_x, thresh_GPD, std.err=FALSE))
  }else{
    fitted=test_fit
  }
#  print(AIC(fitted))
  ### Estimate Fm1_hat
  q_of_x=evd::qgpd(probs_to_estimate, loc =thresh_GPD,
                   scale=fitted$estimate[1], shape=fitted$estimate[2])
  return(q_of_x)
}

boot_param_GPD<-function(data_x, thresh_GPD, points_to_estimate, B, level_){
  ### Fit F_hat with data_x
  test_fit=tryCatch({
    (fitted<- evd::fpot(data_x, thresh_GPD))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- evd::fpot(data_x, thresh_GPD, std.err=FALSE))
  }else{
    fitted=test_fit
  }
  
  prof_68=confint(fitted, level=level_)
  min_scale=prof_68[1] 
  max_scale=prof_68[3] 
  min_shape=prof_68[2] 
  max_shape=prof_68[4] 
  res_boot_p_of_x=matrix(NaN, nrow=length(points_to_estimate), ncol=B) 
  ### Estimate F_hat(points_to_estimate)
  for(b in 1:B){
    b_scale=runif(1, min_scale, max_scale)
    b_shape=runif(1, min_shape, max_shape)
    res_boot_p_of_x[,b]=evd::pgpd(points_to_estimate,loc =thresh_GPD,
              scale=b_scale, shape=b_shape)
  }
  return(res_boot_p_of_x)
}

QQb_new_GPD<-function(Rc,Mc,Mp, thresh_GPD_Rc, thresh_GPD_Mc){ 
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
    fitted_Mc=evd::fpot(Mc[,k], thresh_GPD_Mc, std.err=FALSE)
    FMC=evd::pgpd(Mc[,k],loc =thresh_GPD_Mc,
                 scale=fitted_Mc$estimate[1],
                 shape=fitted_Mc$estimate[2])
    fitted_Rc=evd::fpot(Rc[,k], thresh_GPD_Rc, std.err=FALSE)
    Mch[,k]<-evd::qgpd(FMC,loc =thresh_GPD_Rc,
                       scale=fitted_Rc$estimate[1],
                       shape=fitted_Rc$estimate[2])
  }
  return(list(Mch=Mch))
}


####  Fit gaussian and estimate its inverse
F_Gaussian<-function(data_x, points_to_estimate){
  fitted=fitdistr(data_x, "normal")
  print(AIC(fitted))
  p_of_x=pnorm(points_to_estimate, mean = fitted$estimate[1], sd = fitted$estimate[2])
  return(p_of_x)
}

boot_param_Gaussian<-function(data_x, points_to_estimate, B, level_){
  fitted=fitdistr(data_x, "normal")
  prof_68=confint(fitted, level=level_)
  min_mean=prof_68[1]
  max_mean=prof_68[3]
  min_sd=prof_68[2] 
  max_sd=prof_68[4] 
  res_boot_p_of_x=matrix(NaN, nrow=length(points_to_estimate), ncol=B) 
  ### Estimate F_hat(points_to_estimate)
  for(b in 1:B){
    b_mean=runif(1, min_mean, max_mean)
    b_sd=runif(1, min_sd, max_sd)
    res_boot_p_of_x[,b]=pnorm(points_to_estimate, mean=b_mean,
              sd=b_sd)
  }
  return(res_boot_p_of_x)
}


Fm1_Gaussian<-function(data_x, probs_to_estimate){ 
  fitted=fitdistr(data_x, "normal")
  print(AIC(fitted))
  q_of_x=qnorm(probs_to_estimate, mean = fitted$estimate[1], sd = fitted$estimate[2])
  return(q_of_x)
}

#### Fit GEV and estimate its inverse
F_GEV_for_min<-function(data_x, points_to_estimate){
  #### Reverse data
  minus_data_x=(-1)*data_x
  minus_points_to_estimate=(-1)*points_to_estimate
  ### Fit F_hat with data_x
  # can have pb as: "observed information matrix is singular; use std.err = FALSE"
  test_fit=tryCatch({
    (fitted <- fgev(minus_data_x))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- fgev(minus_data_x, std.err=FALSE))
  }else{
    fitted=test_fit
  }
  ### Estimate 1-F_hat(points_to_estimate) (1- is here to reverse again)
  p_of_x= 1- evd::pgev(minus_points_to_estimate,loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  return(p_of_x)
}

boot_param_GEV_for_min<-function(data_x, points_to_estimate, B, level_){
  #### Reverse data
  minus_data_x=(-1)*data_x
  minus_points_to_estimate=(-1)*points_to_estimate
  ### Fit F_hat with data_x
  # can have pb as: "observed information matrix is singular; use std.err = FALSE"
  test_fit=tryCatch({
    (fitted <- fgev(minus_data_x))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- fgev(minus_data_x, std.err=FALSE))
  }else{
    fitted=test_fit
  }
  prof_68=confint(fitted, level=level_)
  min_loc=prof_68[1]
  max_loc=prof_68[4]
  min_scale=prof_68[2] 
  max_scale=prof_68[5] 
  min_shape=prof_68[3] 
  max_shape=prof_68[6] 
  res_boot_p_of_x=matrix(NaN, nrow=length(points_to_estimate), ncol=B) 
  ### Estimate F_hat(points_to_estimate)
  for(b in 1:B){
    b_loc=runif(1, min_loc, max_loc)
    b_scale=runif(1, min_scale, max_scale)
    b_shape=runif(1, min_shape, max_shape)
    res_boot_p_of_x[,b]=1-evd::pgev(minus_points_to_estimate,loc =b_loc,
              scale=b_scale, shape=b_shape)
  }
  return(res_boot_p_of_x)
}


Fm1_GEV_for_min<-function(data_x, probs_to_estimate){ 
  #### Reverse data
  minus_data_x=(-1)*data_x
  ### Fit F_hat with data_x
  # can have pb as: "observed information matrix is singular; use std.err = FALSE"
  test_fit=tryCatch({
    (fitted <- fgev(minus_data_x))
  }, error=function(e){})
  if(is.null(test_fit)){
      (fitted <- fgev(minus_data_x, std.err=FALSE))
  }else{
    fitted=test_fit
  }

 ### Estimate Fm1_hat
  q_of_x=evd::qgev(1-probs_to_estimate, loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  minus_q_of_x=(-1)*q_of_x
  return(minus_q_of_x)
}
#data_x=evd::rgev(50,loc =0, scale=1, shape=0)
#F_GEV(data_x, data_x)
#Fm1_GEV(data_x, c(0.01,0.99))
#
#
#data_x=rnorm(50,0,4)
#F_Gaussian(data_x, data_x)
#Fm1_Gaussian(data_x, c(0.05,0.5,0.97))
#
#
#F_Kernel(data_x, data_x)
#Fm1_Kernel(data_x, c(0.05,0.5, 0.97))


testIndepCopula<-function(data1, data2, pobs_to_calc=FALSE){
  if(pobs_to_calc==TRUE){
    pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
    pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
  }else{
    pobs_data1<-data1
    pobs_data2<-data2
  }
  a=BiCopIndTest(pobs_data1,pobs_data2)
  return(a)
}

reord_in_list_winter<-function(index, label_sw, nb_index_by_year,
                        length_sw){
  res_list=list()
  k=0
  for(l in 1:length(label_sw)){
    if(l==1){### 29 years of index and JanFeb for the first window
      start_=1
      end_=nb_index_by_year*(length_sw-1)+59
      res_list[[label_sw[l]]]=index[start_:end_]
    }
    if(l==length(label_sw)){ #29 years of index + Dec. 2100
      start_=length(index)-(nb_index_by_year*(length_sw-1))-31+1
      end_= length(index)
      res_list[[label_sw[l]]]=index[start_:end_]
    }
    if(l!=1 && l!= length(label_sw)){ ## 30 years of index
      start_=(k*nb_index_by_year+1-31)
      end_=start_+nb_index_by_year*length_sw-1
      res_list[[label_sw[l]]]=index[start_:end_]
      #res_list[[label_sw[l]]]=index[(k*nb_index_by_year+1-31):(k*nb_index_by_year+nb_index_by_year*length_sw)]
    }
    k=k+1
  } 
  return(res_list)
}



########## COPULA ##################



selectedCopula3<-function(data_U, pobs_to_calc=FALSE){
  AIC=list()
  par=list()
  logLik=list()
  family=c("Clayton", "Gumbel", "Frank", "Joe")
  for(f in family){
    tmp_res=determine_theta3(data_U, f, pobs_to_calc=FALSE)
    AIC[[f]]=2-2*tmp_res$logLik
    par[[f]]=tmp_res$par
    logLik[[f]]=tmp_res$logLik
  }
  res=list()
  res[["family"]]=names(which.min(AIC))
  res[["par"]]=as.numeric(par[[names(which.min(AIC))]])
  res[["logLik"]]=as.numeric(logLik[names(which.min(AIC))])
#  print("par")
#  print(unlist(par))
#  print("AIC")
#  print(unlist(AIC))
  return(res)
}

determine_theta3<-function(data_U, family, pobs_to_calc=FALSE){
  if(pobs_to_calc==TRUE){
    data_U<-pobs(as.matrix(cbind(data_U)))
  }
  res=list()
  res[["family"]]=family
  if(family=="Clayton"){
    cop_object3 <- claytonCopula(dim=2)
    number_f=3
    par_lim_inf=-0.9999
    }
  if(family=="Gumbel"){
    cop_object3 <- gumbelCopula(dim=2)
    number_f=4
    par_lim_inf=1.0001
  }
  if(family=="Frank"){
    cop_object3 <- frankCopula(dim=2)
    number_f=5
    par_lim_inf=-Inf
  }
  if(family=="Joe"){
    cop_object3 <- joeCopula(dim=2)
    number_f=6
    par_lim_inf=1.0001
  }
  
  test_fit_BFGS=tryCatch({
    fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS") #"Nelder-Mead"
  }, error=function(e){})
  
  if(is.null(test_fit_BFGS)){
     test_fit_NelderMead=tryCatch({
      fit=fitCopula(cop_object3, data_U, method="ml", optim.method="Nelder-Mead") #"Nelder-Mead"
    }, error=function(e){})
    if(is.null(test_fit_NelderMead)){
      print("determine_theta done with VineCopula")
      est_BiCop=BiCopEst(data_U[,1],data_U[,2],number_f,method = "mle")
      res[["logLik"]]= est_BiCop$logLik
      res[["par"]]= est_BiCop$par
    }else{
      print("determine_theta done with Nelder-Mead")
      fit=test_fit_NelderMead
      #fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
      res[["par"]]= fit@estimate
      res[["logLik"]]= fit@loglik
    }
  }else{
    print("determine_theta done with BFGS")
    fit=test_fit_BFGS
    res[["par"]]= fit@estimate
    res[["logLik"]]= fit@loglik
  }
  if(res[["par"]]<=par_lim_inf){res[["par"]]<-par_lim_inf}
  return(res)
}

testCopula3<-function(data_U, family,par_, pobs_to_calc=FALSE){
  name_family_lower=tolower(family)
  if(family=="Clayton"){
    cop_object3 <- claytonCopula(par_, dim=2)
    number_f=3
  }
  if(family=="Gumbel"){
    cop_object3 <- gumbelCopula(par_, dim=2)
    number_f=4
  }
  if(family=="Frank"){
    cop_object3 <- frankCopula(par_, dim=2)
    number_f=5
  }
  if(family=="Joe"){
    cop_object3<- joeCopula(par_, dim=2)
    number_f=6
  }
  if(pobs_to_calc==TRUE){
    data_U<-pobs(as.matrix(cbind(data_U)))
  }
  test_gofWhite=tryCatch({
    BiCopGofTest(data_U[,1],data_U[,2], number_f, par = par_, method = "white",max.df = 30, B = 100, obj = NULL)
    # gofWhite(name_family_lower, data_U, param = par_, param.est = FALSE, M = 100)
  },error=function(e){
  },finally={})
  if(is.null(test_gofWhite)){
    print("testCopula for VineCopula failed -> package copula 'BFGS'")
    #### sometimes, random bug due to bootstrap, iteration until it works 
    boolError<-TRUE
    k=0
    while(boolError==TRUE && k<=9)
    {
      k=k+1
      tryCatch({
        res <- gofCopula(cop_object3, data_U, N = 100,  estim.method="ml",verbose=FALSE,
                         optim.method="BFGS")#optimize())#"Nelder-Mead")#,optim.method="BFGS")
        boolError<-FALSE
      },error=function(e){
      },finally={})
    }
    ### If no gof can be computed after several tries,
    if(boolError==TRUE && k>=9){
      print("testCopula for copula 'BFGS' failed -> package copula 'Nelder-Mead'")
      
      boolError<-TRUE
      k=0
      while(boolError==TRUE && k<=9)
      {
        k=k+1
        tryCatch({
          res <- gofCopula(cop_object3, data_U, N = 100,  estim.method="ml",verbose=FALSE,
                           optim.method="Nelder-Mead")#optimize())#"Nelder-Mead")#,optim.method="BFGS")
          boolError<-FALSE
        },error=function(e){
        },finally={})
      }
      ###pval is assumed to be 0
      if(boolError==TRUE && k>=9){
        print("testCopula for all packages failed -> pval=0")
        res=list("p.value"=0)
        }
    }
    #print(res)
  }else{
    #get pval
    #print(test_gofWhite)
    res=list("p.value"=test_gofWhite$p.value)
    # res=list("p.value"=test_gofWhite[[name_family_lower]]$res.tests[1])
  }
  return(res)
}
#### si Profile Likelihood bug, on fait à la main
#Hofert 2011:  Likelihood inference for Archimedean copulas
LogL <- function(theta, acop, u, n.MC=0, ...) { # -(log-likelihood)
  sum(acop@dacopula(u, theta, n.MC=n.MC, log=TRUE, ...))
}

#### Three tries to compute C.I. for theta:
####1.  Package copula, 2. Package VineCopula 3. by hand
determine_ci_theta2<-function(cop_object, data_U, ci_level=0.95){
   print(paste0("CI ", ci_level, ":"))
  l_mle= cop_object[["logLik"]]
  theta_mle= cop_object[["par"]]
  
  d=2
  if(cop_object[["family"]]=="Clayton"){
    name_family="Clayton"
    par_lim_inf=-0.9999
    cop_object3 <- claytonCopula(theta_mle, dim=2)
  }
  if(cop_object[["family"]]=="Gumbel"){
    name_family="Gumbel"
    par_lim_inf=1.0001
    cop_object3 <- gumbelCopula(theta_mle, dim=2)
  }
  if(cop_object[["family"]]=="Frank"){
    name_family="Frank"
    par_lim_inf=theta_mle-3
    cop_object3 <- frankCopula(theta_mle, dim=2)
  }
  if(cop_object[["family"]]=="Joe"){
    name_family="Joe"
    par_lim_inf=1.0001
    cop_object3 <- joeCopula(theta_mle, dim=2)
  }
  cop_object2 <- onacopulaL(name_family, list(theta_mle,1:d))
  
  ### V1 
  test_efm=tryCatch({
    fit.tau <- fitCopula(cop_object3, data_U, method="ml",
                         optim.method="BFGS")
    ci=confint(fit.tau, level=ci_level)
  }, error=function(e){})
  if(is.null(test_efm)){
    ci=c( "MLE"=theta_mle, c(NaN,NaN))
  }else{
    fit.tau <- fitCopula(cop_object3, data_U, method="ml",
                         optim.method="BFGS")
    ci=confint(fit.tau, level=ci_level)
    ci<-c("MLE"=theta_mle, ci)
  }
  print(paste0("v1: ", round(ci[1],4), ", [", round(ci[2],4), ", ",
               round(ci[3],4), "]"))
 
  ### V2 Si ca bug...: 
  if(sum(is.na(ci))>0){
    test_efm=tryCatch({
      efm <- emle(data_U, cop_object2)
      pfm <- profile(efm)
    }, error=function(e){})
    if(is.null(test_efm)){
      ci=c( "MLE"=theta_mle, c(NaN,NaN))
    }else{
      efm <- emle(data_U, cop_object2)
      pfm <- profile(efm)
      ci  <- confint(pfm, level=ci_level)
      if(length(ci)==0){
        ci<-c("MLE"=theta_mle, NaN, NaN)
      }else{
        ci<-c("MLE"=theta_mle, ci)
      }
   }
  }
  print(paste0("v2: ", round(ci[1],4), ", [", round(ci[2],4), ", ",
               round(ci[3],4), "]"))
  
  ### V3: si ci comporte toujours Na...
  if(sum(is.na(ci))>0){
    if(is.na(ci[2])){
      print("ci_low a la main")
      th4=seq(theta_mle,par_lim_inf-0.0001,by=-0.0001) ### on va de theta mle a la borne inf
      Chi=qchisq(ci_level, df=1)
      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
      LR=2*(l_mle-Lt1)
      k=1
      while(LR<=Chi & th4[k]>=par_lim_inf){
        k=k+1
        # print(k)
        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
        LR=2*(l_mle-Lt1)
      }
      if(k==length(th4)){
        ci[2]=par_lim_inf
      }else{
        ci[2]=th4[k]
      }
      print(paste0("CI comp: ", ci[2]))
    }
    if(is.na(ci[3])){
      print("ci_high a la main")
      th4=seq(theta_mle,theta_mle+3,by=0.0001)
      Chi=qchisq(ci_level, df=1)
      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
      LR=2*(l_mle-Lt1)
      k=1
      while(LR<=Chi){
        k=k+1
        # print(k)
        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
        LR=2*(l_mle-Lt1)
      }
      ci[3]=th4[k]
      print(paste0("CI comp: ", ci[3]))
    }
  }
  print(paste0("v3: ", round(ci[1],4), ", [", round(ci[2],4), ", ",
               round(ci[3],4), "]"))
  ### to avoid cases where ci goes beyong lim inf
  if(ci[2]<par_lim_inf & name_family %in% c("Frank", "Joe", "Gumbel", "Clayton")){ci[2]<-par_lim_inf}
  print(paste0("vf: ", round(ci[1],4), ", [", round(ci[2],4), ", ",
               round(ci[3],4), "]"))
  return(ci)
}


### A generaliser et checker
determine_ci_prob_bivar_array2<-function(cop_object, ci_theta, prob_vector1=NaN,
                                        prob_vector2=NaN, zone_bivar=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
  ### zone_bivar: cadran bivarie among "topright" and "topleft"
  name_family=cop_object$family
  d=2
  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe

  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
  length_ci_theta=3
  ### loop over the different value of theta: MLE and ci_low and ci_high
  for(coord_theta in 1:length_ci_theta){
    test_fit_Cop_lowci=tryCatch({
      BiCop(cop_object$family, ci_theta[coord_theta])
    }, error=function(e){})
    if(is.null(test_fit_Cop_lowci)){
      fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[coord_theta],1:d))
      ## Evaluate this copula at prob_bivar
      if(zone_bivar=="topright"){
        tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
                                                       -pCopula(c(x[1], 1), fitted_cop_lowci)
                                                       +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
      }
      if(zone_bivar=="topleft"){
        tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci) 
                                                       -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
      }
        res_prob[,,coord_theta]<-tmp_mat 
    }else{
      fitted_cop_lowci<-test_fit_Cop_lowci
      ## Evaluate this copula at prob_bivar
      if(zone_bivar=="topright"){
        tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
                                                       -BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                       + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
      }
      if(zone_bivar=="topleft"){
        tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                       - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
      }
      res_prob[,,coord_theta]<-tmp_mat 
    }
  }
  #### On inverse low_ci and high_ci for zone_bivar="topleft"
  if(zone_bivar=="topleft"){
    tmp_save_low=res_prob[,, 2] 
    tmp_save_high=res_prob[,, 3] 
    res_prob[,,2]<-tmp_save_high 
    res_prob[,,3]<-tmp_save_low 
  }
  ### For very small values:
  epsilon=.Machine$double.eps
  res_prob[res_prob<=epsilon]<-0
  return(res_prob)
}



determine_bootstrap_prob_bivar<-function(cop_object, theta, prob_vector1=NaN,
                                        prob_vector2=NaN, zone_bivar=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
  ### zone_bivar: cadran bivarie among "topright" and "topleft"
  name_family=cop_object$family
  d=2
  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe

  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 
  res_prob=matrix(NaN, nrow=length(prob_vector1), ncol=length(prob_vector2))
  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
  test_fit_Cop_lowci=tryCatch({
    BiCop(cop_object$family, theta)
  }, error=function(e){})
  if(is.null(test_fit_Cop_lowci)){
    fitted_cop_lowci <- onacopulaL(name_family, list(theta,1:d))
    ## Evaluate this copula at prob_bivar
    if(zone_bivar=="topright"){
      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
                                                     -pCopula(c(x[1], 1), fitted_cop_lowci)
                                                     +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    }
    if(zone_bivar=="topleft"){
      tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci) 
                                                     -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    }
    res_prob<-tmp_mat 
  }else{
    fitted_cop_lowci<-test_fit_Cop_lowci
    ## Evaluate this copula at prob_bivar
    if(zone_bivar=="topright"){
      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
                                                     -BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                     + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    }
    if(zone_bivar=="topleft"){
      tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                     - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    }
    res_prob<-tmp_mat 
  }
#  #### On inverse low_ci and high_ci for zone_bivar="topleft"
#  if(zone_bivar=="topleft"){
#    tmp_save_low=res_prob[,, 2] 
#    tmp_save_high=res_prob[,, 3] 
#    res_prob[,,2]<-tmp_save_high 
#    res_prob[,,3]<-tmp_save_low 
#  }
  ### For very small values:
  epsilon=.Machine$double.eps
  res_prob[res_prob<=epsilon]<-0
  return(res_prob)
}



compute_ToE<-function(proba, low_bound, high_bound, length_SlidingWindow,
                      period_){
    substrRight <- function(x, n){
      substr(x, nchar(x)-n+1, nchar(x))
    }
    period_max=as.numeric(substrRight(period_, 4)) ### is equal to 2100
    k=length(proba) 
    ## particular case for very small values....
    if(low_bound>=high_bound && low_bound<=10*.Machine$double.eps){
      low_bound=high_bound=0
    }
    if(proba[k]>=low_bound && proba[k]<=high_bound){
      ### no ToE
      TOE=NaN
      TOE_text="NaN"
      length_last_seq=NaN
      coord_TOE=NaN
    }else{
      while(proba[k]<low_bound | proba[k]>high_bound){
        k=k-1
      }
      #print(k)
      length_last_seq=length(proba)-(k)#seq_out_of_bound$length[length(seq_out_of_bound$values)]
      coord_TOE=length(proba)-length_last_seq+1#-(length_SlidingWindow-1)
      TOE=period_max-(length_SlidingWindow/2)+1-length_last_seq+1
       TOE_text=paste0(period_max-length_last_seq+1-length_SlidingWindow+1,
                       "-", period_max-length_last_seq+1)
      }
    #print(TOE)
    return(list(TOE=TOE, TOE_text=TOE_text, length_last_seq=length_last_seq, coord_TOE=coord_TOE))
}

compute_ToE_matrix_from_array<-function(array_res_proba, length_SlidingWindow,
                                        period_, coord_baseline){
  nb_index1=dim(array_res_proba)[2]
  nb_index2=dim(array_res_proba)[1]
  mat_TOE=matrix(NaN, ncol=nb_index1, nrow=nb_index2)
  mat_coord_TOE=matrix(NaN, ncol=nb_index1, nrow=nb_index2)
  mat_TOE_text=matrix(NaN, ncol=nb_index1, nrow=nb_index2)
  for(p_of_index2 in 1:nb_index2){
    for(p_of_index1 in 1:nb_index1){
    tmp_proba=array_res_proba[p_of_index2,p_of_index1,,1]
    tmp_low_bound=min(array_res_proba[p_of_index2,p_of_index1,coord_baseline,2])
    tmp_high_bound=max(array_res_proba[p_of_index2,p_of_index1,coord_baseline,3])
    tmp_TOE=compute_ToE(tmp_proba, tmp_low_bound, tmp_high_bound, length_SlidingWindow,
                period_)
    TOE=tmp_TOE$TOE
    TOE_coord_TOE=tmp_TOE$coord_TOE
    TOE_text=tmp_TOE$TOE_text
    #### Fill the matrices
    mat_TOE[p_of_index2, p_of_index1]=TOE
    mat_coord_TOE[p_of_index2, p_of_index1]=TOE_coord_TOE
    mat_TOE_text[p_of_index2, p_of_index1]=TOE_text
    }
  }
  return(list(TOE=mat_TOE, coord_TOE=mat_coord_TOE, TOE_text=mat_TOE_text))
}


postproc_image.plot<-function(mat){ # row and col in mat correspond to what is
#plotted by image.plot
   res_mat<-mat[nrow(mat):1,]
   res_mat=t(res_mat)
   res_mat<-res_mat[,ncol(res_mat):1]
 return(res_mat)
} 


substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}
substrLeft <- function(x, n){
    substr(x, 1, n)
}

compute_PerkinsSS<-function(data_x, data_y,nbins=1000){
  dx=density(data_x)
  dy=density(data_y)
  range=c(min(data_x,data_y), max(data_x,data_y))
  tnew=seq(range[1],range[2],length.out=nbins)
      
  pdf_x=approx(dx$x,dx$y,xout=tnew)$y*diff(tnew)[1]
  pdf_y=approx(dy$x,dy$y,xout=tnew)$y*diff(tnew)[1]
        
  pdf_x[is.na(pdf_x)]<-0
  pdf_y[is.na(pdf_y)]<-0
  # print(sum(pdf_x))
  # print(sum(pdf_y))
  PDF_X_Y=matrix(NaN,ncol=2, nrow=length(tnew))
  PDF_X_Y[,1]=pdf_x
  PDF_X_Y[,2]=pdf_y
             
  Perkins_SS=sum(apply(PDF_X_Y,1,min,na.rm=TRUE))
  return(Perkins_SS)
}



fastplot_toe<-function(array_prob,p_i2, p_i1){
    #tmp_toe=compute_ToE_matrix_from_array(array_prob,222, "1850_2100", 1:22)
    #image.plot(1:19, 1:19,tmp_toe$TOE,zlim=c(1950,2050))
    # max_=0.5
    plot(array_prob[p_i2,p_i1,,1],type='l',ylim=c(0,0.15))
    lines(array_prob[p_i2,p_i1,,2], col="red", lty=2)
    lines(array_prob[p_i2,p_i1,,3], col="red", lty=2)
    abline(h=min(array_prob[p_i2,p_i1,1:22,2]))
    abline(h=max(array_prob[p_i2,p_i1,1:22,3]))
             # abline(v=tmp_toe$coord_TOE[p_i2,p_i1])
}


###################################################################################################################################################
# old but good
#old_equivalent_bootstrap_determine_ponctual_prob_bivar<-function(cop_object, theta, prob_vector1=NaN,
#                                        prob_vector2=NaN, zone_bivar=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
#  ### zone_bivar: cadran bivarie among "topright" and "topleft"
#  name_family=cop_object$family
#  d=2
#  res_prob=NaN
#
#  test_fit_Cop=tryCatch({
#      BiCop(cop_object$family, theta)
#    }, error=function(e){})
#    if(is.null(test_fit_Cop)){
#      fitted_cop <- onacopulaL(name_family, list(theta,1:d))
#      ## Evaluate this copula at prob_bivar
#      if(zone_bivar=="topright"){
#        tmp_mat=(1-pCopula(c(1, prob_vector2),
#                           fitted_cop)-pCopula(c(prob_vector1, 1),
#                           fitted_cop)+pCopula(c(prob_vector1, prob_vector2), fitted_cop))
#      }
#      if(zone_bivar=="topleft"){
#        tmp_mat=(pCopula(c(prob_vector1,1),fitted_cop)-pCopula(c(prob_vector1, prob_vector2), fitted_cop))
#      }
##
##      if(zone_bivar=="topleft"){
##        tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci) 
##                                                       -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
##      }
#       res_prob<-tmp_mat 
#    }else{
#      fitted_cop<-test_fit_Cop
#      ## Evaluate this copula at prob_bivar
#      if(zone_bivar=="topright"){
#        tmp_mat=(1-BiCopCDF(1, prob_vector2, fitted_cop)
#                                                       -BiCopCDF(prob_vector1, 1, fitted_cop)
#                                                       + BiCopCDF(prob_vector1,
#                                                                  prob_vector2, fitted_cop))
#      }
#      if(zone_bivar=="topleft"){
#        tmp_mat=(BiCopCDF(prob_vector1, 1, fitted_cop) - BiCopCDF(prob_vector1, prob_vector2, fitted_cop))
#      }
##
##      
##      if(zone_bivar=="topleft"){
##        tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
##                                                       - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
##      }
#      res_prob<-tmp_mat 
#    }
##  #### On inverse low_ci and high_ci for zone_bivar="topleft"
##  if(zone_bivar=="topleft"){
##    tmp_save_low=res_prob[,, 2] 
##    tmp_save_high=res_prob[,, 3] 
##    res_prob[,,2]<-tmp_save_high 
##    res_prob[,,3]<-tmp_save_low 
##  }
#  ### For very small values:
#  epsilon=.Machine$double.eps
#  res_prob[res_prob<=epsilon]<-0
#  return(res_prob)
#}



#old_before_generalisation_determine_ci_prob_bivar_array2<-function(cop_object, ci_theta, prob_vector1=NaN,
#                                        prob_vector2=NaN, zone_bivar=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
#  ### zone_bivar: cadran bivarie among "topright" and "topleft"
#  name_family=cop_object$family
#  d=2
#  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe
#
#  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
#  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  ####Proba for low_ci
#  test_fit_Cop_lowci=tryCatch({
#    BiCop(cop_object$family, ci_theta[2])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_lowci)){
#    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
#    ## Evaluate this copula at prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_lowci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat 
#  }else{
#    fitted_cop_lowci<-test_fit_Cop_lowci
#    ## Evaluate this copula at prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_lowci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  
#  ####Proba for high_ci
#  test_fit_Cop_highci=tryCatch({
#    BiCop(cop_object$family, ci_theta[3])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_highci)){
#    #print("Use of pCopula")
#    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_highci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_highci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_highci<-test_fit_Cop_highci
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_highci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_highci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
# 
#  test_fit_Cop_mle=tryCatch({
#    BiCop(cop_object$family, ci_theta["MLE"])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_mle)){
#    #print("Use of pCopula")
#    fitted_cop_mle <- onacopulaL(name_family, list(ci_theta["MLE"],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]),
#                                                               fitted_cop_mle)
#                                                     -pCopula(c(x[1], 1),
#                                                              fitted_cop_mle)
#                                                     +pCopula(c(x[1], x[2]),
#                                                              fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_mle<-test_fit_Cop_mle
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2],
#                                                                fitted_cop_mle)
#                                                     -BiCopCDF(x[1], 1,
#                                                               fitted_cop_mle)
#                                                     + BiCopCDF(x[1], x[2],
#                                                                fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  return(res_prob)
#}
#
#
#determine_ci_prob_bivar_array<-function(cop_object, ci_theta, prob_vector1=NaN,
##                                        prob_vector2=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
##  if(cop_object$family==3){name_family="Clayton"}
#  if(cop_object$family==4){name_family="Gumbel"}
#  if(cop_object$family==5){name_family="Frank"}
#  if(cop_object$family==6){name_family="Joe"}
#  d=2
#  
#  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe
#
#  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
#  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  
#  ####Proba for low_ci
#  test_fit_Cop_lowci=tryCatch({
#    BiCop(cop_object$family, ci_theta[2])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_lowci)){
#    #print("Use of pCopula")
#    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_lowci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#    
#  }else{
#    fitted_cop_lowci<-test_fit_Cop_lowci
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_lowci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  
#  ####Proba for high_ci
#  test_fit_Cop_highci=tryCatch({
#    BiCop(cop_object$family, ci_theta[3])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_highci)){
#    #print("Use of pCopula")
#    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_highci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_highci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_highci<-test_fit_Cop_highci
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_highci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_highci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
# 
#  test_fit_Cop_mle=tryCatch({
#    BiCop(cop_object$family, ci_theta["MLE"])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_mle)){
#    #print("Use of pCopula")
#    fitted_cop_mle <- onacopulaL(name_family, list(ci_theta["MLE"],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]),
#                                                               fitted_cop_mle)
#                                                     -pCopula(c(x[1], 1),
#                                                              fitted_cop_mle)
#                                                     +pCopula(c(x[1], x[2]),
#                                                              fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_mle<-test_fit_Cop_mle
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2],
#                                                                fitted_cop_mle)
#                                                     -BiCopCDF(x[1], 1,
#                                                               fitted_cop_mle)
#                                                     + BiCopCDF(x[1], x[2],
#                                                                fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  return(res_prob)
#}
#determine_theta<-function(cop_object, data_U){
#  if(cop_object$family==3){
#    name_family="Clayton"
#    cop_object3 <- claytonCopula(dim=2)
#  }
#  if(cop_object$family==4){
#    name_family="Gumbel"
#    cop_object3 <- gumbelCopula(dim=2)
#  }
#  if(cop_object$family==5){
#    name_family="Frank"
#    cop_object3 <- frankCopula(dim=2)
#  }
#  if(cop_object$family==6){
#    name_family="Joe"
#    cop_object3 <- joeCopula(dim=2)
#  }
#  res_theta=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#  return(res_theta)
#}
#
#old_TardiveFrost_determine_ci_prob_bivar_array<-function(cop_object, ci_theta, prob_vector1=NaN,
#                                        prob_vector2=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
#  if(cop_object$family==3){name_family="Clayton"}
#  if(cop_object$family==4){name_family="Gumbel"}
#  if(cop_object$family==5){name_family="Frank"}
#  if(cop_object$family==6){name_family="Joe"}
#  d=2
#  
#  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe
#
#  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
#  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
#  
#  ####Proba for low_ci
#  test_fit_Cop_lowci=tryCatch({
#    BiCop(cop_object$family, ci_theta[2])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_lowci)){
#    #print("Use of pCopula")
#    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci) 
#                                                     -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#    
#  }else{
#    fitted_cop_lowci<-test_fit_Cop_lowci
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
#                                                     - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  
#  ####Proba for high_ci
#  test_fit_Cop_highci=tryCatch({
#    BiCop(cop_object$family, ci_theta[3])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_highci)){
#    #print("Use of pCopula")
#    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(pCopula(c(x[1], 1), fitted_cop_highci)
#                                                     -pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_highci<-test_fit_Cop_highci
#    ## Evaluate this copula at the vector prob_bivar
#    tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_highci)
#                                                     - BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  
#  fitted_cop_mle<-BiCop(cop_object$family, ci_theta["MLE"])
#  ## Evaluate this copula at the vector prob_bivar
#  tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_mle)
#                                                   - BiCopCDF(x[1], x[2], fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#  res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  return(res_prob)
#}
#
#determine_ci_theta<-function(cop_object, data_U, ci_level=0.95){
#  test_theta=tryCatch({
#    est_BiCop=determine_theta(cop_object, data_U) 
#  }, error=function(e){})
#  if(is.null(test_theta)){
#    est_BiCop=BiCopEst(data_U[,1],data_U[,2],cop_object$family,method = "mle")
#    l_mle= est_BiCop$logLik
#    theta_mle= est_BiCop$par
#  }else{
#    est_BiCop=determine_theta(cop_object, data_U) 
#    l_mle= est_BiCop@loglik
#    theta_mle= est_BiCop@estimate
#  }
#  
#  #est_BiCop=determine_theta(cop_object, data_U) 
#  #l_mle= est_BiCop@loglik
#  #theta_mle= est_BiCop@estimate
# 
#  # est_BiCop=BiCopEst(data_U[,1],data_U[,2],cop_object$family,method = "mle")
#  #l_mle= est_BiCop$logLik
#  #theta_mle= est_BiCop$par
#  d=2
#  if(cop_object$family==3){
#    name_family="Clayton"
#    par_lim_inf=-0.9999
#    cop_object3 <- claytonCopula(theta_mle, dim=2)
#  }
#  if(cop_object$family==4){
#    name_family="Gumbel"
#    par_lim_inf=1.0001
#    cop_object3 <- gumbelCopula(theta_mle, dim=2)
#  }
#  if(cop_object$family==5){
#    name_family="Frank"
#    par_lim_inf=theta_mle-3
#    cop_object3 <- frankCopula(theta_mle, dim=2)
#  }
#  if(cop_object$family==6){
#    name_family="Joe"
#    par_lim_inf=1.0001
#    cop_object3 <- joeCopula(theta_mle, dim=2)
#  }
#  cop_object2 <- onacopulaL(name_family, list(theta_mle,1:d))
#  
#  #### le choix de theta_mle est pas important dans onacopulaL ici (car on le
#  ### recalcule ensuite)
#  ## Estimation
#  #### efm peut soit crasher complet, soit a moitie..
#  test_efm=tryCatch({
#    fit.tau <- fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#    ci=confint(fit.tau, level=ci_level)
#  }, error=function(e){})
#  if(is.null(test_efm)){
#    ci=c( "MLE"=theta_mle, c(NaN,NaN))
#  }else{
#    fit.tau <- fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#    ci=confint(fit.tau, level=ci_level)
#    ci<-c("MLE"=theta_mle, ci)
#  }
#    print(paste0("version1: ", ci[1], ", [", ci[2], ", ", ci[3], "]"))
# 
#  ### Si ca bug...: 
#  if(sum(is.na(ci))>0){
#    test_efm=tryCatch({
#      efm <- emle(data_U, cop_object2)
#      pfm <- profile(efm)
#    }, error=function(e){})
#    if(is.null(test_efm)){
#      ci=c( "MLE"=theta_mle, c(NaN,NaN))
#    }else{
#      efm <- emle(data_U, cop_object2)
#      #summary(efm) # using bblme's 'mle2' method
#      ## Profile likelihood plot [using S4 methods from bbmle/stats4] :
#      pfm <- profile(efm)
#      ci  <- confint(pfm, level=ci_level)
#      if(length(ci)==0){
#        ci<-c("MLE"=theta_mle, NaN, NaN)
#      }else{
#        ci<-c("MLE"=theta_mle, ci)
#      }
#   }
#   }
#    print(paste0("version2: ", ci[1], ", [", ci[2], ", ", ci[3], "]"))
# 
#  ### Parfois, ci comporte toujours Na... on fait a la main:
#  if(sum(is.na(ci))>0){
#    print("ATTENTION, ci avec Na")
#    if(is.na(ci[2])){
#      print("ci_low a la main")
#      th4=seq(theta_mle,par_lim_inf-0.0001,by=-0.0001) ### on va de theta mle a la borne inf
#      Chi=qchisq(ci_level, df=1)
#      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
#      LR=2*(l_mle-Lt1)
#      k=1
#      while(LR<=Chi & th4[k]>=par_lim_inf){
#        k=k+1
#        # print(k)
#        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
#        LR=2*(l_mle-Lt1)
#      }
#      if(k==length(th4)){
#        ci[2]=par_lim_inf
#      }else{
#        ci[2]=th4[k]
#      }
#      print(paste0("CI comp: ", ci[2]))
#    }
#    if(is.na(ci[3])){
#      print("ci_high a la main")
#      th4=seq(theta_mle,theta_mle+3,by=0.0001)
#      Chi=qchisq(ci_level, df=1)
#      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
#      LR=2*(l_mle-Lt1)
#      k=1
#      while(LR<=Chi){
#        k=k+1
#        # print(k)
#        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
#        LR=2*(l_mle-Lt1)
#      }
#      ci[3]=th4[k]
#      print(paste0("CI comp: ", ci[3]))
#    }
#  }
#  ### to avoid cases where ci goes beyong lim inf
#  if(ci[2]<par_lim_inf & name_family %in% c("Frank", "Joe", "Gumbel", "Clayton")){ci[2]<-par_lim_inf}
#  print(paste0("version3: ", ci[1], ", [", ci[2], ", ", ci[3], "]"))
# 
#  return(ci)
#}
#
#selectedCopula<-function(data1,data2, pobs_to_calc=FALSE){ #from Clayton, Gumbel, Frank, Joe, ### not BB1, BB6, BB7, and BB8. 
#  if(pobs_to_calc==TRUE){
#    pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#    pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  }else{
#    pobs_data1<-data1
#    pobs_data2<-data2
#  }
#
#  a=BiCopSelect(pobs_data1,pobs_data2,familyset=c(3:6),indeptest
#                = FALSE,selectioncrit = "AIC", rotations=FALSE)
#  return(a)
#}
#
#testCopula<-function(data1, data2, family, par=0, par2=0, pobs_to_calc=FALSE){
#  if(pobs_to_calc==TRUE){
#    pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#    pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  }else{
#    pobs_data1<-data1
#    pobs_data2<-data2
#  }
#  a=BiCopGofTest(pobs_data1,pobs_data2, family, method = "white",max.df = 30, B = 100, obj = NULL)
#  return(a)
#}
#selectedCopula2<-function(data_U, pobs_to_calc=FALSE){ 
#  if(pobs_to_calc==TRUE){
#    data_U<-pobs(as.matrix(cbind(data_U)))
#  }
#  AIC=list()
#  par=list()
#  logLik=list()
#  family=c("Clayton", "Gumbel", "Frank", "Joe")
#  for(f in family){
#    if(f=="Clayton"){cop_object3 <- claytonCopula(dim=2)}
#    if(f=="Gumbel"){cop_object3 <- gumbelCopula(dim=2)}
#    if(f=="Frank"){cop_object3 <- frankCopula(dim=2)}
#    if(f=="Joe"){cop_object3 <- joeCopula(dim=2)}
#    
#    test_fit=tryCatch({
#      fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#    }, error=function(e){})
#    if(is.null(test_fit)){
#      print("selectedCopula with package copula failed -> VineCopula")
#      if(f=="Clayton"){number_f=3}
#      if(f=="Gumbel"){number_f=4}
#      if(f=="Frank"){number_f=5}
#      if(f=="Joe"){number_f=6}
#      
#      est_BiCop=BiCopEst(data_U[,1],data_U[,2],number_f,method = "mle")
#      l_mle= est_BiCop$logLik
#      theta_mle= est_BiCop$par
#      AIC[[f]]=2-2*l_mle
#      par[[f]]=theta_mle
#      logLik[[f]]=l_mle
#    }else{
#      fit=test_fit
#      #fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#      AIC[[f]]=2-2*fit@loglik
#      par[[f]]=fit@estimate
#      logLik[[f]]=fit@loglik
#    }
#  }
#  res=list()
#  res[["family"]]=names(which.min(AIC))
#  res[["par"]]=as.numeric(par[[names(which.min(AIC))]])
#  res[["logLik"]]=as.numeric(logLik[names(which.min(AIC))])
#  print("par")
#  print(unlist(par))
#  print("AIC")
#  print(unlist(AIC))
# return(res)
#}
#
#determine_theta2<-function(data_U, family, pobs_to_calc=FALSE){  
#  if(pobs_to_calc==TRUE){
#    data_U<-pobs(as.matrix(cbind(data_U)))
#  }
#  res=list()
#  res[["family"]]=family
#  if(family=="Clayton"){cop_object3 <- claytonCopula(dim=2)}
#  if(family=="Gumbel"){cop_object3 <- gumbelCopula(dim=2)}
#  if(family=="Frank"){cop_object3 <- frankCopula(dim=2)}
#  if(family=="Joe"){cop_object3 <- joeCopula(dim=2)}
#  
#  test_fit=tryCatch({
#    fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#  }, error=function(e){})
#  if(is.null(test_fit)){
#    print("determine_theta with package copula failed -> VineCopula")
#    if(family=="Clayton"){number_f=3}
#    if(family=="Gumbel"){number_f=4}
#    if(family=="Frank"){number_f=5}
#    if(family=="Joe"){number_f=6}
#    
#    est_BiCop=BiCopEst(data_U[,1],data_U[,2],number_f,method = "mle")
#    res[["logLik"]]= est_BiCop$logLik
#    res[["par"]]= est_BiCop$par
#   }else{
#      fit=test_fit
#     #fit=fitCopula(cop_object3, data_U, method="ml", optim.method="BFGS")
#      res[["logLik"]]= fit@loglik
#      res[["par"]]= fit@estimate
#  }
#  
#  return(res)
#}
# 
#
#testCopula2<-function(data_U, family, pobs_to_calc=FALSE){  
#  if(pobs_to_calc==TRUE){
#    data_U<-pobs(as.matrix(cbind(data_U)))
#  }
#  if(family=="Clayton"){cop_object3 <- claytonCopula(dim=2)}
#  if(family=="Gumbel"){ cop_object3 <- gumbelCopula(dim=2) }
#  if(family=="Frank"){ cop_object3 <- frankCopula(dim=2)}
#  if(family=="Joe"){cop_object3<- joeCopula(dim=2)}
#
#  #### sometimes, random bug due to bootstrap, iteration until it works 
#  boolError<-TRUE
#  k=0
#  while(boolError==TRUE && k<=20)
#  {
#    k=k+1
#    print(paste0("Essai test", k))
#    tryCatch({
#      res <- gofCopula(cop_object3, data_U, N = 100,  estim.method="ml",
#                optim.method="Nelder-Mead")#,optim.method="BFGS")
#      boolError<-FALSE      
#    },error=function(e){
#                     },finally={})
#  }
#  ### If no gof can be computed after several tries, 
#  ###pval is assumed to be 0
#  if(boolError==TRUE && k>=20){res=list("p.value"=0)}
#  return(res)
#}#################################################################################################################################################
#old_compute_ToE<-function(proba, low_bound, high_bound, length_SlidingWindow=30){
#  seq_out_of_bound=rle(diff(which(proba<=low_bound | proba>=high_bound)))
#  ### Si il y a 1 fois la valeur 1 pour le rle a la fin, c'est qu'il y a 2 valeurs qui excede le thresh d'affilee -> 2050-2099 et 2051-2100
#  bool= seq_out_of_bound$values[length(seq_out_of_bound$values)]==1
#  bool_outside_bound= (proba[length(proba)]>= high_bound
#                       | proba[length(proba)]<= low_bound)
#  if((length(bool)>0) && bool ==TRUE && bool_outside_bound==TRUE){
#    length_last_seq=seq_out_of_bound$length[length(seq_out_of_bound$values)]
#    coord_TOE=length(proba)-length_last_seq#-(length_SlidingWindow-1)
#    TOE=2100-length_last_seq-(length_SlidingWindow-1)
#    TOE_text=paste0(2100-length_last_seq-(length_SlidingWindow-1),
#                            "-",
#                            2100-length_last_seq)
#  }else{
#    TOE=NaN
#    TOE_text="NaN"
#    length_last_seq=NaN
#    coord_TOE=NaN
#  }
#  print(TOE)
# return(list(TOE=TOE, TOE_text=TOE_text, length_last_seq=length_last_seq,
#             coord_TOE=coord_TOE)) 
#}
#
##  test_fit_Cop_lowci=tryCatch({
#    BiCop(cop_object$family, ci_theta[2])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_lowci)){
#    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
#    ## Evaluate this copula at prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_lowci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci) 
#                                                     -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#      res_prob[,,2]<-tmp_mat 
#  }else{
#    fitted_cop_lowci<-test_fit_Cop_lowci
#    ## Evaluate this copula at prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_lowci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
#                                                     - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  
#  ####Proba for high_ci
#  test_fit_Cop_highci=tryCatch({
#    BiCop(cop_object$family, ci_theta[3])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_highci)){
#    #print("Use of pCopula")
#    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_highci)
#                                                     -pCopula(c(x[1], 1), fitted_cop_highci)
#                                                     +pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(pCopula(c(x[1], 1), fitted_cop_highci)
#                                                     -pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_highci<-test_fit_Cop_highci
#    ## Evaluate this copula at the vector prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_highci)
#                                                     -BiCopCDF(x[1], 1, fitted_cop_highci)
#                                                     + BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_highci)
#                                                     - BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
# 
#  test_fit_Cop_mle=tryCatch({
#    BiCop(cop_object$family, ci_theta["MLE"])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_mle)){
#    #print("Use of pCopula")
#    fitted_cop_mle <- onacopulaL(name_family, list(ci_theta["MLE"],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]),
#                                                               fitted_cop_mle)
#                                                     -pCopula(c(x[1], 1),
#                                                              fitted_cop_mle)
#                                                     +pCopula(c(x[1], x[2]),
#                                                              fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(pCopula(c(x[1], 1), fitted_cop_mle)
#                                                     -pCopula(c(x[1], x[2]),
#                                                              fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#      
#    }
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }else{
#    #print("Use of BiCop")
#    fitted_cop_mle<-test_fit_Cop_mle
#    ## Evaluate this copula at the vector prob_bivar
#    if(zone_bivar=="topright"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2],
#                                                                fitted_cop_mle)
#                                                     -BiCopCDF(x[1], 1,
#                                                               fitted_cop_mle)
#                                                     + BiCopCDF(x[1], x[2],
#                                                                fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    if(zone_bivar=="topleft"){
#      tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_mle)
#                                                   - BiCopCDF(x[1], x[2], fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
#    }
#    res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
#  }
#  return(res_prob)
#}
#


#old_determine_ci_theta<-function(cop_object, data_U, ci_level=0.95){
#  est_BiCop=BiCopEst(data_U[,1],data_U[,2],cop_object$family,method = "mle")
#  l_mle= est_BiCop$logLik
#  theta_mle= est_BiCop$par
#  d=2
#  if(cop_object$family==3){
#    name_family="Clayton"
#    par_lim_inf=-0.9999
#  }
#  if(cop_object$family==4){
#    name_family="Gumbel"
#    par_lim_inf=1.0001
#  }
#  if(cop_object$family==5){
#    name_family="Frank"
#    par_lim_inf=theta_mle-3
#  }
#  if(cop_object$family==6){
#    name_family="Joe"
#    par_lim_inf=1.0001
#  }
#  cop_object2 <- onacopulaL(name_family, list(theta_mle,1:d))
#  
#  #### le choix de theta_mle est pas important dans onacopulaL ici (car on le
#  ### recalcule ensuite)
#  ## Estimation
#  #### efm peut soit crasher complet, soit a moitie..
#  test_efm=tryCatch({
#    efm <- emle(data_U, cop_object2)
#    pfm <- profile(efm)
#  }, error=function(e){})
#  if(is.null(test_efm)){
#    ci=c( "MLE"=theta_mle, c(NaN,NaN))
#  }else{
#    efm <- emle(data_U, cop_object2)
#    summary(efm) # using bblme's 'mle2' method
#    ## Profile likelihood plot [using S4 methods from bbmle/stats4] :
#    pfm <- profile(efm)
#    ci  <- confint(pfm, level=ci_level)
#    ci<-c("MLE"=theta_mle, ci)
#  }
#  
#  ### Parfois, ci comporte Na... on fait a la main:
#  if(sum(is.na(ci))>0){
#    print("ATTENTION, ci avec Na")
#    if(is.na(ci[2])){
#      print("ci_low a la main")
#      th4=seq(theta_mle,par_lim_inf-0.0001,by=-0.0001) ### on va de theta mle a la borne inf
#      Chi=qchisq(ci_level, df=1)
#      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
#      LR=2*(l_mle-Lt1)
#      k=1
#      while(LR<=Chi & th4[k]>=par_lim_inf){
#        k=k+1
#        # print(k)
#        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
#        LR=2*(l_mle-Lt1)
#      }
#      if(k==length(th4)){
#        ci[2]=par_lim_inf
#      }else{
#        ci[2]=th4[k]
#      }
#      print(paste0("CI comp: ", ci[2]))
#    }
#    if(is.na(ci[3])){
#      print("ci_high a la main")
#      th4=seq(theta_mle,theta_mle+3,by=0.0001)
#      Chi=qchisq(ci_level, df=1)
#      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
#      LR=2*(l_mle-Lt1)
#      k=1
#      while(LR<=Chi){
#        k=k+1
#        # print(k)
#        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
#        LR=2*(l_mle-Lt1)
#      }
#      ci[3]=th4[k]
#      print(paste0("CI comp: ", ci[3]))
#    }
#  }
#  return(ci)
#}
#

### Brouillon
# 
# 
#
#sampling.median_CI <- function(n_length, n_drop, B = 2000) {
#  median_boot <- rep(NaN, B)
##  ### v1: median theorique
##  for (i in 1:B) {
##    median_boot[i] <- median(runif(n,0,1))
##  }
##  ci_high=quantile(median_boot,0.975)
##  ci_low=quantile(median_boot,0.025)
##  res=c(ci_low,ci_high)
#  ### v2
#  u_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  for (i in 1:B) {
#    median_boot[i] <- median(sample(u_thresh,n_drop))
#  }
#  ci95_high=quantile(median_boot,0.975)
#  ci95_low=quantile(median_boot,0.025)
#  ci90_high=quantile(median_boot,0.95)
#  ci90_low=quantile(median_boot,0.05)
#  ci99_high=quantile(median_boot,0.005)
#  ci99_low=quantile(median_boot,0.995)
#  res_ci95=c(ci95_low,ci95_high)
#  res_ci90=c(ci90_low,ci90_high)
#  res_ci99=c(ci99_low,ci99_high)
#  return(list(res_ci95=res_ci95, res_ci90=res_ci90, res_ci99=res_ci99))
#}
#
#
#
#
##### Evolution
#compute_evolution_Condcorr<-function(data1, data2, col_="black",
#                                     col_pval_="black", name_plot="to def.", label_HW="HW", 
#                                sens=1,recomp_rank=FALSE){
#
#  # pval Fisher 
#  p.val_fisher = function(r1,r2,n1,n2){
#    fisher= ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
#    p.value = (2*(1-pnorm(abs(fisher))))
#    return(p.value)
#  }
#
#  ### compute Condcorr for calib and proj the half of the whole ts
#  n_length=length(data1)/2
# 
#  u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  res_mat_cor_Calib=res_mat_cor_Proj=res_mat_pval_Calib = res_mat_pval_Proj= res_mat_pvalFisher=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#  res_mat_length_id_Calib=res_mat_length_id_Proj=matrix(rep(NaN,n_length*n_length),ncol=n_length)
# 
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  if(recomp_rank==FALSE){
#    pobs_data1_Calib=pobs_data1[1:n_length]
#    pobs_data2_Calib=pobs_data2[1:n_length]
#    pobs_data1_Proj=pobs_data1[(n_length+1):(2*n_length)]
#    pobs_data2_Proj=pobs_data2[(n_length+1):(2*n_length)]
#  }else{
#    pobs_data1_Calib=pobs(as.matrix(cbind(data1[1:n_length])))[,1]
#    pobs_data2_Calib=pobs(as.matrix(cbind(data2[1:n_length])))[,1]
#    pobs_data1_Proj=pobs(as.matrix(cbind(data1[(n_length+1):(2*n_length)])))[,1]
#    pobs_data2_Proj=pobs(as.matrix(cbind(data2[(n_length+1):(2*n_length)])))[,1]
#  }
#  
#  ### Calib first
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      if(sens==1){
#        id_Calib=which(pobs_data1_Calib>=u & pobs_data2_Calib>=v)
#      }
#     if(length(id_Calib)<=5){
#        next
#      }else{
#     ### subselection
#        pobs_data1_Calib_sub<-pobs_data1_Calib[id_Calib]#(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_Calib_sub<-pobs_data2_Calib[id_Calib]#(as.matrix(cbind(data2[id])))[,1]
#        res_mat_cor_Calib[count_u,count_v]=cor(pobs_data1_Calib_sub,pobs_data2_Calib_sub,method="spearman")
#        res_mat_pval_Calib[count_u, count_v]=cor.test(pobs_data1_Calib_sub,pobs_data2_Calib_sub,method="spearman")$p.value
#        res_mat_length_id_Calib[count_u, count_v]=length(id_Calib)
#      }
#    }
#  }
#  ### Proj then
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      if(sens==1){
#        id_Proj=which(pobs_data1_Proj>=u & pobs_data2_Proj>=v)
#      }
#     if(length(id_Proj)<=5){
#        next
#      }else{
#     ### subselection
#        pobs_data1_Proj_sub<-pobs_data1_Proj[id_Proj]#(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_Proj_sub<-pobs_data2_Proj[id_Proj]#(as.matrix(cbind(data2[id])))[,1]
#        res_mat_cor_Proj[count_u,count_v]=cor(pobs_data1_Proj_sub,pobs_data2_Proj_sub,method="spearman")
#        res_mat_pval_Proj[count_u, count_v]=cor.test(pobs_data1_Proj_sub,pobs_data2_Proj_sub,method="spearman")$p.value
#        res_mat_length_id_Proj[count_u, count_v]=length(id_Proj)
#      }
#    }
#  }
#  
#  ### Fisher test between calib and proj
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      res_mat_pvalFisher[count_u, count_v]=p.val_fisher(res_mat_cor_Calib[count_u, count_v],
#                                               res_mat_cor_Proj[count_u, count_v],
#                                               res_mat_length_id_Calib[count_u, count_v],
#                                               res_mat_length_id_Proj[count_u, count_v])
#    }
#  }
#  
#  plot(pobs_data1_Calib, pobs_data2_Calib, col= "blue", ylim=c(0,1),
#       xlim=c(0,1), ylab="", xlab="", main=name_plot)
#  points(pobs_data1_Proj, pobs_data2_Proj, col= "darkorange")
#  title(ylab="F(Drought)",xlab=paste0("F(",label_HW,")"),
#     cex.lab=1.2, line=2)
#  legend("topleft",pch=c(1,1),c("Period 1", "Period 2"),
#       col=c("blue","darkorange"))
#
#  #### Condcorr
#  #### before plotting, which value in Calib is not computed but computed in
#  #Proj period?
#  points_missingCalib=which((is.na(res_mat_cor_Calib)
#                            & !is.na(res_mat_cor_Proj)),arr.ind=TRUE)
#  points_missingCalibx=points_missingCalib[,1]
#  points_missingCaliby=points_missingCalib[,2]
#
#
#  image.plot(1:(n_length), 1:(n_length),res_mat_cor_Calib,
#           col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#  legend("topleft",c("Period 1", "missing cor."),pch=c(NaN,6),
#       text.col=c("blue", "black"))
#  points(points_missingCalibx, points_missingCaliby, pch=6, cex=0.9,
#         col="gray0")
# 
#  #### before plotting, which value in Proj is not computed but computed in
#  #Calib period?
#  points_missingProj=which((is.na(res_mat_cor_Proj)
#                            & !is.na(res_mat_cor_Calib)),arr.ind=TRUE)
#  points_missingProjx=points_missingProj[,1]
#  points_missingProjy=points_missingProj[,2]
#
#
#  image.plot(1:(n_length), 1:(n_length),res_mat_cor_Proj,
#           col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#  legend("topleft",c("Period 2", "missing cor."), pch=c(NaN,6),
#       text.col=c("darkorange", "black"))
#  points(points_missingProjx, points_missingProjy, pch=6, cex=0.9,
#         col="gray0")
#  image.plot(1:(n_length), 1:(n_length),res_mat_pvalFisher,
#           col=col_pval_(seq(0,1, length=100)), xaxt="n", yaxt="n",zlim=c(0,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#  legend("topleft",c("Fisher-test pval"))
# 
#
#  points_pvalFisherx_10=which(res_mat_pvalFisher<=0.10,arr.ind=TRUE)[,1]
#  points_pvalFishery_10=which(res_mat_pvalFisher<=0.10,arr.ind=TRUE)[,2]
#  points(points_pvalFisherx_10, points_pvalFishery_10, pch=".", cex=3,
#         col="gray90")
#  points_pvalFisherx_05=which(res_mat_pvalFisher<=0.05,arr.ind=TRUE)[,1]
#  points_pvalFishery_05=which(res_mat_pvalFisher<=0.05,arr.ind=TRUE)[,2]
#  points(points_pvalFisherx_05, points_pvalFishery_05, pch=".", cex=3,
#         col="gray66")
#  points_pvalFisherx_01=which(res_mat_pvalFisher<=0.01,arr.ind=TRUE)[,1]
#  points_pvalFishery_01=which(res_mat_pvalFisher<=0.01,arr.ind=TRUE)[,2]
#  points(points_pvalFisherx_01, points_pvalFishery_01, pch=".", cex=3, col="gray0")
#  legend("topright", inset=.02, 
#           legend= c("pval10","pval5","pval1"), col=c("gray90", "gray66",
#                                                      "gray0"),
#         pch=c(20,20,20) , horiz=FALSE, cex=0.8)
#
#  diff_cor_Proj_Calib=res_mat_cor_Proj-res_mat_cor_Calib
#  diff_cor_Proj_Calib[which(res_mat_pvalFisher>=0.10)]<-NaN
#  image.plot(1:(n_length), 1:(n_length),diff_cor_Proj_Calib,
#           col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#  legend("topleft",c("Diff. cor"))
# 
#
#  return(pobs_data1_Calib)
#}
#
#
#
#compute_evolution_conditional_med<-function(data1, data2, col_="black",
#                                     col_pval_="black", name_plot="to def.", label_HW="HW", 
#                                sens=1,recomp_rank=FALSE){
#
#  ### compute Condcorr for calib and proj the half of the whole ts
#  n_length=length(data1)/2
# 
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  if(recomp_rank==FALSE){
#    pobs_data1_Calib=pobs_data1[1:n_length]
#    pobs_data2_Calib=pobs_data2[1:n_length]
#    pobs_data1_Proj=pobs_data1[(n_length+1):(2*n_length)]
#    pobs_data2_Proj=pobs_data2[(n_length+1):(2*n_length)]
#  }else{
#    pobs_data1_Calib=pobs(as.matrix(cbind(data1[1:n_length])))[,1]
#    pobs_data2_Calib=pobs(as.matrix(cbind(data2[1:n_length])))[,1]
#    pobs_data1_Proj=pobs(as.matrix(cbind(data1[(n_length+1):(2*n_length)])))[,1]
#    pobs_data2_Proj=pobs(as.matrix(cbind(data2[(n_length+1):(2*n_length)])))[,1]
#  }
#  
# 
#  plot(pobs_data1_Calib, pobs_data2_Calib, col= "blue", ylim=c(0,1),
#       xlim=c(0,1), ylab="", xlab="", main=name_plot)
#  points(pobs_data1_Proj, pobs_data2_Proj, col= "darkorange")
#  title(ylab="F(Drought)",xlab=paste0("F(",label_HW,")"),
#     cex.lab=1.2, line=2)
#  legend("topleft",pch=c(1,1),c("Period 1", "Period 2"),
#       col=c("blue","darkorange"))
#
#  ### Cond. median conditional to PR
#  compute_conditional_median(pobs_data1_Calib,pobs_data2_Calib,col_="black", lab_X="F(Drought)",
#                             lab_Y=paste0("F(",label_HW,")"))
#  compute_conditional_median(pobs_data1_Proj,pobs_data2_Proj,col_="black", lab_X="F(Drought)",
#                             lab_Y=paste0("F(",label_HW,")"))
#  
#  ### Cond. median conditional to TAS
#  compute_conditional_median(pobs_data2_Calib,pobs_data1_Calib,col_="black",
#                             lab_Y="F(Drought)", 
#                             lab_X=paste0("F(", label_HW,")"))
#  compute_conditional_median(pobs_data2_Proj,pobs_data1_Proj,col_="black",
#                             lab_Y="F(Drought)",
#                             lab_X=paste0("F(",label_HW,")"))
#  return(pobs_data1_Calib)
#}
#
#
#
#
#
#
#compute_evolution_edist<-function(data1, data2, col_="black",
#                                  col_pval_="black", name_plot="to def.", label_HW="HW", 
#                                sens=1,recomp_rank=FALSE){
#
#  ### compute Condcorr for calib and proj the half of the whole ts
#  n_length=length(data1)/2
# 
#  u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  res_mat_edist=res_mat_pval_edist =matrix(rep(NaN,n_length*n_length),ncol=n_length)
# 
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  if(recomp_rank==FALSE){
#    pobs_data1_Calib=pobs_data1[1:n_length]
#    pobs_data2_Calib=pobs_data2[1:n_length]
#    pobs_data1_Proj=pobs_data1[(n_length+1):(2*n_length)]
#    pobs_data2_Proj=pobs_data2[(n_length+1):(2*n_length)]
#  }else{
#    pobs_data1_Calib=pobs(as.matrix(cbind(data1[1:n_length])))[,1]
#    pobs_data2_Calib=pobs(as.matrix(cbind(data2[1:n_length])))[,1]
#    pobs_data1_Proj=pobs(as.matrix(cbind(data1[(n_length+1):(2*n_length)])))[,1]
#    pobs_data2_Proj=pobs(as.matrix(cbind(data2[(n_length+1):(2*n_length)])))[,1]
#  }
#  
#  ### Calib and Proj diff
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      if(sens==1){
#        id_Calib=which(pobs_data1_Calib>=u & pobs_data2_Calib>=v)
#        id_Proj=which(pobs_data1_Proj>=u & pobs_data2_Proj>=v)
#      }
#     if((length(id_Calib)<=5) | (length(id_Proj)<=5)){
#        next
#      }else{
#     ### subselection
#        pobs_data1_Calib_sub<-pobs_data1_Calib[id_Calib]#(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_Calib_sub<-pobs_data2_Calib[id_Calib]#(as.matrix(cbind(data2[id])))[,1]
#        pobs_data1_Proj_sub<-pobs_data1_Proj[id_Proj]#(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_Proj_sub<-pobs_data2_Proj[id_Proj]#(as.matrix(cbind(data2[id])))[,1]
#
#        ### Edist
#        x=cbind(pobs_data1_Calib_sub, pobs_data2_Calib_sub)
#        y=cbind(pobs_data1_Proj_sub, pobs_data2_Proj_sub)
#        nrow_x=nrow(x)
#        nrow_y=nrow(y)
#        tmp_edist=eqdist.etest(rbind(x,y), c(nrow_x, nrow_y), R=1000)
#      
#        res_mat_edist[count_u,count_v]=tmp_edist$statistic*(nrow_x+nrow_y)/(nrow_x*nrow_y)
#        res_mat_pval_edist[count_u, count_v]=tmp_edist$p.value
# 
#     }
#    }
#  }
# 
#  plot(pobs_data1_Calib, pobs_data2_Calib, col= "blue", ylim=c(0,1),
#       main=name_plot,
#       xlim=c(0,1), ylab="", xlab="")
#  points(pobs_data1_Proj, pobs_data2_Proj, col= "darkorange")
#  title(ylab="F(Drought)",xlab=paste0("F(",label_HW,")"),
#     cex.lab=1.2, line=2)
#  image.plot(1:(n_length), 1:(n_length),res_mat_edist,
#           col=col_[33:64], xaxt="n", yaxt="n", xlab="", ylab="",
#           zlim=c(0,max(res_mat_edist,na.rm=TRUE)))
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#
#  image.plot(1:(n_length), 1:(n_length),res_mat_pval_edist,
#           col=col_pval_(seq(0,1,length=100)), xaxt="n", yaxt="n",zlim=c(0,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)",cex.lab=1.2, line=2)
#  points_pval_edistx_10=which(res_mat_pval_edist<=0.10,arr.ind=TRUE)[,1]
#  points_pval_edisty_10=which(res_mat_pval_edist<=0.10,arr.ind=TRUE)[,2]
#  points(points_pval_edistx_10, points_pval_edisty_10, pch=".", cex=3,
#         col="gray90")
#  points_pval_edistx_05=which(res_mat_pval_edist<=0.05,arr.ind=TRUE)[,1]
#  points_pval_edisty_05=which(res_mat_pval_edist<=0.05,arr.ind=TRUE)[,2]
#  points(points_pval_edistx_05, points_pval_edisty_05, pch=".", cex=3,
#         col="gray66")
#  points_pval_edistx_01=which(res_mat_pval_edist<=0.01,arr.ind=TRUE)[,1]
#  points_pval_edisty_01=which(res_mat_pval_edist<=0.01,arr.ind=TRUE)[,2]
#  points(points_pval_edistx_01, points_pval_edisty_01, pch=".", cex=3, col="gray0")
#  legend("topright", inset=.02, 
#           legend= c("pval10","pval5","pval1"), col=c("gray90", "gray66",
#                                                      "gray0"),
#         pch=c(20,20,20) , horiz=FALSE, cex=0.8)
#
#  return(pobs_data1_Calib)
#}
##### Running Evolution
#compute_running_edist<-function(data1, data2, sliding_window=30,
#                                recomp_rank=FALSE){
#  n_length=length(data1)
#  res_running_edist_v1= res_running_edistpval_v1=c()
#  res_running_cor_v1= res_running_corpval_v1=res_running_corFisherpval_v1=c()
#
#  p.val_fisher = function(r1,r2,n1,n2){
#    fisher= ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
#    p.value = (2*(1-pnorm(abs(fisher))))
#    return(p.value)
#  }
#
#  if(recomp_rank==TRUE){
#  ### v1: compute edist on ranks by recomputing the ranks
#    x=cbind(pobs(as.matrix(cbind(data1[(1:sliding_window)])))[,1], #rank(data1[(1:sliding_window)])/sliding_window,
#            pobs(as.matrix(cbind(data2[(1:sliding_window)])))[,1]) #rank(data2[(1:sliding_window)])/sliding_window)
#    cor_x=cor(x[,1],x[,2],method="spearman")
#    for(i in 0:(n_length-sliding_window)){
#      y=cbind(pobs(as.matrix(cbind(data1[(1:sliding_window)+i])))[,1], 
#            pobs(as.matrix(cbind(data2[(1:sliding_window)+i])))[,1]) 
#      #cbind(rank(data1[(1:sliding_window)+i])/sliding_window,
#        #      rank(data2[(1:sliding_window)+i])/sliding_window)  
#      ### Edist
#      tmp_edist=eqdist.etest(rbind(x,y), c(nrow(x), nrow(y)), R=1000)
#      res_running_edist_v1=c(res_running_edist_v1,
#                             tmp_edist$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_v1=c(res_running_edistpval_v1,tmp_edist$p.value)
#      ### Cor
#      cor_y=cor(y[,1], y[,2] ,method="spearman")
#      res_running_cor_v1=c(res_running_cor_v1,
#                             cor(y[,1], y[,2] ,method="spearman"))
#      res_running_corpval_v1=c(res_running_corpval_v1,
#                               cor.test(x,y,method="spearman")$p.value)
#      res_running_corFisherpval_v1=c(res_running_corFisherpval_v1, p.val_fisher(cor_x,cor_y,sliding_window,sliding_window))
#    }
#  }else{
#    ### v2: compute edist on ranks without recomputing the ranks
#    rank_data1=pobs(as.matrix(cbind(data1)))[,1] #rank(data1)/n_length
#    rank_data2=pobs(as.matrix(cbind(data2)))[,1] #r]ank(data2)/n_length
#    x=cbind(rank_data1[(1:sliding_window)], rank_data2[(1:sliding_window)])
#    cor_x=cor(x[,1],x[,2],method="spearman")
#    for(i in 0:(n_length-sliding_window)){
#      y=cbind(rank_data1[(1:sliding_window)+i], rank_data2[(1:sliding_window)+i])  
#      
#      ### Edist
#      tmp_edist=eqdist.etest(rbind(x,y), c(nrow(x), nrow(y)), R=1000)
#      res_running_edist_v1=c(res_running_edist_v1,
#                             tmp_edist$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_v1=c(res_running_edistpval_v1,tmp_edist$p.value)
#  
#      ### Cor
#      cor_y=cor(y[,1], y[,2] ,method="spearman")
#      res_running_cor_v1=c(res_running_cor_v1,
#                             cor(y[,1], y[,2] ,method="spearman"))
#      res_running_corpval_v1=c(res_running_corpval_v1,
#                               cor.test(x,y,method="spearman")$p.value)
#      res_running_corFisherpval_v1=c(res_running_corFisherpval_v1, p.val_fisher(cor_x,cor_y,sliding_window,sliding_window))
#    }
#  }
#  color_vec=c(rep("blue",sliding_window),rep("black",10),
#              rep("darkorange",sliding_window))
#  plot(pobs(as.matrix(cbind(data1)))[,1],pobs(as.matrix(cbind(data2)))[,1], col=color_vec)
#  plot(res_running_edist_v1, ylab="Running Edist", xlab="Sliding window",
#       ylim=c(0,0.2))
#  abline(h=0,type='l')
#  plot(res_running_edistpval_v1, ylab= "pval Edist", xlab= "Sliding window",ylim=c(0,1))
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  plot(res_running_cor_v1, ylab="Running Cor", xlab="Sliding window", type='l',
#       ylim=c(-1,1))
#  abline(h=0,type='l')
#  plot(res_running_corpval_v1, ylab= "pval Cor", xlab= "Sliding window", ylim=c(0,1))
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  plot(res_running_corFisherpval_v1, ylab= "Fisher pval Cor", xlab= "Sliding window", ylim=c(0,1))
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  return(list=c(edist=res_running_edist_v1, pval=res_running_edistpval_v1))
#}
#
#
#
#compute_running_edist_between_RefMod<-function(Refdata1, Refdata2, Moddata1, Moddata2, sliding_window=30,
#                                               recomp_rank=FALSE){
#  n_length=length(Refdata1)
#  res_running_edist_Ref= res_running_edistpval_Ref= res_running_edist_Mod= res_running_edistpval_Mod=c()
#  res_running_cor_Ref= res_running_corpval_Ref= res_running_cor_Mod=c()
#  res_running_corpval_Mod=res_running_corFisherpval_Ref=res_running_corFisherpval_Mod=c()
# 
#  p.val_fisher = function(r1,r2,n1,n2){
#    fisher= ((0.5*log((1+r1)/(1-r1)))-(0.5*log((1+r2)/(1-r2))))/((1/(n1-3))+(1/(n2-3)))^0.5
#    p.value = (2*(1-pnorm(abs(fisher))))
#    return(p.value)
#  }
#
# 
#  if(recomp_rank==TRUE){
#    x=cbind(rank(Refdata1[(1:sliding_window)])/sliding_window,
#            rank(Refdata2[(1:sliding_window)])/sliding_window)
#    ### Ref: compute edist on ranks by recomputing the ranks
#    for(i in 0:(n_length-sliding_window)){
#      yRef=cbind(rank(Refdata1[(1:sliding_window)+i])/sliding_window,
#                 rank(Refdata2[(1:sliding_window)+i])/sliding_window)
#      yMod=cbind(rank(Moddata1[(1:sliding_window)+i])/sliding_window,
#              rank(Moddata2[(1:sliding_window)+i])/sliding_window)
#      ### Edist
#      tmp_edist_Ref=eqdist.etest(rbind(x,yRef), c(nrow(x), nrow(yRef)), R=1000)
#      res_running_edist_Ref=c(res_running_edist_Ref,
#                              tmp_edist_Ref$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_Ref=c(res_running_edistpval_Ref,tmp_edist_Ref$p.value)
#      
#      tmp_edist_Mod=eqdist.etest(rbind(x,yMod), c(nrow(x), nrow(yMod)), R=1000)
#      res_running_edist_Mod=c(res_running_edist_Mod,
#                              tmp_edist_Mod$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_Mod=c(res_running_edistpval_Mod,tmp_edist_Mod$p.value)
#      
#      ### Cor
#      cor_x=cor(x[,1], x[,2], method="spearman")
#      cor_yRef= cor(yRef[,1], yRef[,2], method = "spearman")
#      cor_yMod= cor(yMod[,1], yMod[,2], method= "spearman")
#      res_running_cor_Ref=c(res_running_cor_Ref, cor_yRef)
#      res_running_cor_Mod=c(res_running_cor_Mod, cor_yMod)
#
#      res_running_corFisherpval_Ref=c(res_running_corFisherpval_Ref,p.val_fisher(cor_x, cor_yRef, sliding_window,
#                                                 sliding_window))
#      res_running_corFisherpval_Mod=c(res_running_corFisherpval_Mod,p.val_fisher(cor_x, cor_yMod, sliding_window,
#                                                 sliding_window))
#    }
#  }else{
#    ### v2: compute edist on ranks without recomputing the ranks
#    rank_Refdata1=rank(Refdata1)/n_length
#    rank_Refdata2=rank(Refdata2)/n_length
#    rank_Moddata1=rank(Moddata1)/n_length
#    rank_Moddata2=rank(Moddata2)/n_length
#    
#    x=cbind(rank_Refdata1[(1:sliding_window)], rank_Refdata2[(1:sliding_window)])
#    for(i in 0:(n_length-sliding_window)){
#      yRef=cbind(rank_Refdata1[(1:sliding_window)+i], rank_Refdata2[(1:sliding_window)+i])
#      ### Edist
#      tmp_edist_Ref=eqdist.etest(rbind(x,yRef), c(nrow(x), nrow(yRef)), R=1000)
#      res_running_edist_Ref=c(res_running_edist_Ref,
#                              tmp_edist_Ref$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_Ref=c(res_running_edistpval_Ref,tmp_edist_Ref$p.value)
#      
#      yMod=cbind(rank_Moddata1[(1:sliding_window)+i], rank_Moddata2[(1:sliding_window)+i])
#      ### Edist
#      tmp_edist_Mod=eqdist.etest(rbind(x,yMod), c(nrow(x), nrow(yMod)), R=1000)
#      res_running_edist_Mod=c(res_running_edist_Mod,
#                              tmp_edist_Mod$statistic*(sliding_window+sliding_window)/(sliding_window*sliding_window))
#      res_running_edistpval_Mod=c(res_running_edistpval_Mod,tmp_edist_Mod$p.value)      
#
#      ### Cor
#      cor_x=cor(x[,1], x[,2], method="spearman")
#      cor_yRef= cor(yRef[,1], yRef[,2], method = "spearman")
#      cor_yMod= cor(yMod[,1], yMod[,2], method= "spearman")
#      res_running_cor_Ref=c(res_running_cor_Ref, cor_yRef)
#      res_running_cor_Mod=c(res_running_cor_Mod, cor_yMod)
#
#      res_running_corFisherpval_Ref=c(res_running_corFisherpval_Ref,p.val_fisher(cor_x, cor_yRef, sliding_window,
#                                                 sliding_window))
#      res_running_corFisherpval_Mod=c(res_running_corFisherpval_Mod,p.val_fisher(cor_x, cor_yMod, sliding_window,
#                                                 sliding_window))
# 
#
#    }
#  }
#  
#  
#  ###to del
#  rank_Refdata1=rank(Refdata1)/n_length
#  rank_Refdata2=rank(Refdata2)/n_length
#  rank_Moddata1=rank(Moddata1)/n_length
#  rank_Moddata2=rank(Moddata2)/n_length
#  x=cbind(rank_Refdata1, rank_Refdata2)
#  y=cbind(rank_Moddata1, rank_Moddata2)
#  glob_pval=eqdist.etest(rbind(x,y), c(nrow(x), nrow(y)), R=1000)$p.value
#  ###end to del
#  color_vec=c(rep("royalblue",30),rep("black",10), rep("darkorange",30))
#  plot(rank(Refdata1)/n_length,rank(Refdata2)/n_length, col=color_vec,main="ERA5")
#  plot(rank(Moddata1)/n_length,rank(Moddata2)/n_length, col=color_vec, main="CNRMCM6")
#  
#  plot(res_running_edist_Ref, ylab="Running Edist", xlab="Sliding window",
#       type='l',col="blue")
#  lines(res_running_edist_Mod, ylab="Running Edist", xlab="Sliding window",
#       type='l',col="red")
#  plot(res_running_edistpval_Ref, ylab= "pval Edist", xlab= "Sliding window",ylim=c(0,1),
#       main=paste0("global pval: ",round(glob_pval,5)),col="blue")
#  lines(res_running_edistpval_Mod, ylab="Running Edist", xlab="Sliding window",
#        type='l',col="red")
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  
#  plot(res_running_cor_Ref, col="blue", ylim=c(-1,1))
#  lines(res_running_cor_Mod, col="red")
#
#  plot(res_running_corFisherpval_Ref, ylab= "Fisher pval Cor", xlab= "Sliding
#       window", ylim=c(0,1), col="blue")
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  lines(res_running_corFisherpval_Mod, ylab= "Fisher pval Cor", xlab= "Sliding
#        window", ylim=c(0,1), col="red")
#  return(list=c(edist=res_running_edist_Ref, pval=res_running_edistpval_Ref))
#}
#
#
#
#compute_and_plot_probs_bivar<-function(data1, data2, recomp_rank=FALSE, pobs_=0.75,
#                                       name_plot="to_define",label_HW="HW",
#                                       sliding_window=30){ #data1 is tas, data2 is pr
#  n_length=length(data1)/2
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  if(recomp_rank==FALSE){
#    pobs_data1_Calib=pobs_data1[1:n_length]
#    pobs_data2_Calib=pobs_data2[1:n_length]
#    pobs_data1_Proj=pobs_data1[(n_length+1):(2*n_length)]
#    pobs_data2_Proj=pobs_data2[(n_length+1):(2*n_length)]
#  }else{
#    pobs_data1_Calib=pobs(as.matrix(cbind(data1[1:n_length])))[,1]
#    pobs_data2_Calib=pobs(as.matrix(cbind(data2[1:n_length])))[,1]
#    pobs_data1_Proj=pobs(as.matrix(cbind(data1[(n_length+1):(2*n_length)])))[,1]
#    pobs_data2_Proj=pobs(as.matrix(cbind(data2[(n_length+1):(2*n_length)])))[,1]
#  }
#  cross_pobs_Calib<-c(pobs_data1_Calib*pobs_data2_Calib) 
#  cross_pobs_Proj<-c(pobs_data1_Proj*pobs_data2_Proj)#abs((-1)*(1-pobs_data1)*(1-pobs_data2)) 
# 
#### Recompute
#  pobs_data1<-c(pobs_data1_Calib, pobs_data1_Proj)#pobs_data1*pobs_data2#abs((-1)*(1-pobs_data1)*(1-pobs_data2)) 
#  pobs_data2<-c(pobs_data2_Calib, pobs_data2_Proj)#pobs_data1*pobs_data2#abs((-1)*(1-pobs_data1)*(1-pobs_data2)) 
#  cross_pobs<-c(cross_pobs_Calib, cross_pobs_Proj)#pobs_data1*pobs_data2#abs((-1)*(1-pobs_data1)*(1-pobs_data2)) 
#  
#  count_51_86=sum(pobs_data1_Calib>=pobs_)
#  count_87_20=sum(pobs_data1_Proj>=pobs_)
#
#  #### univ. CvM
#  CvM_pval_data1=cvm(pobs_data1_Calib,pobs_data1_Proj)$pval
#  CvM_pval_data2=cvm(pobs_data2_Calib,pobs_data2_Proj)$pval
#  CvM_pval_cross=cvm(cross_pobs_Calib,cross_pobs_Proj)$pval
# 
#  plot(c(pobs_data1_Calib, pobs_data1_Proj), ylab=paste0("prob.", label_HW), xlab="Years", xaxt="n",
#       type='l',col="darkorange", ylim=c(0,1), main=paste0("Period 1: ", count_51_86, ", Period2: ", count_87_20,", CvM: ", round(CvM_pval_data1,3)))
#  axis(1,seq.int(1,70,5),seq.int(1951,2020,5))
#  abline(v=36, lty=2)
#  abline(h=pobs_, col="gray0")
#  xx=(pobs_data1>=pobs_)
#  xx[!xx]=NaN
#  points(xx, col="darkorange")
#
#  count_51_86=sum(pobs_data2_Calib>=pobs_)
#  count_87_20=sum(pobs_data2_Proj>=pobs_)
#  plot(c(pobs_data2_Calib, pobs_data2_Proj), ylab="prob. Drought", xlab="Years", xaxt= "n",
#       type='l',col="deepskyblue", ylim=c(0,1), main=paste0("Period 1: ", count_51_86, ", Period 2: ", count_87_20,", CvM: ", round(CvM_pval_data2,3)))
#
#  axis(1,seq.int(1,70,5),seq.int(1951,2020,5))
#  abline(v=36, lty=2)
#  abline(h=pobs_, col="gray0")
#  xx=(pobs_data2>=pobs_)
#  xx[!xx]=NaN
#  points(xx, col="deepskyblue")
#
#  count_51_86=sum(cross_pobs_Calib>=(pobs_*pobs_))
#  count_87_20=sum(cross_pobs_Proj>=(pobs_*pobs_))
#  plot(c(cross_pobs_Calib, cross_pobs_Proj), ylab="Cross prob.", xlab="Years", xaxt= "n",
#       type='l',col="purple", ylim=c(0,1), main=paste0("Period 1: ", count_51_86, ", Period 2: ", count_87_20,", CvM: ", round(CvM_pval_cross,3)))
#  axis(1,seq.int(1,70,5),seq.int(1951,2020,5))
#  abline(v=36, lty=2)
#  abline(h=(pobs_*pobs_), col="gray0")
#  xx=(cross_pobs>=(pobs_*pobs_))
#  xx[!xx]=NaN
#  points(xx, col="purple")
#
##  ###sliding_window
#  CvM_pval_sliding_data1=c()
#  CvM_pval_sliding_data2=c()
#  CvM_pval_sliding_cross=c()
#  for(i in 0:(2*n_length-sliding_window)){
#    x_data1=pobs_data1[1:sliding_window]
#    y_data1=pobs_data1[(1:sliding_window)+i] 
#    CvM_pval_sliding_data1=c(CvM_pval_sliding_data1, cvm(x_data1,y_data1)$pval)
#    x_data2=pobs_data2[1:sliding_window]
#    y_data2=pobs_data2[(1:sliding_window)+i] 
#    CvM_pval_sliding_data2=c(CvM_pval_sliding_data2, cvm(x_data2,y_data2)$pval)
#    x_cross_pobs=cross_pobs[1:sliding_window]
#    y_cross_pobs=cross_pobs[(1:sliding_window)+i] 
#    CvM_pval_sliding_cross=c(CvM_pval_sliding_cross,
#                             cvm(x_cross_pobs,y_cross_pobs)$pval)
#  }
#  plot(CvM_pval_sliding_data1, ylab= "CvM pval", xlab= "Sliding window",
#       ylim=c(0,1), col="darkorange", type='l')
#  lines(CvM_pval_sliding_data2, col="deepskyblue")
#  lines(CvM_pval_sliding_cross, col="purple")
#  abline(h=0.10, col="gray90")
#  abline(h=0.05, col="gray66")
#  abline(h=0.01, col="gray0")
#  return(n_length)
#}
#

#### Heatwave durant l'ete: Max des Moyennes de T2mean sur 5j sliding window
#compute_max_of_running_mean<-function(seasonal_data,nb_days_by_season=92,bandwidth=5){ #(take into account the time gap btw years)
#  nb_year=length(seasonal_data)/nb_days_by_season
#  # res=rep(NaN,(nb_days_by_season-bandwidth+1)*nb_year)
#  res_running_mean=c()
#  for(k in 1:nb_year){
#    res_running_mean=c(res_running_mean,running_mean(seasonal_data[(1+nb_days_by_season*(k-1)):(nb_days_by_season*k)],bandwidth))
#  }
#  mat_res=matrix(res_running_mean,byrow=FALSE, ncol = nb_year)
#  res_max=apply(mat_res,2,max)
#  #as we divide running_mean by year, length of res is equal to nb_days_by_season-bandwidth+1)*nb_year
#  return(res_max) 
#}
#
#
#### nb jours chauds par ete
##### compute 90th percentile for all the summers and count exceedances by summer
##From Mueller Sonia 2012
#compute_nb_hot_days_alla_Seneveritnae<-function(dt_data, nb_years,
#                                                sliding_days=5,
#                                                nb_days_by_season=92,
#                                     probs_quant=0.90){
#  ### Rearranging data par ete en colonne 
#  mean_spavg_by_summer=matrix(dt_data, byrow=FALSE, ncol = nb_years)
#  
#  quant95=rep(NaN,nb_days_by_season)#quantile(mean_spavg_by_summer,probs=probs_quant)
#  
#  quant95[1]=quantile(mean_spavg_by_summer[1:3,],probs=probs_quant)
#  quant95[2]=quantile(mean_spavg_by_summer[1:4,],probs=probs_quant)
#  
#  for(i in 3:(nb_days_by_season-2)){
#    quant95[i]=quantile(mean_spavg_by_summer[(i-2):(i+2),],probs=probs_quant)
#  }
#  quant95[nb_days_by_season-1]=quantile(mean_spavg_by_summer[(nb_days_by_season-3):nb_days_by_season,],probs=probs_quant)
#  quant95[nb_days_by_season]=quantile(mean_spavg_by_summer[(nb_days_by_season-2):nb_days_by_season,],probs=probs_quant)
# 
#  boolean_mat= matrix(NaN, ncol=nb_years, nrow=nb_days_by_season)# (mean_spavg_by_summer>quant95)
#  
#  for(i in 1:nb_years){
#    boolean_mat[,i]=mean_spavg_by_summer[,i]>=quant95
#  }
#
##  count_sequence<-function(bool,thresh_days=4){
##    count_seq=rle(bool)
##    tmp=count_seq$lengths[which(count_seq$values==TRUE)]
##    return(sum(tmp[which(tmp>=thresh_days)]))
##  }
#
#  res=apply(boolean_mat, 2, sum)
#  #res_persistence_cumule[which(is.nan(res_persistence_cumule))]=0
#  return(res)
#}
#
#
##### persistence cumule: nb of days conditionally to a Heatwave
### example of definition: HW est 4 jours consecutifs > quant95
#### on additionne par été le nb de jours constituant des segments d'au moins 4j consec
#### > quant95
#### From Lorenz and Sonia
#compute_persistence_cumule_alla_Lorenz<-function(dt_data, nb_years,
#                                                 sliding_days=5,
#                                                 nb_days_by_season=92,
#                                     probs_quant=0.90){
#  ### Rearranging data par ete en colonne 
#  mean_spavg_by_summer=matrix(dt_data, byrow=FALSE, ncol = nb_years)
#  
#  quant95=rep(NaN,nb_days_by_season)#quantile(mean_spavg_by_summer,probs=probs_quant)
#  
#  quant95[1]=quantile(mean_spavg_by_summer[1:3,],probs=probs_quant)
#  quant95[2]=quantile(mean_spavg_by_summer[1:4,],probs=probs_quant)
#  
#  for(i in 3:(nb_days_by_season-2)){
#    quant95[i]=quantile(mean_spavg_by_summer[(i-2):(i+2),],probs=probs_quant)
#  }
#  quant95[nb_days_by_season-1]=quantile(mean_spavg_by_summer[(nb_days_by_season-3):nb_days_by_season,],probs=probs_quant)
#  quant95[nb_days_by_season]=quantile(mean_spavg_by_summer[(nb_days_by_season-2):nb_days_by_season,],probs=probs_quant)
# 
#  boolean_mat= matrix(NaN, ncol=nb_years, nrow=nb_days_by_season)# (mean_spavg_by_summer>quant95)
#  
#  for(i in 1:nb_years){
#    boolean_mat[,i]=mean_spavg_by_summer[,i]>=quant95
#  }
#
#  count_sequence<-function(bool,thresh_days=2){
#    count_seq=rle(bool)
#    tmp=count_seq$lengths[which(count_seq$values==TRUE)]
#    return(mean(tmp[which(tmp>=thresh_days)]))
#  }
#
#  res_persistence_cumule=apply(boolean_mat, 2,function(x)count_sequence(x))
#  res_persistence_cumule[which(is.nan(res_persistence_cumule))]=0
#  return(res_persistence_cumule)
#}
#
#
#compute_acf1<-function(dt_data, nb_years){
#  ### Rearranging data par ete en colonne 
#  mean_spavg_by_summer=matrix(dt_data, byrow=FALSE, ncol = nb_years)
#  
#  res_acf1=apply(mean_spavg_by_summer, 2, function(x)return(acf(x, plot=FALSE)$acf[2]))
#  return(res_acf1)
#}
#
### Compute PR standardisation
#### On standardise par site. 
#### Pour chaque site, on calcule la moyenne printanière sur la totalite des annees.
#### On enleve cette moyenne a tous les jours du printemps
#### Ensuite,  On agrège spatialement et temporellement les valeurs pour chaque printemp en moyennant le tout.
#### outputs: 41 valeurs de la moyenne de pluie standardise spatialement.
#scale_pr<-function(pr_data,Ind_season,nb_years){
#  nb_lat=dim(pr_data)[1]
#  nb_lon=dim(pr_data)[2]
#  array_of_spring_means=array(rep(apply(pr_data[,,Ind_season],c(1,2),mean,na.rm=T),length(Ind_season)),dim=c(nb_lat,nb_lon,length(Ind_season)))
#  daily_scaled_pr<-pr_data[,,Ind_season]-array_of_spring_means
#  
#  tmp_daily_scaled_pr <- t(matrix(daily_scaled_pr, nrow = nrow(daily_scaled_pr) * ncol(daily_scaled_pr)))
#  
#  spavg_scaled_pr <- apply(tmp_daily_scaled_pr, 1, mean, na.rm=T) #Spatial average for the season
#  
#  tmp_spavg_scaled_pr <- matrix(spavg_scaled_pr,
#                                byrow=FALSE, ncol = nb_years)
#  #Temporal average by season
#  scaled_mean_pr=apply(tmp_spavg_scaled_pr,2,mean)
#  return(scaled_mean_pr)
#}
#
#scale_bypoint_pr<-function(dt_matrix,nb_years, Ind_info=NaN, show_plot=FALSE){#dt_matrix matrix of nb_days_spring x point
#  res_scaled=matrix(NaN, nrow=dim(dt_matrix)[1], ncol=dim(dt_matrix)[2])
#  res_scaled_mean_pr= matrix(NaN,nrow=nb_years, ncol= dim(dt_matrix)[2])
#  point_max= 1:dim(dt_matrix)[2]
#  for(point in point_max){
#    #i=Ind_info[point,1]
#    #j=Ind_info[point,2]
#    spring_mean=mean(dt_matrix[,point])
#    res_scaled[,point]<-dt_matrix[,point]-spring_mean
#  
#    tmp_avgbyspring_scaled_pr <- matrix(res_scaled[,point],byrow=FALSE, ncol = nb_years)
#    ### compute mean by season/year
#    res_scaled_mean_pr[,point]=apply(tmp_avgbyspring_scaled_pr,2,mean)
#  }
#  return(res_scaled_mean_pr)
#}
#
#
#
#
#### Quel seuil u,v maximise la correlation de Spearman entre HW et scaled_pr?
#compute_matrix_Condcorr<-function(data1, data2, col_="black", 
#                                  sens=1, 
#                                  name_plot="to_define",label_HW="HW",
#                                  boot_array_prop=NaN){ #data1 is tas, data2 is pr
#  n_length=length(data1)
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  u_thresh=v_thresh=sort(pobs_data1)#problem of not equality with seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  res_mat_cor=res_mat_prop=res_mat_pval=res_mat_airetheorique=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#  
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      if(sens==1){id=which(pobs_data1>=u & pobs_data2>=v)}
##      if(sens==2){id=which(pobs_data1>=u & pobs_data2<=v)}
##      if(sens==3){id=which(pobs_data1<=u & pobs_data2<=v)}
##      if(sens==4){id=which(pobs_data1<=u & pobs_data2>=v)}
#      
#      if(length(id)<=5){
#        next
#      }else{
#     ### subselection
#        pobs_data1_sub<-pobs(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_sub<-pobs(as.matrix(cbind(data2[id])))[,1]
#        res_mat_cor[count_u, count_v]=cor(pobs_data1_sub,pobs_data2_sub,method="spearman")
#        res_mat_prop[count_u, count_v]=length(id)/n_length
#        res_mat_airetheorique[count_u, count_v]=(1-u)*(1-v)
#        res_mat_pval[count_u, count_v]=cor.test(pobs_data1_sub,pobs_data2_sub,method="spearman")$p.value
#      }
#    }
#  }
#  
#  ### ATTENTION
#  res_mat_cor[which(res_mat_pval>0.10, arr.ind=TRUE)]<-NaN
#  ### END ATTENTION
#  
#  plot(pobs_data1,pobs_data2, main=name_plot,ylim=c(0,1), xlim=c(0,1),
#     col=c(rep("black",n_length/2),rep("black",n_length/2)), xlab="", ylab="")
#  title(ylab="F(Drought)",xlab=paste0("F(",label_HW,")"),
#     cex.lab=1.2, line=2)
#  image.plot(1:n_length, 1:n_length,res_mat_cor,
#           col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#  axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#  axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#  title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)", cex.lab=1.2, line=2)
#  legend("topright", inset=.02, 
#         legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                    "gray0", "green", "red"),
#       pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.95)
#  points_pvalx_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,1]
#  points_pvaly_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,2]
#  points(points_pvalx_10, points_pvaly_10, pch=".", cex=3, col="gray90")
#  points_pvalx_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,1]
#  points_pvaly_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,2]
#  points(points_pvalx_5, points_pvaly_5, pch=".", cex=3, col="gray66")
#  points_pvalx_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,1]
#  points_pvaly_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,2]
#  points(points_pvalx_1, points_pvaly_1, pch=".", cex=3, col="gray0")
#  points_minpvalx=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,1]
#  points_minpvaly=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,2]
#  points(points_minpvalx, points_minpvaly, pch=0, cex=1, col="green")
#  
#  IND_pval_5=which(res_mat_pval<=0.05,arr.ind=TRUE)
#  if(!(length(IND_pval_5)==0)){
#    maxcor_5=res_mat_cor[IND_pval_5][which.max(abs(res_mat_cor)[IND_pval_5])]
#    points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
#    points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#    points(points_maxcorx, points_maxcory, pch=0, cex=1, col="red")
#  }
#  
#  ### Cond. median conditional to PR
#  compute_conditional_median(pobs_data1,pobs_data2,col_="black", lab_X="F(Drought)",
#                             lab_Y=paste0("F(",label_HW,")"))
#  
#  ### Cond. median conditional to TAS
#  compute_conditional_median(pobs_data2,pobs_data1,col_="black", lab_X=paste0("F(",
#                                                                      label_HW,
#                                                                      ")"),
#                                                                      lab_Y="F(Drought)")
#  if(length(IND_pval_5)==0){
#   thresh_tas=NaN
#   thresh_pr=NaN
#   maxcor_5=NaN
#  }else{
#   thresh_tas=min(u_thresh[points_maxcorx])
#   thresh_pr=min(v_thresh[points_maxcory])
#   print(thresh_tas) ###TAS thresh
#   print(thresh_pr) ### PR thresh
#   print(maxcor_5)
#  }
#  return(list(res_mat_cor=res_mat_cor, thresh_tas=thresh_tas,
#              thresh_pr=thresh_pr, max_cor=maxcor_5))
#}
#
#
#### essai
#
#FDR_compute_matrix_Condcorr<-function(data1, data2, col_="black", 
#                                  sens=1, 
#                                  name_plot="to_define",label_HW="HW",
#                                  boot_array_prop=NaN){ #data1 is tas, data2 is pr
#  n_length=length(data1)
#  pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
#  u_thresh=v_thresh=sort(pobs_data1)#seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  res_mat_cor=res_mat_prop=res_mat_pval=res_mat_airetheorique=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#  
#  list_id=list(c("init"))
#  count_u=0
#  for(u in u_thresh){
#    count_u=count_u+1
#    count_v=0
#    for(v in v_thresh){
#      count_v=count_v+1
#      if(sens==1){
#        id=which(pobs_data1>=u & pobs_data2>=v)
#        ### TO AVOID RECOMPUTING CORR.
#        check_id_in_list_id=sum(Reduce("+",lapply(list_id,function(x)identical(x,id))))
#      }
#      if(length(id)<=5 | check_id_in_list_id==1){ ### TO AVOID RECOMPUTING CORR.
#        next
#      }else{
#        list_id[[length(list_id)+1]]=id
#        ### subselection
#        pobs_data1_sub<-pobs(as.matrix(cbind(data1[id])))[,1]
#        pobs_data2_sub<-pobs(as.matrix(cbind(data2[id])))[,1]
#        res_mat_cor[count_u, count_v]=cor(pobs_data1_sub,pobs_data2_sub,method="spearman")
#        res_mat_prop[count_u, count_v]=length(id)/n_length
#        res_mat_airetheorique[count_u, count_v]=(1-u)*(1-v)
#        res_mat_pval[count_u, count_v]=cor.test(pobs_data1_sub,pobs_data2_sub,method="spearman")$p.value
#      }
#    }
#  }
#  
#  ### ATTENTION
#  # res_mat_cor[which(res_mat_pval>0.10, arr.ind=TRUE)]<-NaN
#  ### END ATTENTION
#  
#  plot(pobs_data1,pobs_data2, main=name_plot,ylim=c(0,1), xlim=c(0,1),
#       col=c(rep("black",n_length/2),rep("black",n_length/2)), xlab="", ylab="")
#  title(ylab="F(Drought)",xlab=paste0("F(",label_HW,")"),
#        cex.lab=1.2, line=2)
#  for(k in 1:2){
#    if(k==2){res_mat_cor[which(res_mat_pval>0.10,arr.ind=TRUE)]<-NaN}
#    image.plot(1:n_length, 1:n_length,res_mat_cor,
#               col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="", 
#               main=paste0("nb_cor: ", length(which(!is.na(res_mat_cor))), ", nb_pval<=0.10: ", length(which((res_mat_pval<=0.1)))))
#    axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#    axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#    title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)", cex.lab=1.2, line=2)
#    legend("topright", inset=.02, 
#           legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                                           "gray0", "green", "red"),
#           pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.95)
#    points_pvalx_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,1]
#    points_pvaly_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,2]
#    points(points_pvalx_10, points_pvaly_10, pch=".", cex=3, col="gray90")
#    points_pvalx_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,1]
#    points_pvaly_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,2]
#    points(points_pvalx_5, points_pvaly_5, pch=".", cex=3, col="gray66")
#    points_pvalx_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,1]
#    points_pvaly_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,2]
#    points(points_pvalx_1, points_pvaly_1, pch=".", cex=3, col="gray0")
#    points_minpvalx=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,1]
#    points_minpvaly=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,2]
#    points(points_minpvalx, points_minpvaly, pch=0, cex=1, col="green")
#    
#    IND_pval_5=which(res_mat_pval<=0.05,arr.ind=TRUE)
#    if(!(length(IND_pval_5)==0)){
#      maxcor_5=res_mat_cor[IND_pval_5][which.max(abs(res_mat_cor)[IND_pval_5])]
#      points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
#      points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#      points(points_maxcorx, points_maxcory, pch=0, cex=1, col="red")
#    }
#  }
#  
#  compute_mat_pval_FDR<-function(mat_pval, mat_cor, alpha_FDR=c(0.2)){
#    IND_pval=which(!(is.na(mat_pval)),arr.ind=TRUE)
#    vec_pval=c()
#    for(k in 1:dim(IND_pval)[1]){
#      i=IND_pval[k,1]
#      j=IND_pval[k,2]
#      vec_pval=c(vec_pval,mat_pval[i,j])
#    }
#    sort_pval=sort(vec_pval)
#    length_pval=length(sort_pval)
#    
#    plot(sort_pval)
#    col_=c("red", "darkorange", "orange", "yellow")
#    for(i in 1:length(alpha_FDR)){
#      vec_comp=(1:length_pval)/length_pval*alpha_FDR[i]
#      id_pval_star=max(which(vec_comp-sort_pval>=0))
#      pval_star=sort_pval[id_pval_star]
#      lines(vec_comp, col=col_[i])
#      abline(h=pval_star, col=col_[i])
#      
#      if(!is.na(pval_star)){
#        mat_cor[which(mat_pval>pval_star,arr.ind=TRUE)]<-NaN
#        image.plot(1:n_length, 1:n_length,mat_cor,
#                   col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(64)), xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="", 
#                   main=paste0("nb_cor: ", length(which(!is.na(mat_cor))), "FDR:", alpha_FDR[i],", pvalFDR: ", round(pval_star,3)))
#      }
#
#    }
#
#  }
#  compute_mat_pval_FDR(res_mat_pval,res_mat_cor)
#  
#  # ### Cond. median conditional to PR
#  # compute_conditional_median(pobs_data1,pobs_data2,col_="black", lab_X="F(Drought)",
#  #                            lab_Y=paste0("F(",label_HW,")"))
#  # 
#  # ### Cond. median conditional to TAS
#  # compute_conditional_median(pobs_data2,pobs_data1,col_="black", lab_X=paste0("F(",
#  #                                                                             label_HW,
#  #                                                                             ")"),
#  #                            lab_Y="F(Drought)")
#  if(length(IND_pval_5)==0){
#    thresh_tas=NaN
#    thresh_pr=NaN
#    maxcor_5=NaN
#  }else{
#    thresh_tas=min(u_thresh[points_maxcorx])
#    thresh_pr=min(v_thresh[points_maxcory])
#    print(thresh_tas) ###TAS thresh
#    print(thresh_pr) ### PR thresh
#    print(maxcor_5)
#  }
#  return(list(res_mat_cor=res_mat_cor, thresh_tas=thresh_tas,
#              thresh_pr=thresh_pr, max_cor=maxcor_5))
#}
#
#### end essai
#
##### Conditional distrib
#compute_conditional_median<-function(pobs_data1,pobs_cond_data2,col_="black",
#                                      lab_X="F(Drought)",
#                                      lab_Y="F(HW)", pobs_to_calc=FALSE,
#                                      plot_pobs=FALSE){
#  if(pobs_to_calc==TRUE){
#    pobs_data1<-pobs(as.matrix(cbind(pobs_data1)))[,1]
#    pobs_cond_data2<-pobs(as.matrix(cbind(pobs_cond_data2)))[,1]
#  }
#  if(plot_pobs==TRUE){
#    plot(pobs_data1, pobs_cond_data2, xlim=c(0,1), ylim=c(0,1))
#  }
#
#
#  n_length=length(pobs_data1)
#  #pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#  #pobs_cond_data2<-pobs(as.matrix(cbind(cond_data2)))[,1]
#  res_ci_high_mean=res_ci_low_mean=c()
#  #plot(pobs_data1,pobs_cond_data2,
#  #             col="black",ylab="F(-PR)",xlab="F(TAS)")
#  
#  res_cond_median=c()
#  boot_ci95_cond_median_high=boot_ci90_cond_median_high=boot_ci99_cond_median_high=c()
#  boot_ci95_cond_median_low=boot_ci90_cond_median_low=boot_ci99_cond_median_low=c()
#  res_cond_CvM_v1=res_cond_CvM_v2=c()
#  u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#
#  #sample_generated_v1=seq(1/(10000+1),1-1/(10000+1),by=1/(10000+1))
#  sample_generated_v2=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  count_v=0
#  for(v in v_thresh){
#    count_v=count_v+1
#    id_dry=which(pobs_cond_data2>=v)
#    if(length(id_dry)<=5){
#      next
#    }else{
#      print(length(id_dry))
#      tmp_cond_median=median(pobs_data1[id_dry])
#      res_cond_median=c(res_cond_median,tmp_cond_median)#fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#      
#      boot_ci_median=sampling.median_CI(n_length, n_length-(count_v-1),B=3000)
#      boot_ci95_cond_median_low=c(boot_ci95_cond_median_low,boot_ci_median$res_ci95[1])
#      boot_ci95_cond_median_high=c(boot_ci95_cond_median_high,boot_ci_median$res_ci95[2])
#      boot_ci90_cond_median_low=c(boot_ci90_cond_median_low,boot_ci_median$res_ci90[1])
#      boot_ci90_cond_median_high=c(boot_ci90_cond_median_high,boot_ci_median$res_ci90[2])
#      boot_ci99_cond_median_low=c(boot_ci99_cond_median_low,boot_ci_median$res_ci99[1])
#      boot_ci99_cond_median_high=c(boot_ci99_cond_median_high,boot_ci_median$res_ci99[2])
#                   
#      ### subselection
##      tmp_cond_CvM_v1=cvm(pobs_data1[id_dry],sample_generated_v1)
#      tmp_cond_CvM_v2=cvm(pobs_data1[id_dry],sample_generated_v2)
##      res_cond_CvM_v1=c(res_cond_CvM_v1,tmp_cond_CvM_v1$pval)#fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#      res_cond_CvM_v2=c(res_cond_CvM_v2,tmp_cond_CvM_v2$pval)#fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#    }
#  }
#
#  color_vec_boot=rep("black",length(res_cond_median))
#  color_vec_boot[which(res_cond_median<boot_ci90_cond_median_low
#                       | res_cond_median>boot_ci90_cond_median_high)]="red"
# 
#  plot(res_cond_median,ylim=c(0,1),xaxt="n", col=color_vec_boot, ylab="",
#       xlab="")
#  title(ylab=paste0("med(",lab_Y,"|", lab_X, ")"),xlab=lab_X, cex.lab=1.2, line=2)
#  abline(h=0.5,lty=2,lwd=2)#lines(res_ci_high_median)
#  lines(boot_ci95_cond_median_low, col="gray66")
#  lines(boot_ci95_cond_median_high, col="gray66")
#  lines(boot_ci90_cond_median_low, col="gray90")
#  lines(boot_ci90_cond_median_high, col="gray90")
#  lines(boot_ci99_cond_median_low, col="gray0")
#  lines(boot_ci99_cond_median_high, col="gray0")
#  axis(1,seq.int(1,length(res_cond_median),10),round(v_thresh[seq.int(1,length(res_cond_median),10)],2))
#  #color_vec_pval_v1=color_vec_pval_v2=rep("black",length(res_cond_median))
#  #color_vec_pval_v1[which(res_cond_CvM_v1<=0.05)]="red"
# # color_vec_pval_v2[which(res_cond_CvM_v2<=0.10)]="red"
#  legend("topleft", inset=.02, 
#           legend= c("CI10","CI5","CI1"), col=c("gray90", "gray66",
#                                                      "gray0"),
#         pch=c(20,20,20) , horiz=FALSE, cex=1)
#
#
# #### CvM 
#  #plot(res_cond_CvM_v1,ylim=c(0,1),xaxt="n",ylab="pval",xlab=lab_X,col="blue")#color_vec_pval_v1)
##  plot(res_cond_CvM_v2,ylim=c(0,1),xaxt="n",col=color_vec_pval_v2, xlab="",
##       ylab="")
##  title(ylab="pval",xlab="F(Drought)",
##       cex.lab=1.2, line=2)
##  axis(1,seq.int(1,length(res_cond_median),10),round(v_thresh[seq.int(1,length(res_cond_median),10)],2))
##  abline(h=0.05, lty=2, lwd=2, col="gray66")
##  abline(h=0.01, lty=2, lwd=2, col="gray0")
##  abline(h=0.10, lty=2, lwd=2, col="gray90")
# 
#  return(res_cond_median)
#}
#
#
#
#compute_CvM_evolution_conditional_distrib<-function(pobs_data1_Ref,pobs_cond_data2_Ref,
#                                                    pobs_data1_Evol,
#                                                    pobs_cond_data2_Evol,col_="black",
#                                      lab_X="F(Drought)",
#                                      lab_Y="F(HW)", pobs_to_calc=FALSE){
#  if(pobs_to_calc==TRUE){
#    pobs_data1_Ref<-pobs(as.matrix(cbind(pobs_data1_Ref)))[,1]
#    pobs_cond_data2_Ref<-pobs(as.matrix(cbind(pobs_cond_data2_Ref)))[,1]
#    pobs_data1_Evol<-pobs(as.matrix(cbind(pobs_data1_Evol)))[,1]
#    pobs_cond_data2_Evol<-pobs(as.matrix(cbind(pobs_cond_data2_Evol)))[,1]
#  }
#
#  n_length=length(pobs_data1_Ref)
#  #res_ci_high_mean=res_ci_low_mean=c()
#  
#  #res_cond_median=c()
#  #boot_ci95_cond_median_high=boot_ci90_cond_median_high=boot_ci99_cond_median_high=c()
#  #boot_ci95_cond_median_low=boot_ci90_cond_median_low=boot_ci99_cond_median_low=c()
#  res_cond_CvM=c()
#  u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#
#  #sample_generated_v1=seq(1/(10000+1),1-1/(10000+1),by=1/(10000+1))
#  #sample_generated_v2=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#  count_v=0
#  for(v in v_thresh){
#    count_v=count_v+1
#    id_dry_Ref=which(pobs_cond_data2_Ref>=v)
#    id_dry_Evol=which(pobs_cond_data2_Evol>=v)
#    if(length(id_dry_Ref)<=5 | length(id_dry_Evol)<=5){
#      next
#    }else{
#      print(length(id_dry_Ref))
#      cond_distrib_Ref=pobs_data1_Ref[id_dry_Ref]
#      cond_distrib_Evol=pobs_data1_Evol[id_dry_Evol]
#      
#      ### subselection
#      tmp_cond_CvM=cvm(cond_distrib_Ref, cond_distrib_Evol)
#      res_cond_CvM=c(res_cond_CvM,tmp_cond_CvM$pval)#fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#    }
#  }
#
#  color_vec_boot=rep("black",length(res_cond_CvM))
#  color_vec_boot[which(res_cond_CvM< 0.10)]="red"
# 
#  plot(res_cond_CvM,ylim=c(0,1),xaxt="n", col=color_vec_boot, ylab="",
#       xlab="")
#  title(ylab=paste0("med(",lab_Y,"|", lab_X, ")"),xlab=lab_X, cex.lab=1.2, line=2)
#  abline(h=0.10,lty=2,lwd=2)#lines(res_ci_high_median)
#  axis(1,seq.int(1,length(res_cond_CvM),10),round(v_thresh[seq.int(1,length(res_cond_CvM),10)],2))
#
#  return(res_cond_CvM)
#}
#
#
#
#
#
#
#compute_CvM_1d<-function(ts_data, sliding_window_size=50, col_="black"){
#  ts_data_Ref=ts_data[1:sliding_window_size]
#  
#  nb_sliding_window=length(ts_data)
#  # nb_sliding_window=n_length-sliding_window_size+1
#  res_CvM=c()
#  for(k in 1:nb_sliding_window){
#    ts_data_Evol=ts_data[(1+(k-1)):(sliding_window_size+(k-1))]
#    
#    ### subselection
#    tmp_CvM=cvm(ts_data_Ref, ts_data_Evol)
#    res_CvM=c(res_CvM,tmp_CvM$pval)#fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#  }
#
#  plot(res_CvM,ylim=c(0,1),xaxt="n", col="black", ylab="CvM",
#       xlab="")
#  abline(h=0.10,lty=2,lwd=2)#lines(res_ci_high_median)
#  abline(h=0.05,lty=2,lwd=2)#lines(res_ci_high_median)
#  axis(1,seq.int(1,length(res_CvM),10),seq.int(1,length(res_CvM),10))
#  
#  return(res_CvM)
#}
#
#
#

#determine_ci_prob_bivar<-function(cop_object, ci_theta, prob_bivar=c(0.9,0.9)){
#  if(cop_object$family==3){name_family="Clayton"}
#  if(cop_object$family==4){name_family="Gumbel"}
#  if(cop_object$family==5){name_family="Frank"}
#  if(cop_object$family==6){name_family="Joe"}
#  d=2
#  res_prob=c("MLE"=NaN, rep(NaN,2))
#  
#  ####Proba for theta_mle
#  fitted_cop_mle<-BiCop(cop_object$family, ci_theta["MLE"])
#  res_prob["MLE"]=(1-BiCopCDF(1, prob_bivar[2], fitted_cop_mle)
#               -BiCopCDF(prob_bivar[1], 1, fitted_cop_mle)
#               + BiCopCDF(prob_bivar[1], prob_bivar[2], fitted_cop_mle))
#    
#  ####Proba for low_ci
#  test_fit_Cop_lowci=tryCatch({
#    BiCop(cop_object$family, ci_theta[2])
#        }, error=function(e){})
#  if(is.null(test_fit_Cop_lowci)){
#    print("Use of pCopula")
#    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    res_prob[2]=(1-pCopula(c(1, prob_bivar[2]), fitted_cop_lowci)
#                 -pCopula(c(prob_bivar[1], 1), fitted_cop_lowci)
#                 +pCopula(c(prob_bivar[1], prob_bivar[2]), fitted_cop_lowci))
#  }else{
#    fitted_cop_lowci<-test_fit_Cop_lowci
#    res_prob[2]=(1-BiCopCDF(1, prob_bivar[2], fitted_cop_lowci)
#                 -BiCopCDF(prob_bivar[1], 1, fitted_cop_lowci)
#                 + BiCopCDF(prob_bivar[1], prob_bivar[2], fitted_cop_lowci))
#  }
#
#  ####Proba for high_ci
#  test_fit_Cop_highci=tryCatch({
#    BiCop(cop_object$family, ci_theta[3])
#  }, error=function(e){})
#  if(is.null(test_fit_Cop_highci)){
#    print("Use of pCopula")
#    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
#    ## Evaluate this copula at the vector prob_bivar
#    res_prob[3]=(1-pCopula(c(1, prob_bivar[2]), fitted_cop_highci)
#                 -pCopula(c(prob_bivar[1], 1), fitted_cop_highci)
#                 +pCopula(c(prob_bivar[1], prob_bivar[2]), fitted_cop_highci))
#  }else{
#    print("Use of BiCop")
#    fitted_cop_highci<-test_fit_Cop_highci
#    res_prob[3]=(1-BiCopCDF(1, prob_bivar[2], fitted_cop_highci)
#                 -BiCopCDF(prob_bivar[1], 1, fitted_cop_highci)
#                 + BiCopCDF(prob_bivar[1], prob_bivar[2], fitted_cop_highci))
#  }
#  return(res_prob)
#}
#
## determine_ci_prob_bivar(cop, estim_ci_theta, prob_bivar=c(0.9,0.9))
#
# 
## 
# 
# 
# 
# 
# old2_compute_evolution_Condcorn<-function(data1, data2, col_="black", 
#                                   sens=1, 
#                                   name_plot="to_define",label_HW="HW",
#                                   boot_array_prop=NaN){ #data1 is tas, data2 is pr
#   n_length=length(data1)
#   #### Compute the pobs before dividing in period
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#  
# 
#   pobs_data1_Calib=pobs_data1[1:35]
#   pobs_data1_Proj=pobs_data1[36:70]
#   pobs_data2_Calib=pobs_data2[1:35]
#   pobs_data2_Proj=pobs_data2[36:70]
#   
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   res_mat_cor_Calib=res_mat_pval_Calib=res_mat_cor_Proj=res_mat_pval_Proj=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#   
#   #### Calib 
#   count_u=0
#   for(u in u_thresh){
#     count_u=count_u+1
#     count_v=0
#     for(v in v_thresh){
#       count_v=count_v+1
#       if(sens==1){id_Calib=which(pobs_data1_Calib>=u & pobs_data2_Calib>=v)}
#       
#       if(length(id_Calib)<=9){
#         next
#       }else{
#      ### subselection
#         pobs_data1_Calib_sub<-pobs_data1_Calib[id_Calib]
#         pobs_data2_Calib_sub<-pobs_data2_Calib[id_Calib]
#         res_mat_cor_Calib[count_u, count_v]=cor(pobs_data1_Calib_sub,pobs_data2_Calib_sub,method="spearman")
#         res_mat_pval_Calib[count_u,count_v]=cor.test(pobs_data1_Calib_sub,pobs_data2_Calib_sub,method="spearman")$p.value
#       }
#     }
#   }
#  
#   #### Proj 
#   count_u=0
#   for(u in u_thresh){
#     count_u=count_u+1
#     count_v=0
#     for(v in v_thresh){
#       count_v=count_v+1
#       if(sens==1){id_Proj=which(pobs_data1_Proj>=u & pobs_data2_Proj>=v)}
#       if(length(id_Proj)<=9){
#         next
#       }else{
#      ### subselection
#         pobs_data1_Proj_sub<-pobs_data1_Proj[id_Proj]
#         pobs_data2_Proj_sub<-pobs_data2_Proj[id_Proj]
#         res_mat_cor_Proj[count_u, count_v]=cor(pobs_data1_Proj_sub,pobs_data2_Proj_sub,method="spearman")
#         res_mat_pval_Proj[count_u,count_v]=cor.test(pobs_data1_Proj_sub,pobs_data2_Proj_sub,method="spearman")$p.value
#       }
#     }
#   }
#   
#   IND_pval_5=which(res_mat_pval_Calib<=0.05,arr.ind=TRUE)
#   if(length(IND_pval_5)==0){
#     thresh_tas=NaN
#     thresh_pr=NaN
#     maxcor_5=NaN
#   }else{
# #    maxcor_5=res_mat_cor[IND_pval_5][which.max(abs(res_mat_cor)[IND_pval_5])]
# #    points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
# #    points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#     plot(pobs_data1_Calib,pobs_data2_Calib, main=name_plot,ylim=c(0,1), xlim=c(0,1),
#        col=c(rep("blue",n_length/2),rep("red",n_length/2)),ylab="F(Drought)",xlab=paste0("F(",label_HW,")"))
#     image.plot(1:n_length, 1:n_length,res_mat_cor_Calib,
#              col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#     axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#     axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#     title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)")
#     legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=1)
# 
# 
#     ### Proj
#     plot(pobs_data1_Proj,pobs_data2_Proj, main=name_plot,ylim=c(0,1), xlim=c(0,1),
#        col=c(rep("blue",n_length/2),rep("red",n_length/2)),ylab="F(Drought)",xlab=paste0("F(",label_HW,")"))
#     image.plot(1:n_length, 1:n_length,res_mat_cor_Proj,
#              col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#     axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#     axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#     title(xlab=paste0("F(",label_HW,")"), ylab="F(Drought)")
#     legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=1)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############################################################################################################################
# ###############################################################################################################################
# ##### A la retraite    
# ###############################################################################################################################
# ###############################################################################################################################
#old_and_long_determine_ci_theta<-function(cop_object, data_U){
#  est_BiCop=BiCopEst(data_U[,1],data_U[,2],cop_object$family,method = "mle")
#  l_mle= est_BiCop$logLik
#  theta_mle=  est_BiCop$par
#  d=2
#  if(cop_object$family==3){
#    name_family="Clayton"
#    par_lim_inf=-0.9999
#  }
#  if(cop_object$family==4){
#    name_family="Gumbel"
#    par_lim_inf=1.0001
#  } 
#  if(cop_object$family==5){
#    name_family="Frank"
#    par_lim_inf=theta_mle-3
#  }
#  if(cop_object$family==6){
#    name_family="Joe"
#    par_lim_inf=1.0001
#  }
#
#  cop_object2 <- onacopulaL(name_family, list(theta_mle,1:d))
#
#
#  #### le choix de theta_mle est pas important dans onacopulaL ici (car on le
#  ### recalcule ensuite)
#  ## Estimation
#  #### efm peut soit crasher complet, soit a moitie..
#  test_efm=tryCatch({
#    efm <- emle(data_U, cop_object2)
#    pfm <- profile(efm)
#  }, error=function(e){})
#  if(is.null(test_efm)){
#    ci=c( "MLE"=theta_mle, c(NaN,NaN))
#  }else{
#    efm <- emle(data_U, cop_object2)
#    summary(efm) # using bblme's 'mle2' method
#    ## Profile likelihood plot [using S4 methods from bbmle/stats4] :
#    pfm <- profile(efm)
#    ci  <- confint(pfm, level=0.95)
#    ci<-c("MLE"=theta_mle, ci)
#  }
#
#  ### Parfois, ci comporte Na... on fait a la main:
#  if(sum(is.na(ci))>0){
#    print("ATTENTION, ci avec Na")
#    if(is.na(ci[2])){
#      print("ci_low a la main")
#      th4=seq(par_lim_inf,theta_mle,by=0.001)
#      Lt1 <- sapply(th4, function(th) LogL(th, cop_object2@copula, data_U))#, method="log.poly")) # default
#      th4[which.max(Lt1)]
#      # plot(Lt1~th4, main=paste0("theta_mle: ", round(theta_mle,2)))
#      # abline(v=theta_mle)
#      Chi=qchisq(.95, df=1)
#      LR=2*(l_mle-Lt1)
#      # plot(th4,LR,type='l')
#      # abline(h=Chi,col='red')
#      ci[2]=th4[which.min(abs(LR-Chi))]
#      print(paste0("CI comp: ", ci[2]))
#    }
#    if(is.na(ci[3])){
#      print("ci_high a la main")
#      th4=seq(theta_mle,theta_mle+1,by=0.001)
#      Lt1 <- sapply(th4, function(th) LogL(th, cop_object2@copula, data_U))#, method="log.poly")) # default
#      th4[which.max(Lt1)]
#      # plot(Lt1~th4, main=paste0("theta_mle: ", round(theta_mle,2)))
#      # abline(v=theta_mle)
#      Chi=qchisq(.95, df=1)
#      LR=2*(l_mle-Lt1)
#      # plot(th4,LR,type='l')
#      # abline(h=Chi,col='red')
#      ci[3]=th4[which.min(abs(LR-Chi))]
#      print(paste0("CI comp: ", ci[3]))
#    }
#  }
#  return(ci)
#}
#

 #
# ####TAS detrending:
# ### Option 1: dtrend1
# ### Pour chaque jour de l'été, on calcule les moyennes spatiales journalières
# ### Pour chaque année (ou été), on calcule la moyenne estivale des moyennes spatiale journalières ->41 years = 41 valeurs.
# ### on enleve la tendance de ces 41 valeurs sur les moyennes spatiales journalières de l'été.
# #### output: vector nb_day_by_season*nb_years de moyenne spatial detrendé
# dtrend1<-function(tas_array,Ind_season,nb_years,show_plot=FALSE){
#   tas_mat <- t(matrix(tas_array, nrow = nrow(tas_array) * ncol(tas_array)))
#   spavg_summer <- apply(tas_mat[Ind_season, ], 1, mean,na.rm=T) #Spatial average for the season
#   
#   ### Rearranging data
#   mean_spavg_by_summer=matrix(spavg_summer, byrow=FALSE, ncol = nb_years)
#   ### Compute spatial mean by season of T2
#   tendancy_spavg=apply(mean_spavg_by_summer, 2, mean)
#   ### Temporal trend 
#   lm_tendancy_spavg= lm(tendancy_spavg~c(1:nb_years))
#   
#   ### Detrend on daily spatial mean data
#   tmp_detrended_tas=mean_spavg_by_summer
#   for(i in 1:nb_years){
#     tmp_detrended_tas[,i]=tmp_detrended_tas[,i]-(lm_tendancy_spavg$coefficients[1]+lm_tendancy_spavg$coefficients[2]*i)
#   }
#   
#   ### Rearranging detrended daily spatial mean data
#   dt_spavg_summer=as.vector(tmp_detrended_tas)
#   ### Verification: Compute spatial mean by season of detrended T2
#   tendancy_dt_spavg=apply(tmp_detrended_tas, 2, mean)
#   ### Temporal trend 
#   lm_tendancy_dt_spavg= lm(tendancy_dt_spavg~c(1:nb_years))
#   
#   if(show_plot){
#     par(mfrow=c(2,2))
#     plot(spavg_summer,ylab="Summer spatial mean  (T2)",xlab="Summer days")
#     abline(lm(spavg_summer~c(1:length(spavg_summer))),col="red")
#     plot(tendancy_spavg,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy_spavg,col="red")
#     axis(1,seq.int(1,nb_years,9),seq.int(1951,2020,9))
#     plot(as.vector(dt_spavg_summer),ylab="Summer spatial mean  (Det. T2)",xlab="Summer days")
#     abline(lm(as.vector(dt_spavg_summer)~c(1:length(as.vector(dt_spavg_summer)))),col="red")
#     plot(tendancy_dt_spavg,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy_dt_spavg,col="red")
#     axis(1,seq.int(1,nb_years,9),seq.int(1951,2020,9))
#   }
#   
#   return(dt_spavg_summer)
# }
# 
# ### Same but for a matrix of points to detrend
# dtrend_by_points<-function(dt_matrix,nb_years, Ind_info=NaN, show_plot=FALSE){#dt_matrix matrix of nb_days_summer x point
#   res_detrended=matrix(NaN, nrow=dim(dt_matrix)[1], ncol=dim(dt_matrix)[2])
#   point_max= 1:dim(dt_matrix)[2]
#   for(point in point_max){
#     #i=Ind_info[point,1]
#     #j=Ind_info[point,2]
#     tmp_dt=matrix(dt_matrix[,point], byrow=FALSE, ncol=nb_years)
#     ### Compute mean by season/year of T2
#     tendancy=apply(tmp_dt, 2, mean)
#     ### Temporal trend 
#     lm_tendancy= lm(tendancy~c(1:nb_years))
# 
#     ### Detrend on daily spatial mean data
#     tmp_detrended_tas=tmp_dt
#     for(i in 1:nb_years){
#       tmp_detrended_tas[,i]=tmp_detrended_tas[,i]-(lm_tendancy$coefficients[1]+lm_tendancy$coefficients[2]*i)
#     }
#     
#     ### Rearranging detrended daily spatial mean data
#     res_detrended[,point]=as.vector(tmp_detrended_tas)
#     ### Verification: Compute spatial mean by season of detrended T2
#     tendancy_detrended=apply(matrix(res_detrended[,point], byrow=FALSE,
#                                     ncol=nb_years), 2, mean)
#     ### Temporal trend 
#     lm_tendancy_detrended= lm(tendancy_detrended~c(1:nb_years))
#      
#   }
# 
#   if(show_plot){
#     par(mfrow=c(2,2))
#     plot(dt_matrix[,1],ylab="Summer spatial mean  (T2)",xlab="Summer days")
#     abline(lm(dt_matrix[,1]~c(1:length(dt_matrix[,1]))),col="red")
#     plot(tendancy,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy,col="red")
#     axis(1,seq.int(1,nb_years,5),seq.int(1951,2020,5))
#     plot(as.vector(res_detrended[,1]),ylab="Summer spatial mean  (Det. T2)",xlab="Summer days")
#     abline(lm(as.vector(res_detrended[,1])~c(1:length(as.vector(res_detrended[,1])))),col="red")
#     plot(tendancy_detrended,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy_detrended,col="red")
#     axis(1,seq.int(1,nb_years,5),seq.int(1951,2020,5))
#   }
#   
#   return(res_detrended)
# }
# 
# dtrend_by_patch<-function(dt_matrix,nb_years, Ind_info=NaN, show_plot=FALSE){#dt_matrix matrix of nb_days_summer x point
#   res_detrended=matrix(NaN, nrow=dim(dt_matrix)[1], ncol=dim(dt_matrix)[2])
#  
#   spavg_summer <- apply(dt_matrix, 1, mean,na.rm=T) #Spatial average by nb_days_summer for the season
#   
#   ### Rearranging data by years in col
#   mean_spavg_by_summer=matrix(spavg_summer, byrow=FALSE, ncol = nb_years)
#   ### Compute spatial mean by season of T2
#   tendancy_spavg=apply(mean_spavg_by_summer, 2, mean)
#   ### Temporal trend 
#   lm_tendancy_spavg= lm(tendancy_spavg~c(1:nb_years))
#   point_max= 1:dim(dt_matrix)[2]
#   for(point in point_max){
#     tmp_dt=matrix(dt_matrix[,point], byrow=FALSE, ncol=nb_years)
#     ### Compute mean by season/year of T2
#     #tendancy=apply(tmp_dt, 2, mean)
#     ### Temporal trend 
#     #lm_tendancy= lm(tendancy~c(1:nb_years))
# 
#     ### Detrend on daily spatial mean data
#     tmp_detrended_tas=tmp_dt
#     for(i in 1:nb_years){
#       tmp_detrended_tas[,i]=tmp_detrended_tas[,i]-(lm_tendancy_spavg$coefficients[1]+lm_tendancy_spavg$coefficients[2]*i)
#     }
#     
#     ### Rearranging detrended daily spatial mean data
#     res_detrended[,point]=as.vector(tmp_detrended_tas)
#     ### Verification: Compute spatial mean by season of detrended T2
#     tendancy_detrended=apply(matrix(res_detrended[,point], byrow=FALSE,
#                                     ncol=nb_years), 2, mean)
#     ### Temporal trend 
#     lm_tendancy_detrended= lm(tendancy_detrended~c(1:nb_years))
#      
#   }
# 
#   if(show_plot){
#     par(mfrow=c(2,2))
#     plot(dt_matrix[,1],ylab="Summer spatial mean  (T2)",xlab="Summer days")
#     abline(lm(dt_matrix[,1]~c(1:length(dt_matrix[,1]))),col="red")
#     plot(tendancy_spavg,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy_spavg,col="red")
#     axis(1,seq.int(1,nb_years,5),seq.int(1951,2020,5))
#     plot(as.vector(res_detrended[,1]),ylab="Summer spatial mean  (Det. T2)",xlab="Summer days")
#     abline(lm(as.vector(res_detrended[,1])~c(1:length(as.vector(res_detrended[,1])))),col="red")
#     plot(tendancy_detrended,ylab="Summer spatial mean (T2)", xlab="Years",xaxt="n")
#     abline(lm_tendancy_detrended,col="red")
#     axis(1,seq.int(1,nb_years,5),seq.int(1951,2020,5))
#   }
#   
#   return(res_detrended)
# }
# 
#   }
#   return(list(res_mat_cor_Calib=res_mat_cor_Calib))
# }
# 
# 
# old_compute_nb_hot_days<-function(dt_data, nb_years){ #(take into account the time gap btw years
#   
#   ### Rearranging data par ete en colonne 
#   mean_spavg_by_summer=matrix(dt_data, byrow=FALSE, ncol = nb_years)
# 
#   quant90=quantile(mean_spavg_by_summer,probs=0.9)
#   res_nb_hot_days=apply(mean_spavg_by_summer, 2, function(x)sum(x>quant90))
#  
#   return(res_nb_hot_days) 
# }
# 
# 
# old_compute_persistence_cumule<-function(dt_data, nb_years, thresh_days=4,
#                                      probs_quant=0.9){
#   ### Rearranging data par ete en colonne 
#   mean_spavg_by_summer=matrix(dt_data, byrow=FALSE, ncol = nb_years)
#   
#   quant95=quantile(mean_spavg_by_summer,probs=probs_quant)
#   boolean_mat= (mean_spavg_by_summer>quant95)
#     
#   count_sequence<-function(bool,thresh_days=4){
#     count_seq=rle(bool)
#     tmp=count_seq$lengths[which(count_seq$values==TRUE)]
#     return(sum(tmp[which(tmp>=thresh_days)]))
#   }
#   res_persistence_cumule=apply(boolean_mat, 2,function(x)count_sequence(x,thresh_days=thresh_days))
#   return(res_persistence_cumule)
# }
# 
# 
# ### Quel seuil u,v maximise la correlation de Spearman entre HW et scaled_pr?
# old_compute_evolution_Condcorr<-function(data1, data2,  col_="black", name_plot="to_define"){ #data1 is tas, data2 is pr
#   data1_Calib=data1[1:35]
#   data1_Proj=data1[36:70]
#   data2_Calib=data2[1:35]
#   data2_Proj=data2[36:70]
#   
#   point=1
#   compute_matrix_Condcorr(data1_Calib,
#                         data2_Calib,  col_=col_RdBu, x3_version=FALSE,
#                                   sens=1, nloop=point,
#                                   name_plot=name_plot,
#                                   boot_array_prop= NaN) #data1 is tas, data2 is pr
# 
#   compute_matrix_Condcorr(data1_Proj,
#                         data2_Proj,  col_=col_RdBu, x3_version=FALSE,
#                                   sens=1, nloop=point,
#                                   name_plot=name_plot,
#                                   boot_array_prop= NaN) #data1 is tas, data2 is pr
# 
# 
#   pobs_data1_Calib<-pobs(as.matrix(cbind(data1_Calib)))[,1]
#   pobs_data1_Proj<-pobs(as.matrix(cbind(data1_Proj)))[,1]
#   pobs_data2_Calib<-pobs(as.matrix(cbind(data2_Calib)))[,1]
#   pobs_data2_Proj<-pobs(as.matrix(cbind(data2_Proj)))[,1]
# 
# 
#   plot(pobs_data1_Calib,pobs_data2_Calib, main=name_plot,
#          col="blue",ylab="F(Dry)",xlab="F(HW)")
#   plot(pobs_data1_Proj,pobs_data2_Proj, main=name_plot,
#          col="blue",ylab="F(Dry)",xlab="F(HW)")
# } 
# 
# 
# 
# 
# 
# 
# ### Quel seuil u,v maximise la correlation de Spearman entre HW et scaled_pr?
# Analysis_Copula_thresh_data<-function(data1,data2,col_="black",
#                                       x3_version=FALSE,sens=1){
#   n_length=length(data1)
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   res_mat_cor=res_mat_count=res_mat_pval=res_mat_selectCopula=res_mat_RP=matrix(rep(NaN,n_length*n_length),ncol=n_length)
# 
#   count_u=0
#   for(u in u_thresh){
#     count_u=count_u+1
#     count_v=0
#     for(v in v_thresh){
#       count_v=count_v+1
#       if(sens==1){id=which(pobs_data1>=u & pobs_data2>=v)}
#       if(sens==2){id=which(pobs_data1>=u & pobs_data2<=v)}
#       if(sens==3){id=which(pobs_data1<=u & pobs_data2<=v)}
#       if(sens==4){id=which(pobs_data1<=u & pobs_data2>=v)}
#       
#       #print(paste0("count_u: ",count_u, ", count_v: ", count_v, 
#       #             ", length id: ", length(id)))
#       #print(id)
#       if(length(id)<=9 & x3_version==FALSE){
#         next
#       }
#       if(length(id)<=9 & x3_version==TRUE){
#         next
#       }else{
#      ### subselection
#         pobs_data1_sub<-pobs(as.matrix(cbind(data1[id])))[,1]
#         pobs_data2_sub<-pobs(as.matrix(cbind(data2[id])))[,1]
#         res_mat_cor[count_u, count_v]=cor(pobs_data1_sub,pobs_data2_sub,method="spearman")
#         res_mat_count[count_u, count_v]=length(id)
#         res_mat_pval[count_u, count_v]=cor.test(pobs_data1_sub,pobs_data2_sub,method="spearman")$p.value
#         res_cop=selectedCopula(pobs_data1_sub, pobs_data2_sub)
#         #print(res_cop)
#         res_mat_selectCopula[count_u,count_v]=(res_cop$family!=0)
#         #### Compute RP
#         if(res_cop$family==0){
# #          mu=mean(diff(id))
# #          F_tas_seuil=ecdf(data1[id])
# #          F_pr_seuil=ecdf(data2[id])
# #          u_tas90=F_tas_seuil(quantile(data1,0.9))
# #          print(u_tas90)
# #          u_pr90=F_pr_seuil(quantile(data2,0.9))
# #          print(u_pr90)
#           res_mat_RP[count_u, count_v]=100#mu/(1-u_tas90-u_pr90+u_tas90*u_pr90)
# #          print(mu)
# #          print(res_mat_RP[count_u, count_v])
#         }else{
#           print(paste0("family_nuber",res_cop$family))
#           if(res_cop$family>=20){ #si c'est rotated, c'est bizarre...
#             next
#           }else{
# 
#             theta_cop=res_cop$par
#             if(res_cop$family %in% c(6,16)){
#               family_name="Joe"} ### Joe ou survival
#             if(res_cop$family %in% c(5,15)){
#               family_name="F"} ### Frank
#             if(res_cop$family %in% c(4,14)){
#               family_name="G"} ### Gumbel
#             if(res_cop$family %in% c(3,13)){
#               family_name="C"} ### Clayton
#  
#             if(res_cop$family<=9){
#               U=matrix(c(pobs_data1_sub,pobs_data2_sub),byrow=FALSE,ncol=2)
#             }
#             ## Si survival, on flip x and y
#             if(res_cop$family>=10 & res_cop$family<=19){ #survival 180
#               U=matrix(c(1-pobs_data1_sub,1-pobs_data2_sub),byrow=FALSE,ncol=2)
#             }
#             if(res_cop$family>=20 & res_cop$family<=29){ #90 degrees, on flip x
#               U=matrix(c(1-pobs_data1_sub,pobs_data2_sub),byrow=FALSE,ncol=2)
#             }
#             if(res_cop$family>=30 & res_cop$family<=39){ #270 degrees, on flip y
#               U=matrix(c(pobs_data1_sub,1-pobs_data2_sub),byrow=FALSE,ncol=2)
#             }
#             cop <- onacopulaL(family_name, list(theta_cop,1:2))
#             
#         #    #### Evaluate RP for tas>0.9 and pr>0.9
#             print(id)
#             mu=mean(diff(id))
#             if(x3_version){
#               x3_time_elapsing=rep(1:41, each=3)+rep(c(0,1/12,2/12),41)
#               print(x3_time_elapsing[id])
#               mu=mean(diff(x3_time_elapsing[id]))
#             }
#             print(paste0("Average time elapsing btw selected pairs: ", mu))
#         #    
#             F_tas_seuil=ecdf(data1[id])
#             F_pr_seuil=ecdf(data2[id])
#             u_tas90=F_tas_seuil(quantile(data1,0.9))
#             # print(u_tas90)
#             u_pr90=F_pr_seuil(quantile(data2,0.9))
#             # print(u_pr90)
#         
#             res_mat_RP[count_u, count_v]=mu/(1-u_tas90-u_pr90+BiCopCDF(u_tas90, u_pr90, res_cop$family,
#                                                  theta_cop))
#             print(paste0("essaiRP",res_mat_RP[count_u, count_v]))
#         #  abline(v=theta,col='red')
#         #  abline(h=res_RP,col='red')
#         #  print(res_RP)
#           }
#         }
#       }
#     }
#   }
# 
#   plot(pobs_data1,pobs_data2,
#        col="black",ylab="F(-PR)",xlab="F(TAS)")
# #  image.plot(1:n_length, 1:n_length,res_mat_cor,
# #             col=col_, xaxt="n", yaxt="n",zlim=c(-1,1))
#   image.plot(1:n_length, 1:n_length,res_mat_cor,
#              col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#   legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.8)
#  
#   points_pvalx_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,1]
#   points_pvaly_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,2]
#   points(points_pvalx_10, points_pvaly_10, pch=".", cex=3, col="gray90")
#   points_pvalx_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,1]
#   points_pvaly_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,2]
#   points(points_pvalx_5, points_pvaly_5, pch=".", cex=3, col="gray66")
#   points_pvalx_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,1]
#   points_pvaly_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,2]
#   points(points_pvalx_1, points_pvaly_1, pch=".", cex=3, col="gray0")
#   points_minpvalx=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,1]
#   points_minpvaly=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,2]
#   points(points_minpvalx, points_minpvaly, pch=0, cex=1, col="green")
#   IND_pval_5=which(res_mat_pval<=0.05,arr.ind=TRUE)
#   maxcor_5=max(res_mat_cor[IND_pval_5],na.rm=TRUE)
#   points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
#   points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#   points(points_maxcorx, points_maxcory, pch=0, cex=1, col="red")
#   image.plot(1:n_length, 1:n_length,res_mat_selectCopula,
#              col=col_, xaxt="n", yaxt="n", xlab="", ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#   image.plot(1:n_length, 1:n_length,res_mat_RP,
#              col=col_[33:64], xaxt="n", yaxt="n", zlim=c(0,100), xlab="",
#              ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#  
#   return(res_mat_cor)
# }
# 
# 
# 
# 
# ### Quel seuil u,v maximise la correlation de Spearman entre HW et scaled_pr?
# Analysis_4directions_Copula<-function(data1,data2,col_="black",
#                                       x3_version=FALSE){
#   n_length=length(data1)
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   
#   plot(pobs_data1,pobs_data2,
#        col="black",ylab="F(-PR)",xlab="F(TAS)")
#   plot.new()
#   for(sens in 1:4){
#     res_mat_cor=res_mat_count=res_mat_pval=res_mat_selectCopula=res_mat_RP=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#     count_u=0
#     for(u in u_thresh){
#       count_u=count_u+1
#       count_v=0
#       for(v in v_thresh){
#         count_v=count_v+1
#         if(sens==1){id=which(pobs_data1>=u & pobs_data2>=v)}
#         if(sens==2){id=which(pobs_data1>=u & pobs_data2<=v)}
#         if(sens==3){id=which(pobs_data1<=u & pobs_data2<=v)}
#         if(sens==4){id=which(pobs_data1<=u & pobs_data2>=v)}
#         
#         if(length(id)<=9 & x3_version==FALSE){
#           next
#         }
#         if(length(id)<=9 & x3_version==TRUE){
#           next
#         }else{
#           ### subselection
#           pobs_data1_sub<-pobs(as.matrix(cbind(data1[id])))[,1]
#           pobs_data2_sub<-pobs(as.matrix(cbind(data2[id])))[,1]
#           res_mat_cor[count_u, count_v]=cor(pobs_data1_sub,pobs_data2_sub,method="spearman")
#           res_mat_count[count_u, count_v]=length(id)
#           res_mat_pval[count_u, count_v]=cor.test(pobs_data1_sub,pobs_data2_sub,method="spearman")$p.value
#         }
#       }
#     }
#     image.plot(1:n_length, 1:n_length,res_mat_cor,
#                col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#     axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#     axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#     title(xlab="F(TAS)", ylab="F(-PR)")
#     legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                                            "gray0", "green", "red"),
#            pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.8)
#     
#     points_pvalx_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,1]
#     points_pvaly_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,2]
#     points(points_pvalx_10, points_pvaly_10, pch=".", cex=3, col="gray90")
#     points_pvalx_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,1]
#     points_pvaly_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,2]
#     points(points_pvalx_5, points_pvaly_5, pch=".", cex=3, col="gray66")
#     points_pvalx_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,1]
#     points_pvaly_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,2]
#     points(points_pvalx_1, points_pvaly_1, pch=".", cex=3, col="gray0")
#     points_minpvalx=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,1]
#     points_minpvaly=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,2]
#     points(points_minpvalx, points_minpvaly, pch=0, cex=1, col="green")
#     IND_pval_5=which(res_mat_pval<=0.05,arr.ind=TRUE)
#     maxcor_5=max(res_mat_cor[IND_pval_5],na.rm=TRUE)
#     points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
#     points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#     points(points_maxcorx, points_maxcory, pch=0, cex=1, col="red")
#   }
#   return(res_mat_cor)
# }
# 
# 
# 
# 
# 
# Bootstrap_Copula<-function(data1,data2,B=100,col_="black",
#                                       x3_version=FALSE){
#   n_length=length(data1)
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   res_mat_cor=res_mat_count=res_mat_pval=res_Bootmat_pval=matrix(rep(NaN,n_length*n_length),ncol=n_length)
#   res_Bootmat_cor=array(rep(NaN,n_length*n_length*B),dim=c(n_length, n_length, B))
# 
#   for(b in 1:B){
#     Bdata1=sample(data1, n_length, replace=F)
#     Bdata2=sample(data2, n_length, replace=F)
#     pobs_Bdata1<-pobs(as.matrix(cbind(Bdata1)))[,1]
#     pobs_Bdata2<-pobs(as.matrix(cbind(Bdata2)))[,1]
#     count_u=0
#     for(u in u_thresh){
#       count_u=count_u+1
#       count_v=0
#       for(v in v_thresh){
#         count_v=count_v+1
#         id=which(pobs_Bdata1>=u & pobs_Bdata2>=v)
#         print(paste0("count_u: ",count_u, ", count_v: ", count_v, 
#                      ", length id: ", length(id)))
#         print(id)
#         if(length(id)<=9 & x3_version==FALSE){
#           next
#         }
#         if(length(id)<=9 & x3_version==TRUE){
#           next
#         }else{
#        ### subselection
#           pobs_Bdata1_sub<-pobs(as.matrix(cbind(Bdata1[id])))[,1]
#           pobs_Bdata2_sub<-pobs(as.matrix(cbind(Bdata2[id])))[,1]
#           res_Bootmat_cor[count_u, count_v,b]=cor(pobs_Bdata1_sub,pobs_Bdata2_sub,method="spearman")
#           print(paste0("Boot ", b, ", cor computed",cor(pobs_Bdata1_sub,pobs_Bdata2_sub,method="spearman")))
#           }
#         }
#       }
#   }
# 
#   count_u=0
#   for(u in u_thresh){
#     count_u=count_u+1
#     count_v=0
#     for(v in v_thresh){
#       count_v=count_v+1
#       id=which(pobs_data1>=u & pobs_data2>=v)
#   #    print(paste0("count_u: ",count_u, ", count_v: ", count_v, 
#   #                 ", length id: ", length(id)))
#   #    print(id)
#       if(length(id)<=9 & x3_version==FALSE){
#         next
#       }
#       if(length(id)<=9 & x3_version==TRUE){
#         next
#       }else{
#      ### subselection
#         pobs_data1_sub<-pobs(as.matrix(cbind(data1[id])))[,1]
#         pobs_data2_sub<-pobs(as.matrix(cbind(data2[id])))[,1]
#         res_mat_cor[count_u, count_v]=cor(pobs_data1_sub,pobs_data2_sub,method="spearman")
#         res_mat_count[count_u, count_v]=length(id)
#         res_mat_pval[count_u, count_v]=cor.test(pobs_data1_sub,pobs_data2_sub,method="spearman")$p.value
#         res_Bootmat_pval[count_u, count_v] = ecdf(res_Bootmat_cor[count_u,
#                                                   count_v,])(res_mat_cor[count_u,count_v])
# # ou se place la vraie valeur de correlation dans la distribution du bootstrap?  
#         }
#       }
#     }
# 
# 
#   plot(pobs_data1,pobs_data2,
#        col="black",ylab="F(-PR)",xlab="F(TAS)")
#   image.plot(1:n_length, 1:n_length,res_mat_cor,
#              col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#   legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.8)
#  
#   points_pvalx_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,1]
#   points_pvaly_10=which(res_mat_pval<=0.10,arr.ind=TRUE)[,2]
#   points(points_pvalx_10, points_pvaly_10, pch=".", cex=3, col="gray90")
#   points_pvalx_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,1]
#   points_pvaly_5=which(res_mat_pval<=0.05,arr.ind=TRUE)[,2]
#   points(points_pvalx_5, points_pvaly_5, pch=".", cex=3, col="gray66")
#   points_pvalx_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,1]
#   points_pvaly_1=which(res_mat_pval<=0.01,arr.ind=TRUE)[,2]
#   points(points_pvalx_1, points_pvaly_1, pch=".", cex=3, col="gray0")
#   points_minpvalx=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,1]
#   points_minpvaly=which(res_mat_pval==min(res_mat_pval,na.rm=T),arr.ind=TRUE)[,2]
#   points(points_minpvalx, points_minpvaly, pch=0, cex=1, col="green")
#   IND_pval_5=which(res_mat_pval<=0.05,arr.ind=TRUE)
#   maxcor_5=max(res_mat_cor[IND_pval_5],na.rm=TRUE)
#   points_maxcorx=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,1]
#   points_maxcory=which(res_mat_cor==maxcor_5,arr.ind=TRUE)[,2]
#   points(points_maxcorx, points_maxcory, pch=0, cex=1, col="red")
# 
#   image.plot(1:n_length, 1:n_length,apply(res_Bootmat_cor,c(1,2),mean,na.rm=T),
#              col=col_, xaxt="n", yaxt="n",zlim=c(-1,1), xlab="", ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#   legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.8)
#   image.plot(1:n_length, 1:n_length,res_Bootmat_pval,
#              col=col_[33:64], xaxt="n", yaxt="n",zlim=c(0,1), xlab="", ylab="")
#   axis(1, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))#,seq.int(0,1,length.out=n_length),
#   axis(2, seq.int(1,n_length,length.out=6), seq.int(0,1,length.out=6))
#   title(xlab="F(TAS)", ylab="F(-PR)")
#   legend("topright", inset=.02, 
#            legend= c("pval10","pval5","pval1", "minpval","maxcorr"), col=c("gray90", "gray66",
#                                                       "gray0", "green", "red"),
#          pch=c(20,20,20, 0, 0), horiz=FALSE, cex=0.8)
#  
#   print(paste0("ok"))
#   print(dim(res_Bootmat_cor))
#   print(res_Bootmat_cor[20,20,]) 
#   plot(density(res_Bootmat_cor[28,1,],na.rm=T))
#   abline(v=res_mat_cor[28,1], col="red")
#   return(res_mat_cor)
# }
# 
# 
# 
# 
# 
# 
# #from function_EEA.R
# #For application, see EEA_Project/EEA_Method_of_Moments.R
# fit_GEV_with_Method_of_Momentsb<-function(data_gev){ #according to Hosking et al. 1985
#   n=length(data_gev)
#   data_gev=sort(data_gev) #sorting for order statistics 
#   #Computing b0
#   b0=mean(data_gev)
#   
#   #Computing b1
#   sum_gini=0
#   for(i in 2:n){
#     for(j in 1:(i-1)){
#       sum_gini=sum_gini+abs(data_gev[i]-data_gev[j])
#     }
#   }
#   twob1mb0=(2/(n*(n-1)))*0.5*sum_gini
#   b1=(twob1mb0+b0)/2
#   
#   #Computing b2
#   triple_sum=0
#   for(i in 3:n){
#     # if(i%%77==0){print(i)}
#     for(j in 2:(i-1)){
#       for(k in 1:(j-1)){
#         triple_sum=triple_sum+(data_gev[i]-2*data_gev[j]+data_gev[k])
#       }
#     }
#   }
#   sixb2m6b1pb0=(1/3)*(6/(n*(n-1)*(n-2)))*triple_sum
#   b2=(sixb2m6b1pb0-b0+6*b1)/6
#   
#   #estimateur shape:
#   c=(2*b1-b0)/(3*b2-b0)-(log(2)/log(3))
#   estim_shape=7.8590*c+2.9554*c*c
#   #estimateur scale:
#   estim_scale=((2*b1-b0)*estim_shape)/(gamma(1+estim_shape)*(1-2^(-estim_shape)))
#   #estimateur loc:
#   estim_loc=b0+estim_scale*(gamma(1+estim_shape)-1)/estim_shape
#   
#   #On retourne shape_h*(-1) for convention
#   return(list(shape_h=(-1)*estim_shape,scale_h=estim_scale,loc_h=estim_loc))
# }
# 
# 
# 
# #### Conditional distrib
# Analysis_CondDistrib<-function(data1,data2,col_="black",
#                                       x3_version=FALSE){
#   n_length=length(data1)
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
# 
#   plot(pobs_data1,pobs_data2,
#        col="black",ylab="F(-PR)",xlab="F(TAS)")
#   
#   res_mat_pargev=matrix(rep(NaN,3*n_length),ncol=n_length)
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   
#   data1_dry=data1[which(pobs_data2>=0.5)]
#   data1_wet=data1[which(pobs_data2<0.5)]
#   # plot(density(data1_dry),col='orange')
#   # lines(density(data1_wet),col='blue')
# 
#   count_v=0
#   for(v in v_thresh){
#     count_v=count_v+1
#     id_dry=which(pobs_data2>=v)
#     if(length(id_dry)<=5 & x3_version==FALSE){
#         next
#     }
#     if(length(id_dry)<=5 & x3_version==TRUE){
#         next
#       }else{
#         ### subselection
#         data1_dry_sub=data1[id_dry]
#         fit=fit_GEV_with_Method_of_Momentsb(data1_dry_sub)
#         res_mat_pargev[1, count_v]=fit$loc_h
#         res_mat_pargev[2, count_v]=fit$scale_h
#         res_mat_pargev[3, count_v]=fit$shape_h
#       }
#   }
#   
#   image.plot(1:n_length,1, as.matrix(res_mat_pargev[1,]),main="loc",col=col_)
#   image.plot(1:n_length, 1,as.matrix(res_mat_pargev[2,]),main="scale",col=col_)
#   image.plot(1:n_length,1,as.matrix(res_mat_pargev[3,]),main="shape",col=col_)
#   
#   return(data1_dry)
# }
# 
# 
# 
# Analysis_CondDistrib2<-function(data1,data2,col_="black",
#                                x3_version=FALSE){
#   n_length=length(data1)
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#   
#   plot(pobs_data1,pobs_data2,
#        col="black",ylab="F(-PR)",xlab="F(TAS)")
#   
#   res_mat_pargev=matrix(rep(NaN,3*n_length),ncol=n_length)
#   u_thresh=v_thresh=seq(1/(n_length+1),1-1/(n_length+1),by=1/(n_length+1))
#   
#   for(thresh in c(0.25,0.5,0.75,0.9)){
#     data1_dry=data1[which(pobs_data2>=thresh)]
#     data1_wet=data1[which(pobs_data2<thresh)]
#     plot(density(data1_dry),col='orange', main=paste0("dry PR>= ",thresh))
#     lines(density(data1_wet),col='blue')
#   }
# 
#   return(data1_dry)
# }
# 
# old_selectedCopula<-function(data1,data2){
#   pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
#   pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
#   a=BiCopSelect(pobs_data1,pobs_data2,familyset=c(3:6),indeptest = TRUE,selectioncrit = "BIC")
#   return(a)
# }


