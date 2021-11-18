library(fields)
library(igraph)
library(ismev)
library(energy)
library(MASS)
library(evd)
library(VineCopula)
library(goftest)

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
  (fitted<- evd::fpot(data_x, thresh_GPD))
  print(AIC(fitted))
  ### Estimate F_hat(points_to_estimate)
  p_of_x=evd::pgpd(points_to_estimate,loc =thresh_GPD,
              scale=fitted$estimate[1], shape=fitted$estimate[2])
  return(p_of_x)
}

Fm1_GPD<-function(data_x, thresh_GPD, probs_to_estimate){
  ### Fit F_hat with data_x
  (fitted<- evd::fpot(data_x, thresh_GPD))
  print(AIC(fitted))
  ### Estimate Fm1_hat
  q_of_x=evd::qgpd(probs_to_estimate, loc =thresh_GPD,
                   scale=fitted$estimate[1], shape=fitted$estimate[2])
  return(q_of_x)
}
####  Fit gaussian and estimate its inverse
F_Gaussian<-function(data_x, points_to_estimate){
  fitted=fitdistr(data_x, "normal")
  p_of_x=pnorm(points_to_estimate, mean = fitted$estimate[1], sd = fitted$estimate[2])
  return(p_of_x)
}

Fm1_Gaussian<-function(data_x, probs_to_estimate){
  fitted=fitdistr(data_x, "normal")
  q_of_x=qnorm(probs_to_estimate, mean = fitted$estimate[1], sd = fitted$estimate[2])
  return(q_of_x)
}

#### Fit GEV and estimate its inverse
F_GEV_for_min<-function(data_x, points_to_estimate){
  #### Reverse data
  minus_data_x=(-1)*data_x
  minus_points_to_estimate=(-1)*points_to_estimate
  ### Fit F_hat with data_x
  (fitted <- fgev(minus_data_x))
  ### Estimate 1-F_hat(points_to_estimate) (1- is here to reverse again)
  p_of_x= 1- evd::pgev(minus_points_to_estimate,loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  return(p_of_x)
}

Fm1_GEV_for_min<-function(data_x, probs_to_estimate){
  #### Reverse data
  minus_data_x=(-1)*data_x
  ### Fit F_hat with data_x
  (fitted <- fgev(minus_data_x))
  ### Estimate Fm1_hat
  q_of_x=evd::qgev(1-probs_to_estimate, loc =fitted$estimate[1], scale=fitted$estimate[2], shape=fitted$estimate[3])
  minus_q_of_x=(-1)*q_of_x
  return(minus_q_of_x)
}
selectedCopula<-function(data1,data2, pobs_to_calc=FALSE){ #from Clayton, Gumbel, Frank, Joe, ### not BB1, BB6, BB7, and BB8. 
  if(pobs_to_calc==TRUE){
    pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
    pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
  }else{
    pobs_data1<-data1
    pobs_data2<-data2
  }

  a=BiCopSelect(pobs_data1,pobs_data2,familyset=c(3:6),indeptest
                = FALSE,selectioncrit = "BIC", rotations=FALSE)
  return(a)
}

testCopula<-function(data1, data2, family, par=0, par2=0, pobs_to_calc=FALSE){
  if(pobs_to_calc==TRUE){
    pobs_data1<-pobs(as.matrix(cbind(data1)))[,1]
    pobs_data2<-pobs(as.matrix(cbind(data2)))[,1]
  }else{
    pobs_data1<-data1
    pobs_data2<-data2
  }
  a=BiCopGofTest(pobs_data1,pobs_data2, family, method = "white",max.df = 30, B = 100, obj = NULL)
  return(a)
}

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

#### si Profile Likelihood bug, on fait à la main
#Hofert 2011:  Likelihood inference for Archimedean copulas
LogL <- function(theta, acop, u, n.MC=0, ...) { # -(log-likelihood)
  sum(acop@dacopula(u, theta, n.MC=n.MC, log=TRUE, ...))
}
determine_ci_theta<-function(cop_object, data_U, ci_level=0.95){
  est_BiCop=BiCopEst(data_U[,1],data_U[,2],cop_object$family,method = "mle")
  l_mle= est_BiCop$logLik
  theta_mle=  est_BiCop$par
  d=2
  if(cop_object$family==3){
    name_family="Clayton"
    par_lim_inf=-0.9999
  }
  if(cop_object$family==4){
    name_family="Gumbel"
    par_lim_inf=1.0001
  }
  if(cop_object$family==5){
    name_family="Frank"
    par_lim_inf=theta_mle-3
  }
  if(cop_object$family==6){
    name_family="Joe"
    par_lim_inf=1.0001
  }
  cop_object2 <- onacopulaL(name_family, list(theta_mle,1:d))

  #### le choix de theta_mle est pas important dans onacopulaL ici (car on le
  ### recalcule ensuite)
  ## Estimation
  #### efm peut soit crasher complet, soit a moitie..
  test_efm=tryCatch({
    efm <- emle(data_U, cop_object2)
    pfm <- profile(efm)
  }, error=function(e){})
  if(is.null(test_efm)){
    ci=c( "MLE"=theta_mle, c(NaN,NaN))
  }else{
    efm <- emle(data_U, cop_object2)
    summary(efm) # using bblme's 'mle2' method
    ## Profile likelihood plot [using S4 methods from bbmle/stats4] :
    pfm <- profile(efm)
    ci  <- confint(pfm, level=ci_level)
    ci<-c("MLE"=theta_mle, ci)
  }

  ### Parfois, ci comporte Na... on fait a la main:
  if(sum(is.na(ci))>0){
    print("ATTENTION, ci avec Na")
    if(is.na(ci[2])){
      print("ci_low a la main")
      th4=seq(theta_mle,par_lim_inf-0.0001,by=-0.0001) ### on va de theta mle a la borne inf
      Chi=qchisq(ci_level, df=1)
      Lt1=LogL(theta_mle, cop_object2@copula, data_U)
      LR=2*(l_mle-Lt1)
      k=1
      while(LR<=Chi & th4[k]>par_lim_inf){
        k=k+1
        # print(k)
        Lt1 <- LogL(th4[k], cop_object2@copula, data_U)
        LR=2*(l_mle-Lt1)
      }
      ci[2]=th4[k]
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
  return(ci)
}


determine_ci_prob_bivar_array<-function(cop_object, ci_theta, prob_vector1=NaN,
                                        prob_vector2=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
  if(cop_object$family==3){name_family="Clayton"}
  if(cop_object$family==4){name_family="Gumbel"}
  if(cop_object$family==5){name_family="Frank"}
  if(cop_object$family==6){name_family="Joe"}
  d=2

  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe

  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)

  ####Proba for low_ci
  test_fit_Cop_lowci=tryCatch({
    BiCop(cop_object$family, ci_theta[2])
  }, error=function(e){})
  if(is.null(test_fit_Cop_lowci)){
    #print("Use of pCopula")
    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_lowci)
                                                     -pCopula(c(x[1], 1), fitted_cop_lowci)
                                                     +pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 

    # for(p in 1:length(prob_vector1)){
    #   prob_bivar=c(prob_vector1[p], prob_vector2[p])
    #   print((1-pCopula(c(1, prob_bivar[2]), fitted_cop_lowci)
    #                  -pCopula(c(prob_bivar[1], 1), fitted_cop_lowci)
    #                  +pCopula(c(prob_bivar[1], prob_bivar[2]), fitted_cop_lowci)))#res_prob[p,1]=
    # }
  }else{
    fitted_cop_lowci<-test_fit_Cop_lowci
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_lowci)
                                                     -BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                     + BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
    # for(p in 1:length(prob_vector1)){
    #   prob_bivar=c(prob_vector1[p], prob_vector2[p])
    #   res_prob[p,1]=(1-BiCopCDF(1, prob_bivar[2], fitted_cop_lowci)
    #                  -BiCopCDF(prob_bivar[1], 1, fitted_cop_lowci)
    #                  + BiCopCDF(prob_bivar[1], prob_bivar[2], fitted_cop_lowci))
    # }
  }
  ####Proba for high_ci
  test_fit_Cop_highci=tryCatch({
    BiCop(cop_object$family, ci_theta[3])
  }, error=function(e){})
  if(is.null(test_fit_Cop_highci)){
    #print("Use of pCopula")
    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-pCopula(c(1, x[2]), fitted_cop_highci)
                                                     -pCopula(c(x[1], 1), fitted_cop_highci)
                                                     +pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
  }else{
    #print("Use of BiCop")
    fitted_cop_highci<-test_fit_Cop_highci
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_highci)
                                                     -BiCopCDF(x[1], 1, fitted_cop_highci)
                                                     + BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
  }

  fitted_cop_mle<-BiCop(cop_object$family, ci_theta["MLE"])
  ## Evaluate this copula at the vector prob_bivar
  tmp_mat=matrix(apply(expand_proba,1,function(x){(1-BiCopCDF(1, x[2], fitted_cop_mle)
                                                   -BiCopCDF(x[1], 1, fitted_cop_mle)
                                                   + BiCopCDF(x[1], x[2], fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
  res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
  return(res_prob)
}

TardiveFrost_determine_ci_prob_bivar_array<-function(cop_object, ci_theta, prob_vector1=NaN,
                                        prob_vector2=NaN){ #prob_matrix format:prob_vector2 x prob_vector1 x CI
  if(cop_object$family==3){name_family="Clayton"}
  if(cop_object$family==4){name_family="Gumbel"}
  if(cop_object$family==5){name_family="Frank"}
  if(cop_object$family==6){name_family="Joe"}
  d=2

  expand_proba=expand.grid(x=prob_vector1, y=prob_vector2) #x and y en colonne, x qui defile avec y fixe

  ### res_prob of format probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)
  res_prob=array(NaN, dim=c(length(prob_vector1), length(prob_vector2), 3))
  ### ATTENTION, quand plot avec image.plot c'est probs_to_eval_1 x probs_to_eval_2 x (MLE, 0.025, 0.975)

  ####Proba for low_ci
  test_fit_Cop_lowci=tryCatch({
    BiCop(cop_object$family, ci_theta[2])
  }, error=function(e){})
  if(is.null(test_fit_Cop_lowci)){
    #print("Use of pCopula")
    fitted_cop_lowci <- onacopulaL(name_family, list(ci_theta[2],1:d))
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){( pCopula(c(x[1], 1), fitted_cop_lowci)
                                                     -pCopula(c(x[1], x[2]), fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 

  }else{
    fitted_cop_lowci<-test_fit_Cop_lowci
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_lowci)
                                                     - BiCopCDF(x[1], x[2], fitted_cop_lowci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,2]<-tmp_mat#apply(tmp_mat, 2, rev) 
  }

  ####Proba for high_ci
  test_fit_Cop_highci=tryCatch({
    BiCop(cop_object$family, ci_theta[3])
  }, error=function(e){})
  if(is.null(test_fit_Cop_highci)){
    #print("Use of pCopula")
    fitted_cop_highci <- onacopulaL(name_family, list(ci_theta[3],1:d))
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(pCopula(c(x[1], 1), fitted_cop_highci)
                                                     -pCopula(c(x[1], x[2]), fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
  }else{
    #print("Use of BiCop")
    fitted_cop_highci<-test_fit_Cop_highci
    ## Evaluate this copula at the vector prob_bivar
    tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_highci)
                                                     - BiCopCDF(x[1], x[2], fitted_cop_highci))}),byrow=TRUE,ncol=length(prob_vector1))
    res_prob[,,3]<-tmp_mat#apply(tmp_mat, 2, rev) 
  }

  fitted_cop_mle<-BiCop(cop_object$family, ci_theta["MLE"])
  ## Evaluate this copula at the vector prob_bivar
  tmp_mat=matrix(apply(expand_proba,1,function(x){(BiCopCDF(x[1], 1, fitted_cop_mle)
                                                   - BiCopCDF(x[1], x[2], fitted_cop_mle))}),byrow=TRUE,ncol=length(prob_vector1))
  res_prob[,,1]<-tmp_mat#apply(tmp_mat, 2, rev) 
  return(res_prob)
}


postproc_image.plot<-function(mat){ # row and col in mat correspond to what is
#plotted by image.plot
   res_mat<-mat[nrow(mat):1,]
   res_mat=t(res_mat)
   res_mat<-res_mat[,ncol(res_mat):1]
 return(res_mat)
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

