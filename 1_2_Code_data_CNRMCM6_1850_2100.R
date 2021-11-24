### 2. Models
#### To launch from Ciclad
# # To get ERA5 grid at 0.25
rm(list=ls())
library(timeDate)
library(ncdf4)
library(fields)
#### Convert to RData
convert_netcdf_RData<-function(ncfile,var){
        ncname <- ncfile
        ncfname <- paste(ncname,".nc", sep="")
        ncin <- nc_open(ncfname)
        res=ncvar_get(ncin, var)
        nc_close(ncin)
        return(res)
}
##
system(paste0("cdo -seltimestep,1/10 /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/t2m_era5_5020.nc* /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/target_grid_0_25deg.nc"))
#
Mod="MIROC6" #IPSL CNRMCM6
##### To launch from Ciclad: /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/
Grid_ERA5_0_25=paste("/home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/target_grid_0_25deg.nc")

for(name_varphy in c("sfcWindmax","pr","tasmin", "tasmax", "tas")){#"tas", "pr", "sfcWindmax", "tasmin", "tasmax")){
        if(name_varphy %in%c("tas", "tasmin", "tasmax")){
                op_cdo="-addc,-273.15"
        }
        if(name_varphy %in%c("pr")){
                op_cdo="-mulc,86400"
        }
        if(name_varphy %in%c("sfcWindmax")){
                op_cdo=""
        }

        print(paste0(name_varphy, ": ", op_cdo))

        if(Mod=="CNRMCM6"){
                Mod_folder="CNRM-CERFACS"
                Mod_name="CNRM-CM6-1"
                Run_name="r1i1p1f2"
                G_name="gr"
        }
        if(Mod=="IPSL"){
                Mod_folder="IPSL"
                Mod_name="IPSL-CM6A-LR"
                Run_name="r1i1p1f1i"
                G_name="gr"
        }
        if(Mod=="MIROC6"){
                Mod_folder="MIROC"
                Mod_name="MIROC6"
                Run_name="r1i1p1f1"
                G_name="gn"
       }
       path_historical=paste0("/bdd/CMIP6/CMIP/", Mod_folder, "/", Mod_name, "/historical/", Run_name, "/day/", name_varphy, "/", G_name , "/latest/")
       path_ssp585=paste0("/bdd/CMIP6/ScenarioMIP/", Mod_folder, "/", Mod_name, "/ssp585/", Run_name, "/day/", name_varphy, "/", G_name, "/latest/")

       if(Mod=="IPSL"){
               system(paste0("cdo mergetime ",
                       "-remapbil,", Grid_ERA5_0_25 ,
                       " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                       path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_18500101-20141231.nc ",
                       "-remapbil,", Grid_ERA5_0_25 ,
                       " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                       path_ssp585, name_varphy, "_day_", Mod_name, "_ssp585_", Run_name, "_", G_name, "_20150101-21001231.nc ",
                       "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
       }

       if(Mod=="CNRMCM6"){
               if(name_varphy %in% c("tas", "tasmin", "tasmax")){
                       system(paste0("cdo mergetime ",
                               "-remapbil,", Grid_ERA5_0_25 ,
                               " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_18500101-20141231.nc ",
                               "-remapbil,", Grid_ERA5_0_25 ,
                               " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               path_ssp585, name_varphy, "_day_", Mod_name, "_ssp585_", Run_name, "_", G_name, "_20150101-21001231.nc ",
                               "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
               }
               if(name_varphy %in% c("pr", "sfcWindmax",)){
                       system(paste0("cdo mergetime ",
                               "-remapbil,", Grid_ERA5_0_25 ,
                               " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_18500101-19491231.nc ",
                               "-remapbil,", Grid_ERA5_0_25 ,
                               " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_19500101-20141231.nc ",
                               "-remapbil,", Grid_ERA5_0_25 ,
                               " -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               path_ssp585, name_varphy, "_day_", Mod_name, "_ssp585_", Run_name, "_", G_name, "_20150101-21001231.nc ",
                               "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
               }
       }

       if(Mod=="MIROC6"){
               files_MIROC_hist=list.files(path_historical)
               remapbil_cmd_hist=""
               for(f in files_MIROC_hist){
                       remapbil_cmd_hist=paste0(remapbil_cmd_hist, "-remapbil,", Grid_ERA5_0_25, " -selyear,1850/2100 ", op_cdo, " -del29feb ", path_historical, f, " ")
               }
               files_MIROC_ssp585=list.files(path_ssp585)
               remapbil_cmd_ssp585=""
               for(f in files_MIROC_ssp585){
                       remapbil_cmd_ssp585=paste0(remapbil_cmd_ssp585, "-remapbil,", Grid_ERA5_0_25, " -selyear,1850/2100 ", op_cdo, " -del29feb ", path_ssp585, f, " ")
               }
               if(name_varphy %in% c("tas", "tasmin", "tasmax", "pr", "sfcWindmax")){
                       system(paste0("cdo mergetime ",
                               remapbil_cmd_hist,
                               remapbil_cmd_ssp585,
                               "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
                               #"-remapbil,", Grid_ERA5_0_25 ,
                               #" -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               #path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_18500101-18591231.nc ",
                               #"-remapbil,", Grid_ERA5_0_25 ,
                               #" -selyear,1850/2100 ", op_cdo, " -del29feb ",
                               #path_ssp585, name_varphy, "_day_", Mod_name, "_ssp585_", Run_name, "_", G_name, "_20150101-21001231.nc ",
                               #"/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
               }
       }
       tmp_varphy_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/",name_varphy,"_day_", Mod, "_1850_2100_Europe"),name_varphy)

       LON_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/",name_varphy,"_day_", Mod, "_1850_2100_Europe"),"longitude")
       #### ATTENTION ROTATION LAT and so need to rotate  on row axis
       LAT_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe"),"latitude")

       print(LON_Europe)
       print(LAT_Europe)
       #### Rotation
       rotate <- function(x) t(apply(x, 1, rev))

       good_tmp_varphy_Europe=array(NaN,dim=c(dim(tmp_varphy_Europe)))

       #### si LAT in decreasing order
       if(LAT_Europe[1]>=LAT_Europe[2]){
               LAT_Europe=rev(LAT_Europe)
               for(i in 1:dim(tmp_varphy_Europe)[3]){
                       if(i%%7777==1){print(round(i/dim(tmp_varphy_Europe)[3]*100,2))}
                       good_tmp_varphy_Europe[,,i]<-rotate(tmp_varphy_Europe[,,i])
               }
       }else{
               good_tmp_varphy_Europe<- tmp_varphy_Europe
       }

       if(LON_Europe[1]>=LON_Europe[2]){
               print(pasbon)
       }

       setwd("/scratchu/bfrancois/")
      ##### Preprocess Mod

      rm(tmp_varphy_Europe)
      gc()

      #load commune information to extract a specific region in France
      ##obtained from base=read.csv("http://freakonometrics.free.fr/popfr19752010.csv",header=TRUE)
      #rm(list=ls())

      load("/home/bfrancois/Compound_Event/TCE5/Data/commune_France.RData")
      print(LON_Europe)
      print(LAT_Europe)
      print(dim(good_tmp_varphy_Europe))

      #### France
      library(maps)

      #### Define France coordinates
      lon_coord_france=which(LON_Europe>=min(commune_France$long)-0.5
                                   & LON_Europe<=max(commune_France$long) )
      lat_coord_france=which(LAT_Europe>=min(commune_France$lat)& LAT_Europe<=max(commune_France$lat)+0.5)

      LON_France=LON_Europe[lon_coord_france]
      LAT_France=LAT_Europe[lat_coord_france]

      good_tmp_varphy_France<-good_tmp_varphy_Europe[lon_coord_france, lat_coord_france,]


      old.par <- par(mar = c(0, 0, 0, 0))
      par(old.par)

      commune_France_mat=matrix(NaN,ncol=2, nrow=length(commune_France$long))
      commune_France_mat[,1]=commune_France$long
      commune_France_mat[,2]=commune_France$lat

      good_tmp_varphy_France_cut<- array(NaN, dim=c(dim(good_tmp_varphy_France)))

      #ess plus rapide
      COUNT=0
      for(i in 1:length(LON_France)){
        print(i)
        for(j in 1:length(LAT_France)){
          if(LON_France[i]>8){next}else{
            coord=c(LON_France[i],LAT_France[j])
            subset=which((commune_France_mat[,1]>LON_France[i]-0.2) & (commune_France_mat[,1]<LON_France[i]+0.2) &
               (commune_France_mat[,2]>LAT_France[j]-0.2) & (commune_France_mat[,2]<LAT_France[j]+0.2),arr.ind=TRUE)

              if(length(subset)<=1){
              next}
      else{
                diff_square=apply(commune_France_mat[subset,],1,function(x) (x-coord)^2)
              sum_diff_square=apply(diff_square,2,sum)
              sqrt_sum_diff_square=sqrt(sum_diff_square)
              #print(min(sqrt_sum_diff_square))
              if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
                #print("ok")
                COUNT=COUNT+1
#               print(COUNT)
                good_tmp_varphy_France_cut[i,j,]=good_tmp_varphy_France[i,j,]
              }
            }
          }
        }
        }

        #end essai

        setwd("/scratchu/bfrancois/")
        pdf(paste0("newok", name_varphy,"_", Mod, "_Francenew.pdf"))
        par(mfrow=c(2,2))
        image.plot(LON_France,LAT_France,good_tmp_varphy_France_cut[,,1])
        world(add=T)
        image.plot(LON_France,LAT_France,good_tmp_varphy_France[,,1])
        world(add=T)
        dev.off()


        IND_France=which(!(is.na(good_tmp_varphy_France_cut)),arr.ind=T)

        setwd("/scratchu/bfrancois/")

        assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_France"), good_tmp_varphy_France_cut)



        ##### Bretagne
        good_tmp_varphy_Bretagne_cut=array(NaN, dim=c(dim(good_tmp_varphy_France)))

        count=0
        for(i in 1:length(LON_France)){
        #  print(i)
          for(j in 1:length(LAT_France)){
            if(LON_France[i]>0 | LAT_France[j]<46 | LAT_France[j]>=49 | count>=80){next}else{
              coord=c(LON_France[i],LAT_France[j])
              diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
              sum_diff_square=apply(diff_square,2,sum)
              sqrt_sum_diff_square=sqrt(sum_diff_square)
         #     print(min(sqrt_sum_diff_square))
              if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
                if(LON_France[i]>8){next}else{
                  count=count+1
         #         print("ok")
                  good_tmp_varphy_Bretagne_cut[i,j,] = good_tmp_varphy_France[i,j,]
                }
              }
            }
          }
        }

        #### convert to small arrays
        coord_lat_bretagne=which(LAT_France >46 & LAT_France <=49)
        coord_lon_bretagne=which(LON_France>=(-5.25) & LON_France<=(-1.5))
        LON_Bretagne=LON_France[coord_lon_bretagne]
        LAT_Bretagne=LAT_France[coord_lat_bretagne]


        good_tmp_varphy_Bretagne_cut=good_tmp_varphy_Bretagne_cut[coord_lon_bretagne,coord_lat_bretagne,]

        pdf(paste0("ok",name_varphy,"_", Mod, "_Bretagnenew.pdf"))
        image.plot(LON_Bretagne,LAT_Bretagne,good_tmp_varphy_Bretagne_cut[,,1])
        dev.off()

        IND_Bretagne=which(!(is.na(good_tmp_varphy_Bretagne_cut)),arr.ind=T)

        setwd(paste0("/home/bfrancois/Compound_Event/TCE5/Data/RefMod_Data/", Mod, "/"))
        assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_Bretagne"), good_tmp_varphy_Bretagne_cut)
        eval(parse(text=(paste0("save(", name_varphy, "_day_", Mod, "_1850_2100_Bretagne, LON_Bretagne, LAT_Bretagne, IND_Bretagne, file='", name_varphy, "_day_", Mod, "_1850_2100_Bretagne.RData')"))))
        rm(good_tmp_varphy_Bretagne_cut)
        gc()
        print("ok_final")

        ##### Burgundy
        good_tmp_varphy_Burgundy_cut=array(NaN, dim=c(dim(good_tmp_varphy_France)))

        count=0
        for(i in 1:length(LON_France)){
          print(i)
          for(j in 1:length(LAT_France)){
            if(LON_France[i]>=(5) | LON_France[i]<=(-1) | LAT_France[j]<=46 | LAT_France[j]>=49){next}else{
              coord=c(LON_France[i],LAT_France[j])
              diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
              sum_diff_square=apply(diff_square,2,sum)
              sqrt_sum_diff_square=sqrt(sum_diff_square)
              print(min(sqrt_sum_diff_square))
              if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
                if(LON_France[i]>8){next}else{
                 count=count+1
                 print("ok")
                 good_tmp_varphy_Burgundy_cut[i,j,] = good_tmp_varphy_France[i,j,]
               }
             }
           }
         }
       }

       #### convert to small arrays
       coord_lat_Burgundy=which(LAT_France >=46 & LAT_France <=49)
       coord_lon_Burgundy=which(LON_France>=(-1) & LON_France<=(5))
       LON_Burgundy=LON_France[coord_lon_Burgundy]
       LAT_Burgundy=LAT_France[coord_lat_Burgundy]


       good_tmp_varphy_Burgundy_cut=good_tmp_varphy_Burgundy_cut[coord_lon_Burgundy,coord_lat_Burgundy,]

       pdf(paste0("ok",name_varphy, "_", Mod,"_Burgundynew.pdf"))
       image.plot(LON_Burgundy,LAT_Burgundy,good_tmp_varphy_Burgundy_cut[,,1])
       dev.off()

       IND_Burgundy=which(!(is.na(good_tmp_varphy_Burgundy_cut)),arr.ind=T)

       setwd(paste0("/home/bfrancois/Compound_Event/TCE5/Data/RefMod_Data/", Mod, "/"))
       assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_Burgundy"), good_tmp_varphy_Burgundy_cut)
       eval(parse(text=(paste0("save(", name_varphy, "_day_", Mod, "_1850_2100_Burgundy, LON_Burgundy, LAT_Burgundy, IND_Burgundy, file='", name_varphy, "_day_", Mod, "_1850_2100_Burgundy.RData')"))))
       rm(good_tmp_varphy_Europe, good_tmp_varphy_France, good_tmp_varphy_France_cut, good_tmp_varphy_Burgundy_cut)
       gc()
       print("ok_final")
}

print(ok)
