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
CMIP6_Mod=c("MRIESM2", "MPIESM1LR", "MIROC6")
#CanESM5 CMCCESM2 "CNRMCM6", "CNRMCM6HR", "ECEARTH3", "FGOALSG3","GFDLCM4", "INMCM48", "INMCM50","IPSL", 
path_grid=paste0("/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1-HR/historical/r1i1p1f2/day/tas/gr/latest/tas_day_CNRM-CM6-1-HR_historical_r1i1p1f2_gr_18500101-18991231.nc")
system(paste0("cdo -sellonlatbox,-10,30,30,70 -seltimestep,1/10 ", path_grid, " /home/bfrancois/Compound_Event/TCE6/Data/target_grid_0.5deg_CNRMCM6HR.nc"))
#system(paste0("cdo -seltimestep,1/10 /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/t2m_era5_5020.nc* /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/tas/target_grid_0_25deg.nc"))
#
list_Mod=c("IPSL")
##### To launch from Ciclad: /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/
Grid_CNRMCM6HR_0_5=paste("/home/bfrancois/Compound_Event/TCE6/Data/target_grid_0.5deg_CNRMCM6HR.nc")

for(Mod in list_Mod){
	print(paste0("###########################", Mod, "#############################"))
	for(name_varphy in c("sfcWindmax", "pr", "tasmin", "tas")){#"tas", "pr", "sfcWindmax", "tasmin", "tasmax")){
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
			Run_name="r1i1p1f1"
			G_name="gr"
		}
		if(Mod=="MIROC6"){
			Mod_folder="MIROC"
			Mod_name="MIROC6"
			Run_name="r1i1p1f1"
			G_name="gn"
		}
		if(Mod=="MIROCES2L"){
			Mod_folder="MIROC"
			Mod_name="MIROC-ES2L"
			Run_name="r1i1p1f2"
			G_name="gn"
		}
		if(Mod=="INMCM48"){
			Mod_folder="INM"
			Mod_name="INM-CM4-8"
			Run_name="r1i1p1f1"
			G_name="gr1"
		}
		if(Mod=="INMCM50"){
			Mod_folder="INM"
			Mod_name="INM-CM5-0"
			Run_name="r1i1p1f1"
			G_name="gr1"
		}
		if(Mod=="NorESM2"){
			Mod_folder="NCC"
			Mod_name="NorESM2-LM"
			Run_name="r1i1p1f1"
			G_name="gn"
		}
		if(Mod=="CanESM5"){
			Mod_folder="CCCma"
			Mod_name="CanESM5"
			Run_name="r10i1p1f1"
			G_name="gn"
		}
		if(Mod=="BCCCSM2MR"){
			Mod_folder="BCC"
			Mod_name="BCC-CSM2-MR"
			Run_name="r1i1p1f1"
			G_name="gn"
		}
		if(Mod=="GFDLCM4"){
			Mod_folder="NOAA-GFDL"
			Mod_name="GFDL-CM4"
			Run_name="r1i1p1f1"
			G_name="gr1"
		}
		if(Mod=="MPIESM1LR"){
			Mod_folder="MPI-M"
			Mod_name="MPI-ESM1-2-LR"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="MRIESM2"){
			Mod_folder="MRI"
			Mod_name="MRI-ESM2-0"
			Run_name="r1i2p1f1"
			G_name="gn"
		}	
		if(Mod=="ACCESSESM1"){
			Mod_folder="CSIRO"
			Mod_name="ACCESS-ESM1-5"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="CNRMCM6HR"){
			Mod_folder="CNRM-CERFACS"
			Mod_name="CNRM-CM6-1-HR"
			Run_name="r1i1p1f2"
			G_name="gr"
		}	
		if(Mod=="AWIESM1"){
			Mod_folder="AWI"
			Mod_name="AWI-ESM-1-1-LR"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="FGOALSF3"){
			Mod_folder="CAS"
			Mod_name="FGOALS-f3-L"
			Run_name="r1i1p1f1"
			G_name="gr"
		}	
		if(Mod=="FGOALSG3"){
			Mod_folder="CAS"
			Mod_name="FGOALS-g3"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="CMCCESM2"){
			Mod_folder="CMCC"
			Mod_name="CMCC-ESM2"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="CMCCCM2"){
			Mod_folder="CMCC"
			Mod_name="CMCC-CM2-SR5"
			Run_name="r1i1p1f1"
			G_name="gn"
		}	
		if(Mod=="ECEARTH3"){
			Mod_folder="EC-Earth-Consortium"
			Mod_name="EC-Earth3"
			Run_name="r1i1p1f1"
			G_name="gr"
		}	
		
		path_historical=paste0("/bdd/CMIP6/CMIP/", Mod_folder, "/", Mod_name, "/historical/", Run_name, "/day/", name_varphy, "/", G_name , "/latest/")
		path_ssp585=paste0("/bdd/CMIP6/ScenarioMIP/", Mod_folder, "/", Mod_name, "/ssp585/", Run_name, "/day/", name_varphy, "/", G_name, "/latest/")
	

	        if(Mod=="IPSL"){
	                system(paste0("cdo mergetime ",
	                        "-remapbil,", Grid_CNRMCM6HR_0_5 ,
	                        " -selyear,1850/2100 ", op_cdo, " -del29feb ",
	                        path_historical, name_varphy, "_day_", Mod_name, "_historical_", Run_name, "_", G_name, "_18500101-20141231.nc ",
	                        "-remapbil,", Grid_CNRMCM6HR_0_5 ,
	                        " -selyear,1850/2100 ", op_cdo, " -del29feb ",
	                        path_ssp585, name_varphy, "_day_", Mod_name, "_ssp585_", Run_name, "_", G_name, "_20150101-21001231.nc ",
	                        "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
	        }

	
		if(Mod %in% CMIP6_Mod){
			files_Mod_hist=list.files(path_historical)
			remapbil_cmd_hist=""
			for(f in files_Mod_hist){
				remapbil_cmd_hist=paste0(remapbil_cmd_hist, "-remapbil,", Grid_CNRMCM6HR_0_5, " -selyear,1850/2100 ", op_cdo, " -del29feb ", path_historical, f, " ")
			}		
			files_Mod_ssp585=list.files(path_ssp585)
			remapbil_cmd_ssp585=""
			for(f in files_Mod_ssp585){
				remapbil_cmd_ssp585=paste0(remapbil_cmd_ssp585, "-remapbil,", Grid_CNRMCM6HR_0_5, " -selyear,1850/2100 ", op_cdo, " -del29feb ", path_ssp585, f, " ")
			}		
			if(name_varphy %in% c("tas", "tasmin", "tasmax", "pr", "sfcWindmax")){
				system(paste0("cdo mergetime ",
		        		remapbil_cmd_hist,
					remapbil_cmd_ssp585, 
			       		"/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
			}
		}
		
		if(Mod %in% c("ECEARTH3","FGOALSG3")){
			#### need to do tmp files before final merging
			cmd_final_merging=""
			length_subset=30
			
			files_Mod_hist=list.files(path_historical)
			if(Mod=="FGOALSG3"){ #Case for FGOALSG3 where historical goes to 2016
				coord_2015=which(grepl("2015", files_Mod_hist, fixed = TRUE))
				coord_2016=which(grepl("2016", files_Mod_hist, fixed = TRUE))
				files_Mod_hist=files_Mod_hist[-c(coord_2015, coord_2016)]
			}
			nb_subset_hist=trunc(length(files_Mod_hist)/length_subset)
			list_subset_hist=list()
			for(sub in 1:nb_subset_hist){
				list_subset_hist[[sub]]=(1+(sub-1)*length_subset):(length_subset*sub)
			}
			rest=length(files_Mod_hist) %% length_subset
			if(rest!=0){
				list_subset_hist[[nb_subset_hist+1]]=(1+(nb_subset_hist)*length_subset):((nb_subset_hist)*length_subset+rest)
			}
			for(l in list_subset_hist){
				remapbil_cmd_hist=""
				for(f in files_Mod_hist[l]){
					remapbil_cmd_hist=paste0(remapbil_cmd_hist, "-remapbil,", Grid_CNRMCM6HR_0_5, " -selyear,1850/2014 ", op_cdo, " -del29feb ", path_historical, f, " ")
				}		
				system(paste0("cdo mergetime ",
		        		remapbil_cmd_hist,
			       		"/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_historical_",l[1], "_",l[length(l)],"_Europe.nc"))
				cmd_final_merging=paste0(cmd_final_merging, "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_historical_",l[1], "_",l[length(l)],"_Europe.nc ")
			}

			files_Mod_ssp585=list.files(path_ssp585)
			nb_subset_ssp585=trunc(length(files_Mod_ssp585)/length_subset)
			list_subset_ssp585=list()
			for(sub in 1:nb_subset_ssp585){
				list_subset_ssp585[[sub]]=(1+(sub-1)*length_subset):(length_subset*sub)
			}
			rest=length(files_Mod_ssp585) %% length_subset
			if(rest!=0){
				list_subset_ssp585[[nb_subset_ssp585+1]]=(1+(nb_subset_ssp585)*length_subset):((nb_subset_ssp585)*length_subset+rest)
			}
			for(l in list_subset_ssp585){
				remapbil_cmd_ssp585=""
				for(f in files_Mod_ssp585[l]){
					remapbil_cmd_ssp585=paste0(remapbil_cmd_ssp585, "-remapbil,", Grid_CNRMCM6HR_0_5, " -selyear,2015/2100 ", op_cdo, " -del29feb ", path_ssp585, f, " ")
				}		
				system(paste0("cdo mergetime ",
		        		remapbil_cmd_ssp585,
			       		"/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_ssp585_",l[1], "_",l[length(l)],"_Europe.nc"))
				cmd_final_merging=paste0(cmd_final_merging, "/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_ssp585_",l[1], "_",l[length(l)],"_Europe.nc ")
			}
			#### Final merging:
			system(paste0("cdo mergetime ", cmd_final_merging,
			       		"/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe.nc"))
	
		}
		
		tmp_varphy_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/",name_varphy,"_day_", Mod, "_1850_2100_Europe"),name_varphy)
		
		LON_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/",name_varphy,"_day_", Mod, "_1850_2100_Europe"),"lon")
		#### ATTENTION ROTATION LAT and so need to rotate  on row axis
		LAT_Europe=convert_netcdf_RData(paste0("/scratchu/bfrancois/", name_varphy, "_day_", Mod, "_1850_2100_Europe"),"lat")
		
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
		
		load("/home/bfrancois/Compound_Event/TCE6/Data/commune_France.RData")
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
	#	  print(i)
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
		      if(min(sqrt_sum_diff_square)<=sqrt(0.3^2+0.3^2)){
		        #print("ok")
			COUNT=COUNT+1
	#		print(COUNT)
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
		    if(LON_France[i]>(-1.6) | LAT_France[j]<46.2 | LAT_France[j]>=49 | count>=26){next}else{
		      coord=c(LON_France[i],LAT_France[j])
	      		subset=which((commune_France_mat[,1]>LON_France[i]-0.2) & (commune_France_mat[,1]<LON_France[i]+0.2) &
			 (commune_France_mat[,2]>LAT_France[j]-0.2) & (commune_France_mat[,2]<LAT_France[j]+0.2),arr.ind=TRUE)
			if(length(subset)<=1){
			next}
		else{
		      diff_square=apply(commune_France_mat[subset,],1,function(x) (x-coord)^2)
		      sum_diff_square=apply(diff_square,2,sum)
		      sqrt_sum_diff_square=sqrt(sum_diff_square)
		 #     print(min(sqrt_sum_diff_square))
		      if(min(sqrt_sum_diff_square)<=sqrt(0.3^2+0.3^2)){
		        if(LON_France[i]>8){next}else{
		          count=count+1
		 #         print("ok")
		          good_tmp_varphy_Bretagne_cut[i,j,] = good_tmp_varphy_France[i,j,]
		        }
		      }
		    }
		  }
		}
		}

		#### convert to small arrays
		coord_lat_bretagne=which(LAT_France >46.2 & LAT_France <=49)
		coord_lon_bretagne=which(LON_France>=(-5.25) & LON_France<=(-1.6))
		LON_Bretagne=LON_France[coord_lon_bretagne]
		LAT_Bretagne=LAT_France[coord_lat_bretagne]
		
		
		good_tmp_varphy_Bretagne_cut=good_tmp_varphy_Bretagne_cut[coord_lon_bretagne,coord_lat_bretagne,]
		
		pdf(paste0("ok",name_varphy,"_", Mod, "_Bretagnenew.pdf"))
		image.plot(LON_Bretagne,LAT_Bretagne,good_tmp_varphy_Bretagne_cut[,,1])
		world(add=TRUE)
		dev.off()
		
		IND_Bretagne=which(!(is.na(good_tmp_varphy_Bretagne_cut)),arr.ind=T)
		
		setwd(paste0("/home/bfrancois/Compound_Event/TCE6/Data/RefMod_Data/Bretagne/"))
		assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_Bretagne"), good_tmp_varphy_Bretagne_cut)
		eval(parse(text=(paste0("save(", name_varphy, "_day_", Mod, "_1850_2100_Bretagne, LON_Bretagne, LAT_Bretagne, IND_Bretagne, file='", name_varphy, "_day_", Mod, "_1850_2100_Bretagne.RData')"))))
		rm(good_tmp_varphy_Bretagne_cut)
		gc()	
		print("ok_final")
#	
		##### Burgundy
		good_tmp_varphy_Burgundy_cut=array(NaN, dim=c(dim(good_tmp_varphy_France)))
		
		#### convert to small arrays
		coord_lat_Burgundy=which(LAT_France >=46 & LAT_France <=49)
		coord_lon_Burgundy=which(LON_France>=(-1) & LON_France<=(5))
		LON_Burgundy=LON_France[coord_lon_Burgundy]
		LAT_Burgundy=LAT_France[coord_lat_Burgundy]
		
		good_tmp_varphy_Burgundy_cut=good_tmp_varphy_France[coord_lon_Burgundy,coord_lat_Burgundy,]
		
		setwd("/scratchu/bfrancois/")
	
		pdf(paste0("ok",name_varphy, "_", Mod,"_Burgundynew.pdf"))
		image.plot(LON_Burgundy,LAT_Burgundy,good_tmp_varphy_Burgundy_cut[,,1])
		dev.off()
		
		IND_Burgundy=which(!(is.na(good_tmp_varphy_Burgundy_cut)),arr.ind=T)
		
		setwd(paste0("/home/bfrancois/Compound_Event/TCE6/Data/RefMod_Data/Burgundy/"))
		assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_Burgundy"), good_tmp_varphy_Burgundy_cut)
		eval(parse(text=(paste0("save(", name_varphy, "_day_", Mod, "_1850_2100_Burgundy, LON_Burgundy, LAT_Burgundy, IND_Burgundy, file='", name_varphy, "_day_", Mod, "_1850_2100_Burgundy.RData')"))))
		rm(good_tmp_varphy_Europe, good_tmp_varphy_France, good_tmp_varphy_France_cut, good_tmp_varphy_Burgundy_cut)
		gc()	
		print("ok_final")

	}
}



print(ok)



	
	#	##### Paris
	#	good_tmp_varphy_Paris_cut=array(NaN, dim=c(dim(good_tmp_varphy_France)))
	#	
	#	count=0
	#	for(i in 1:length(LON_France)){
	#	  print(i)
	#	  for(j in 1:length(LAT_France)){
	#	    if(LON_France[i]<1.1 | LON_France[i]>3.6 | LAT_France[j]<47.8
	#	       | LAT_France[j]>49.8 | count>=80){next}else{
	#	      coord=c(LON_France[i],LAT_France[j])
	#	      diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
	#	      sum_diff_square=apply(diff_square,2,sum)
	#	      sqrt_sum_diff_square=sqrt(sum_diff_square)
	#	      print(min(sqrt_sum_diff_square))
	#	      if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
	#	        if(LON_France[i]>8){next}else{
	#	          count=count+1
	#	          print("ok")
	#	          good_tmp_varphy_Paris_cut[i,j,] = good_tmp_varphy_France[i,j,]
	#	        }
	#	      }
	#	    }
	#	  }
	#	}
	#	
	#	#### convert to small arrays
	#	coord_lat_Paris=which(LAT_France >=47.75 & LAT_France <=50)
	#	coord_lon_Paris=which(LON_France>=1 & LON_France<(4))
	#	LON_Paris=LON_France[coord_lon_Paris]
	#	LAT_Paris=LAT_France[coord_lat_Paris]
	#	
	#	
	#	#tas_day_ERA5_1850_2100_Paris=tas_day_ERA5_1850_2100_Paris[coord_lon_Paris,coord_lat_Paris,]
	#	good_tmp_varphy_Paris_cut=good_tmp_varphy_Paris_cut[coord_lon_Paris,coord_lat_Paris,]
	#	
	#	pdf(paste0("ok",name_varphy,"Parisnew.pdf"))
	#	image.plot(LON_Paris,LAT_Paris,good_tmp_varphy_Paris_cut[,,1])
	#	dev.off()
	#	
	#	IND_Paris=which(!(is.na(good_tmp_varphy_Paris_cut)),arr.ind=T)
	#	
	#	setwd(paste0("/home/bfrancois/Compound_Event/TCE6/Data/RefMod_Data/", Mod, "/")
	#	assign(paste0(name_varphy, "_day_", Mod, "_1850_2100_Paris"), good_tmp_varphy_Paris_cut)
	#	eval(parse(text=(paste0("save(", name_varphy, "_day_", Mod, "_1850_2100_Paris, LON_Paris, LAT_Paris, IND_Paris, file='", name_varphy, "_day_", Mod, "_1850_2100_Paris.RData')"))))
	#	rm(good_tmp_varphy_Europe, good_tmp_varphy_France, good_tmp_varphy_France_cut, good_tmp_varphy_Paris_cut)
	#	gc()	
	#	print("ok_final")
	
#
## # #### Function calendrier Mathieu voir mail cmip6 du 27.04.21
## # MOD = c("BCC-CSM2-MR_r1i1p1f1", "CanESM5_r10i1p1f1", "CESM2_r1i1p1f1",
## #         "CESM2-WACCM_r1i1p1f1", "CNRM-CM6-1-HR_r1i1p1f2",
## #         "CNRM-CM6-1_r1i1p1f2", "CNRM-ESM2-1_r1i1p1f2",
## #         "GFDL-CM4_r1i1p1f1", "INM-CM4-8_r1i1p1f1",  "INM-CM5-0_r1i1p1f1",
## #         "IPSL-CM6A-LR_r14i1p1f1", "MIROC6_r1i1p1f1",
## #         "MRI-ESM2-0_r1i1p1f1", "UKESM1-0-LL_r1i1p1f2")
## # CAL =c("365_day", "365_day", "365_day", "365_day", "standard",
## #        "standard", "standard", "365_day", "365_day", "365_day",
## #        "standard", "standard", "proleptic_gregorian", "360_day")
## # 
## # 
## # 
## # 
## # create_calendar_CMIP = function(namecal, start, end){
## #         library(timeDate)
## #         if(namecal=="365_day"){
## #                 DATES_all =atoms(timeSequence(from=paste(start,"-01-01",sep=""),to=paste(end,"-12-31",sep=""),by='day'))
## #                 cat(dim(DATES_all)," ")
## #                 Ind_29Feb = which(DATES_all[,"m"]==2 & DATES_all[,"d"]==29)
## #                 DATES = DATES_all[-Ind_29Feb,] # remove rows with 29 Feb.
## #                 cat("-->",dim(DATES),"\n")
## #         }
## #         else{
## #                 if(namecal=="standard"){
## #                         DATES =atoms(timeSequence(from=paste(start,"-01-01",sep=""),to=paste(end,"-12-31",sep=""),by='day'))
## #                         cat(dim(DATES),"\n")
## #                 }
## #                 else{
## #                         if(namecal=="proleptic_gregorian"){
## #                                 DATES =atoms(timeSequence(from=paste(start,"-01-01",sep=""),to=paste(end,"-12-31",sep=""),by='day'))
## #                                 cat(dim(DATES),"\n")
## #                         }
## #                         else{
## #                                 if(namecal=="360_day"){
## #                                         DATES = array(NaN,
## #                                                       dim=c((360*(as.numeric(end)-as.numeric(start)+1)),3))
## #                                         colnames(DATES) = c("Y", "m", "d")
## #                                         cpt = 0
## #                                         
## #                                         for(y in (as.numeric(start)):(as.numeric(end))){
## #                                                 for(m in 1:12){
## #                                                         for(d in 1:30){
## #                                                                 cpt = cpt+1
## #                                                                 DATES[cpt,1] = y
## #                                                                 DATES[cpt,2] = m
## #                                                                 DATES[cpt,3] = d
## #                                                         }
## #                                                 }
## #                                         }
## # 
## # cat(dim(DATES),"\n")
## #                                 }
## #                         }
## #                 }
## #         }
## #         return(DATES)
## # }
## 
##



### Brouillon

#/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/sfcWindmax/gr/latest
	### old but good
	#COUNT=0
	#for(i in 1:length(LON_France)){
	#  print(i)
	#  for(j in 1:length(LAT_France)){
	#    if(LON_France[i]>8){next}else{
	#      coord=c(LON_France[i],LAT_France[j])
	#      diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
	#      sum_diff_square=apply(diff_square,2,sum)
	#      sqrt_sum_diff_square=sqrt(sum_diff_square)
	#      #print(min(sqrt_sum_diff_square))
	#      if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
	#          print("ok")
	#          COUNT=COUNT+1
	#	  #good_tmp_varphy_France_cut[i,j,]=good_tmp_varphy_France[i,j,]
	#      }
	#    }
	#  }
	#}
	### end old but good
	

### PR
#system(paste0("cdo mergetime ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -mulc,86400 -del29feb ",
#	      "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/pr/gr/latest/pr_day_CNRM-CM6-1_historical_r1i1p1f2_gr_18500101-19491231.nc ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -mulc,86400 -del29feb ",
#              "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/pr/gr/latest/pr_day_CNRM-CM6-1_historical_r1i1p1f2_gr_19500101-20141231.nc ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -mulc,86400 -del29feb ",
#              "/bdd/CMIP6/ScenarioMIP/CNRM-CERFACS/CNRM-CM6-1/ssp585/r1i1p1f2/day/pr/gr/latest/pr_day_CNRM-CM6-1_ssp585_r1i1p1f2_gr_20150101-21001231.nc ",
#              "/scratchu/bfrancois/pr_day_CNRMCM6_1850_2100_Europe.nc"))
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
#
#
#
#
#
#print(ok3)
######TAS at 0.25deg
#tmp_tas=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_1850_2100_Europe","tas")
#
#
#LON_Europe=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_1850_2100_Europe","longitude")
##### ATTENTION ROTATION LAT and so need to rotate tmp_tas on row axis
#LAT_Europe=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_1850_2100_Europe","latitude")
#
#print(LON_Europe)
#print(LAT_Europe)
##### Rotation
#rotate <- function(x) t(apply(x, 1, rev))
#
#tmp_tas_rot=array(NaN,dim=c(dim(tmp_tas)))
#
##### si LAT in decreasing order
#if(LAT_Europe[1]>=LAT_Europe[2]){
#	LAT_Europe=rev(LAT_Europe)
#	for(i in 1:dim(tmp_tas)[3]){
#	        if(i%%777==1){round(print(i)/dim(tmp_tas)[3]*100,2)}
#	        tmp_tas_rot[,,i]<-rotate(tmp_tas[,,i])
#	}
#}
#
#
##dim(tmp_tas_rot)
##dim(tmp_tas)
#setwd("/scratchu/bfrancois/")
#
##pdf("essai.pdf")
##par(mfrow=c(3,3))
##mean(tmp_tas[,,1])
##mean(tmp_tas_rot[,,1])
##mean(tmp_tas[,,500])
##mean(tmp_tas_rot[,,500])
##mean(tmp_tas[,,10000])
##mean(tmp_tas_rot[,,10000])
##dev.off()
#
###### Preprocess CNRMCM6
#tas_day_CNRMCM6_1850_2100_Europe<-tmp_tas_rot
##pr_day_CNRMCM6_1850_2100_Europe<-tmp_pr_rot
#
##par(mfrow=c(2,2))
##image.plot(LON_Europe,LAT_Europe,apply(tas_day_CNRMCM6_1850_2100_Europe,c(1,2),mean))
##world(add=T)
##image.plot(LON_Europe,LAT_Europe,apply(pr_day_CNRMCM6_1850_2100_Europe,c(1,2),mean))
##world(add=T)
#
#
#
#rm(tmp_tas_rot, tmp_tas)
#gc()
#
#
#
#
#
##load commune information to extract a specific region in France
###obtained from base=read.csv("http://freakonometrics.free.fr/popfr19752010.csv",header=TRUE)
##rm(list=ls())
#
#load("/home/bfrancois/Compound_Event/TCE5/Data/commune_France.RData")
##load("/scratchu/bfrancois/tas_pr_day_CNRMCM6_1850_2100_Europe.RData")
#print(LON_Europe)
#print(LAT_Europe)
#print(dim(tas_day_CNRMCM6_1850_2100_Europe))
##print(dim(pr_day_CNRMCM6_1850_2100_Europe))
#
##### France
#library(maps)
##points(base$long,base$lat,cex=.1,col="red",pch=19)
##points(base$long,base$lat,cex=2*base$pop_2010/
#
##### Define France coordinates
#lon_coord_france=which(LON_Europe>=min(commune_France$long)-0.5
#                             & LON_Europe<=max(commune_France$long) )
#lat_coord_france=which(LAT_Europe>=min(commune_France$lat)& LAT_Europe<=max(commune_France$lat)+0.5)
#
#LON_France=LON_Europe[lon_coord_france]
#LAT_France=LAT_Europe[lat_coord_france]
#
#tas_day_CNRMCM6_1850_2100_France_Big=tas_day_CNRMCM6_1850_2100_Europe[lon_coord_france,lat_coord_france,]
##pr_day_CNRMCM6_1850_2100_France_Big=pr_day_CNRMCM6_1850_2100_Europe[lon_coord_france,lat_coord_france,]
#
#
#old.par <- par(mar = c(0, 0, 0, 0))
#par(old.par)
#
#par(mfrow=c(2,2))
#image.plot(LON_France,LAT_France,tas_day_CNRMCM6_1850_2100_France_Big[,,1])
##image.plot(LON_France,LAT_France,pr_day_CNRMCM6_1850_2100_France_Big[,,1])
#
#
#commune_France_mat=matrix(NaN,ncol=2, nrow=length(commune_France$long))
#commune_France_mat[,1]=commune_France$long
#commune_France_mat[,2]=commune_France$lat
#
#
#tas_day_CNRMCM6_1850_2100_France=array(NaN,dim=c(dim(tas_day_CNRMCM6_1850_2100_France_Big)))
##pr_day_CNRMCM6_1850_2100_France=array(NaN,dim=c(dim(pr_day_CNRMCM6_1850_2100_France_Big)))
#
#
#
#for(i in 1:length(LON_France)){
#  print(i)
#  for(j in 1:length(LAT_France)){
#    if(LON_France[i]>8){next}else{
#      coord=c(LON_France[i],LAT_France[j])
#      diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
#      sum_diff_square=apply(diff_square,2,sum)
#      sqrt_sum_diff_square=sqrt(sum_diff_square)
#      print(min(sqrt_sum_diff_square))
#      if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
#        if(LON_France[i]>8){next}else{
#          print("ok")
#          tas_day_CNRMCM6_1850_2100_France[i,j,]=tas_day_CNRMCM6_1850_2100_France_Big[i,j,]
# #         pr_day_CNRMCM6_1850_2100_France[i,j,]=pr_day_CNRMCM6_1850_2100_France_Big[i,j,]
#        }
#      }
#    }
#  }
#}
#
#setwd("/scratchu/bfrancois/")
#pdf("oktasFrancenew.pdf")
#par(mfrow=c(2,2))
#image.plot(LON_France,LAT_France,tas_day_CNRMCM6_1850_2100_France[,,1])
#world(add=T)
#image.plot(LON_France,LAT_France,tas_day_CNRMCM6_1850_2100_France_Big[,,1])
#world(add=T)
#dev.off()
#
##pdf("okprFrancenew.pdf")
##par(mfrow=c(2,2))
###image.plot(LON_France,LAT_France,pr_day_ERA5_1850_2100_France[,,1])
###world(add=T)
###image.plot(LON_France,LAT_France,pr_day_ERA5_1850_2100_France_Big[,,1])
###world(add=T)
##image.plot(LON_France,LAT_France,pr_day_CNRMCM6_1850_2100_France[,,1])
##world(add=T)
##image.plot(LON_France,LAT_France,pr_day_CNRMCM6_1850_2100_France_Big[,,1])
##world(add=T)
##dev.off()
#
#
#
#IND_France=which(!(is.na(tas_day_CNRMCM6_1850_2100_France)),arr.ind=T)
#
##setwd("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/RefMod_Data/ERA5/")
##save(tas_day_ERA5_1850_2100_France,pr_day_ERA5_1850_2100_France,LON_France,LAT_France,IND_France,
##                              file="tas_pr_day_ERA5_1850_2100_France.RData")
#setwd("/scratchu/bfrancois/")
#save(tas_day_CNRMCM6_1850_2100_France,
#	#pr_day_CNRMCM6_1850_2100_France,
#	LON_France,LAT_France,IND_France,
#                              file="tas_day_CNRMCM6_1850_2100_France.RData")
#
#
##
###### Bretagne
#
#tas_day_CNRMCM6_1850_2100_Bretagne=array(NaN,dim=c(dim(tas_day_CNRMCM6_1850_2100_France_Big)))
##pr_day_CNRMCM6_1850_2100_Bretagne=array(NaN,dim=c(dim(pr_day_CNRMCM6_1850_2100_France_Big)))
#
#count=0
#for(i in 1:length(LON_France)){
#  print(i)
#  for(j in 1:length(LAT_France)){
#    if(LON_France[i]>0 | LAT_France[j]<46 | LAT_France[j]>=49 | count>=80){next}else{
#      coord=c(LON_France[i],LAT_France[j])
#      diff_square=apply(commune_France_mat,1,function(x) (x-coord)^2)
#      sum_diff_square=apply(diff_square,2,sum)
#      sqrt_sum_diff_square=sqrt(sum_diff_square)
#      print(min(sqrt_sum_diff_square))
#      if(min(sqrt_sum_diff_square)<=sqrt(0.1^2+0.1^2)){
#        if(LON_France[i]>8){next}else{
#          count=count+1
#          print("ok")
#          tas_day_CNRMCM6_1850_2100_Bretagne[i,j,]=tas_day_CNRMCM6_1850_2100_France_Big[i,j,]
##          pr_day_CNRMCM6_1850_2100_Bretagne[i,j,]=pr_day_CNRMCM6_1850_2100_France_Big[i,j,]
#        }
#      }
#    }
#  }
#}
#
#par(mfrow=c(2,2))
#image.plot(LON_France,LAT_France,tas_day_CNRMCM6_1850_2100_Bretagne[,,1])
#world(add=T)
#image.plot(LON_France,LAT_France,tas_day_CNRMCM6_1850_2100_France_Big[,,1])
#world(add=T)
#
#
##par(mfrow=c(2,2))
##image.plot(LON_France,LAT_France,pr_day_CNRMCM6_1850_2100_Bretagne[,,1])
##world(add=T)
##image.plot(LON_France,LAT_France,pr_day_CNRMCM6_1850_2100_France_Big[,,1])
##world(add=T)
#
##### convert to small arrays
#coord_lat_bretagne=which(LAT_France >46 & LAT_France <=49)
#coord_lon_bretagne=which(LON_France>=(-5.25) & LON_France<=(-1.5))
#LON_Bretagne=LON_France[coord_lon_bretagne]
#LAT_Bretagne=LAT_France[coord_lat_bretagne]
#
#
##tas_day_ERA5_1850_2100_Bretagne=tas_day_ERA5_1850_2100_Bretagne[coord_lon_bretagne,coord_lat_bretagne,]
#tas_day_CNRMCM6_1850_2100_Bretagne=tas_day_CNRMCM6_1850_2100_Bretagne[coord_lon_bretagne,coord_lat_bretagne,]
##pr_day_ERA5_1850_2100_Bretagne=pr_day_ERA5_1850_2100_Bretagne[coord_lon_bretagne,coord_lat_bretagne,]
##pr_day_CNRMCM6_1850_2100_Bretagne=pr_day_CNRMCM6_1850_2100_Bretagne[coord_lon_bretagne,coord_lat_bretagne,]
#
#
#
#
##setwd("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/")
#pdf("oktasBretagnenew.pdf")
#image.plot(LON_Bretagne,LAT_Bretagne,tas_day_ERA5_1850_2100_Bretagne[,,1])
#dev.off()
#
#
#IND_Bretagne=which(!(is.na(tas_day_CNRMCM6_1850_2100_Bretagne)),arr.ind=T)
#
##setwd("/home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/RefMod_Data/ERA5/")
##save(tas_day_ERA5_1850_2100_Bretagne,pr_day_ERA5_1850_2100_Bretagne,LON_Bretagne,LAT_Bretagne,IND_Bretagne,
##                              file="tas_pr_day_ERA5_1850_2100_Bretagne.RData")
#setwd("/home/bfrancois/Compound_Event/TCE5/Data/RefMod_Data/CNRMCM6/")
#save(tas_day_CNRMCM6_1850_2100_Bretagne,
#	#pr_day_CNRMCM6_1850_2100_Bretagne,
#	LON_Bretagne,LAT_Bretagne,IND_Bretagne,
#                              file="tas_day_CNRMCM6_1850_2100_Bretagne.RData")
#
#print("ok_final")
#print(ok)
#
#
#
#
#
#
#tmp_pr=convert_netcdf_RData("/scratchu/bfrancois/pr_day_CNRMCM6_1850_2100_Europe","pr")
#tmp_pr_rot=array(NaN,dim=c(dim(tmp_pr)))
#
#for(i in 1:dim(tmp_pr)[3]){
#        if(i%%777==1){round(print(i)/dim(tmp_pr)[3]*100,2)}
#        #tmp_tas_rot[,,i]<-rotate(tmp_tas[,,i])
#        tmp_pr_rot[,,i]<-rotate(tmp_pr[,,i])
#}
#
#pr_day_CNRMCM6_1850_2100_Europe<-tmp_pr_rot
#
#rm(tmp_pr_rot, tmp_pr)
#
##setwd("/scratchu/bfrancois/")
##save(tas_day_CNRMCM6_1850_2100_Europe,
##     pr_day_CNRMCM6_1850_2100_Europe,
##     LON_Europe,
##     LAT_Europe,
##     file="tas_pr_day_CNRMCM6_1850_2100_Europe.RData")
#
#print("ok_save_scratch")
#
#
#
#
#### sfcWindmax
##/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/sfcWindmax/gr/latest
#
##sfcWindmax_day_CNRM-CM6-1_historical_r1i1p1f2_gr_18500101-19491231.nc
#### sfcWindmax
#system(paste0("cdo mergetime ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -del29feb ",
#	      "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/sfcWindmax/gr/latest/sfcWindmax_day_CNRM-CM6-1_historical_r1i1p1f2_gr_18500101-19491231.nc ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -del29feb ",
#              "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/sfcWindmax/gr/latest/sfcWindmax_day_CNRM-CM6-1_historical_r1i1p1f2_gr_19500101-20141231.nc ",
#              "-remapbil,", Grid_ERA5_0_25 ,
#              " -selyear,1850/2100 -del29feb ",
#              "/bdd/CMIP6/ScenarioMIP/CNRM-CERFACS/CNRM-CM6-1/ssp585/r1i1p1f2/day/sfcWindmax/gr/latest/sfcWindmax_day_CNRM-CM6-1_ssp585_r1i1p1f2_gr_20150101-21001231.nc ",
#              "/scratchu/bfrancois/sfcWindmax_day_CNRMCM6_1850_2100_Europe.nc"))
#
#
##
##
##
#
####CNRMCM6 at 1deg
###system(paste0("cdo -seltimestep,1/10 /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/t2m_era5_5020_1deg.nc* /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/target_grid_1deg.nc"))
####
#######CNRM-CMIP6
######## To launch from Ciclad: /home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/
###
###Grid_ERA5_1=paste("/home/bfrancois/1_Data/Vrac/dependencies_CC/ERA5/target_grid_1deg.nc")
######## TAS
###system(paste0("cdo mergetime ",
###        "-remapbil,", Grid_ERA5_1 ,
###        " -selyear,1951/2020 -addc,-273.15 -del29feb ",
###        "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/tas/gr/latest/tas_day_CNRM-CM6-1_historical_r1i1p1f2_gr_18500101-20141231.nc ",
###        "-remapbil,", Grid_ERA5_1 ,
###        " -selyear,1951/2020 -addc,-273.15 -del29feb ",
###        "/bdd/CMIP6/ScenarioMIP/CNRM-CERFACS/CNRM-CM6-1/ssp585/r1i1p1f2/day/tas/gr/latest/tas_day_CNRM-CM6-1_ssp585_r1i1p1f2_gr_20150101-21001231.nc ",
###        "/scratchu/bfrancois/tas_day_CNRMCM6_51_20_Europe_1deg.nc"))
####
####
###### PR
###system(paste0("cdo mergetime ",
###              "-remapbil,", Grid_ERA5_1 ,
###              " -selyear,1951/2020 -mulc,86400 -del29feb ",
###              "/bdd/CMIP6/CMIP/CNRM-CERFACS/CNRM-CM6-1/historical/r1i1p1f2/day/pr/gr/latest/pr_day_CNRM-CM6-1_historical_r1i1p1f2_gr_19500101-20141231.nc ",
###              "-remapbil,", Grid_ERA5_1 ,
###              " -selyear,1951/2020 -mulc,86400 -del29feb ",
###              "/bdd/CMIP6/ScenarioMIP/CNRM-CERFACS/CNRM-CM6-1/ssp585/r1i1p1f2/day/pr/gr/latest/pr_day_CNRM-CM6-1_ssp585_r1i1p1f2_gr_20150101-21001231.nc ",
###              "/scratchu/bfrancois/pr_day_CNRMCM6_51_20_Europe_1deg.nc"))
####
### 
######TAS at 1deg
###tmp_tas=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_51_20_Europe_1deg","tas")
###tmp_pr=convert_netcdf_RData("/scratchu/bfrancois/pr_day_CNRMCM6_51_20_Europe_1deg","pr")
###
###
###
###LON_Europe_1deg=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_51_20_Europe_1deg","lon")
###LAT_Europe_1deg=convert_netcdf_RData("/scratchu/bfrancois/tas_day_CNRMCM6_51_20_Europe_1deg","lat")
###
###print(LON_Europe_1deg)
###print(LAT_Europe_1deg)
###
###dim(tmp_tas)
###dim(tmp_pr)
###
###
####### Preprocess CNRMCM6
###tas_day_CNRMCM6_51_20_Europe_1deg<-tmp_tas
###pr_day_CNRMCM6_51_20_Europe_1deg<-tmp_pr
###
###image.plot(LON_Europe_1deg,LAT_Europe_1deg,apply(tas_day_CNRMCM6_51_20_Europe_1deg,c(1,2),mean))
###world(add=T)
###image.plot(LON_Europe_1deg,LAT_Europe_1deg,apply(pr_day_CNRMCM6_51_20_Europe_1deg,c(1,2),mean))
###world(add=T)
###
###setwd("/scratchu/bfrancois/")
###save(tas_day_CNRMCM6_51_20_Europe_1deg,
###     pr_day_CNRMCM6_51_20_Europe_1deg,
###     LON_Europe_1deg,
###     LAT_Europe_1deg,
###     file="tas_pr_day_CNRMCM6_51_20_Europe_1deg.RData")
###
##
#
##
## 
## ### Copy the file from Ciclad to Obelix
## ### TAS and PR
## # scp -r bfrancois@ciclad.ipsl.jussieu.fr:/scratchu/bfrancois/*2100* /home/starmip/bfran/LSCE_These/Compound_Event/TCE5/Data/RefMod_Data/CNRMCM6/
##
#
#

