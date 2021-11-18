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
Mod="IPSL"
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
        }
        if(Mod=="IPSL"){
                Mod_folder="IPSL"
                Mod_name="IPSL-CM6A-LR"
                Run_name="r1i1p1f1"
        }
        path_historical=paste0("/bdd/CMIP6/CMIP/", Mod_folder, "/", Mod_name, "/historical/", Run_name, "/day/", name_varphy, "/gr/latest/")
        path_ssp585=paste0("/bdd/CMIP6/ScenarioMIP/", Mod_folder, "/", Mod_name, "/ssp585/", Run_name, "/day/", name_varphy, "/gr/latest/")

