import xarray as xr
import pandas as pd
import re
import glob
import os

# Cambiamos de directorio
os.chdir("../../data/cmip6")

# Lista de GCM
gcms = ["CESM2-WACCM","NorESM2-LM","GFDL-ESM4"]

# Lista de variables
var_cmip = ["co2mass","tas"]

# Lista de macroregiones
macroregiones = ["america","eurafrica","asia"]
agrupa_todo = pd.DataFrame()

for gcm in gcms:
    print("********************************************************")
    print("********************************************************")
    print("%%%%%%%%%%%%%%%%%       {}          %%%%%%%%%%%%%%%%%%%".format(gcm))
    print("********************************************************")
    print("********************************************************")

    # Nos cambiamos al directorio donde se encuetran los datos de los gcm
    os.chdir(gcm)
    # Listamos los directorios dentro de la ruta
    rcps = os.listdir()
    gcm_path = os.getcwd()

    # Por cada directorio listamos los archivos
    for rcp in rcps:
        for var in var_cmip:
            # Agrupamos los NetCFD individuales
            try:
                print("Abriendo los nc historicos de la variable climática \nRCP: {}--- Clim Var: {}/*.nc".format(rcp,var))
                historical = xr.open_mfdataset("{}/{}/{}/*.nc".format(gcm_path,"historical",var))
                historical_xr_var = historical[var].groupby('time.year').mean()
                del historical
                print("Abriendo los nc de la ruta\n{}/{}/{}/*.nc".format(gcm_path,rcp,var))
                rcp_xr = xr.open_mfdataset("{}/{}/{}/*.nc".format(gcm_path,rcp,var))
                rcp_xr_var = rcp_xr[var].groupby('time.year').mean()
                del rcp_xr
                historical_rcp_xr = xr.merge([historical_xr_var,rcp_xr_var])
                del rcp_xr_var
                del historical_xr_var
                print("Para la variable  climática {} Calculamos para cada region".format(var))
                for region in macroregiones:
                    print("Region : {}".format(region))
                    if var == 'tas':
                        historical_rcp_xr_sub = historical_rcp_xr.isel(year=slice(0,240))
                        historical_rcp_xr_sub_c = historical_rcp_xr_sub[var] - 273.15
                        del historical_rcp_xr_sub
                        # Calculamos las anomalías
                        if region == "america":
                            historical_rcp_mean_values = historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >190) & (historical_rcp_xr_sub_c.lon <340)).mean(dim=("lat","lon")).values - historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >190) & (historical_rcp_xr_sub_c.lon <340)).mean(dim=("lat","lon","year")).values
                            anios = historical_rcp_xr_sub[var].year.values
                        elif region == "eurafrica":
                            historical_rcp_mean_values = historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >330) | (historical_rcp_xr_sub_c.lon <60)).mean(dim=("lat","lon")).values - historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >330) | (historical_rcp_xr_sub_c.lon <60)).mean(dim=("lat","lon","year")).values
                            anios = historical_rcp_xr_sub[var].year.values
                        else:
                            historical_rcp_mean_values = historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >60) & (historical_rcp_xr_sub_c.lon <190)).mean(dim=("lat","lon")).values - historical_rcp_xr_sub_c.where( (historical_rcp_xr_sub_c.lon >60) & (historical_rcp_xr_sub_c.lon <190)).mean(dim=("lat","lon","year")).values
                            anios = historical_rcp_xr_sub[var].year.values
                    elif var =="co2":
                        if region=="america":
                            historical_rcp_mean_values = historical_rcp_xr[var].where( (historical_rcp_xr_sub_c.lon >190) & (historical_rcp_xr_sub_c.lon <340)).mean(dim=("plev","lat","lon")).values
                            anios = historical_rcp_xr[var].year.values
                            historical_rcp_mean_values = [v*1000000 for v in historical_rcp_mean_values]
                        elif region =="eurafrica":
                            historical_rcp_mean_values = historical_rcp_xr[var].where( (historical_rcp_xr.lon >330) | (historical_rcp_xr.lon <60)).mean(dim=("plev","lat","lon")).values
                            anios = historical_rcp_xr[var].year.values
                            historical_rcp_mean_values = [v*1000000 for v in historical_rcp_mean_values]
                        else:
                            historical_rcp_mean_values = historical_rcp_xr[var].where( (historical_rcp_xr.lon >60) & (historical_rcp_xr.lon <190)).mean(dim=("plev","lat","lon")).values
                            anios = historical_rcp_xr[var].year.values
                            historical_rcp_mean_values = [v*1000000 for v in historical_rcp_mean_values]

                    else:
                        historical_rcp_mean_values = historical_rcp_xr[var].values
                        anios = historical_rcp_xr[var].year.values

                    n = len(anios)
                    agrupa_rcp = pd.DataFrame.from_dict({'year':anios,
                                                         'climate_model':[gcm for i in range(n)],
                                                         'scenario':[rcp for i in range(n)],
                                                         'region' : [region for i in range(n)],
                                                         'variable':[var for i in range(n)],
                                                         'value':historical_rcp_mean_values})

                    agrupa_todo = pd.concat([agrupa_todo,agrupa_rcp])

            except:
                print("No hay archivos NetCFD en el directorio")

    os.chdir("..")
agrupa_todo.to_csv("gcm_year_cmip6.csv",index=False)
