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
var_cmip = ["co2","co2mass","tas"]

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
                print("Abriendo los nc de la ruta\n{}/{}/{}/*.nc".format(gcm_path,rcp,var))
                ds = xr.open_mfdataset("{}/{}/{}/*.nc".format(gcm_path,rcp,var))
                ds = ds[var].groupby('time.year').mean()
                year_init,year_final =  int(ds.year[0].values), int(ds.year[-1].values)
                ds.to_netcdf("{}_Amon_{}_{}_{}_{}.nc".format(var,gcm,rcp,year_init,year_final))
            except:
                print("No hay archivos NetCFD en el directorio")


    for var in var_cmip:
        try:
            r = re.compile("{}_".format(var))
            lista_filtrada = list(filter(r.match, glob.glob('*.nc')))
            lista_filtrada.sort()

            historical = lista_filtrada[0]
            lista_filtrada.remove(historical)
            historical_xr = xr.open_dataset(historical, chunks={'year': 120})

            for rcp in lista_filtrada:
                print("Variable {}\nModelo {}".format(var,rcp))
                rcp_xr = xr.open_dataset(rcp, chunks={'year': 120})
                historical_xr_var = historical_xr[var]
                rcp_xr_var = rcp_xr[var]
                historical_rcp_xr = xr.merge([historical_xr_var,rcp_xr_var])

                for region in macroregiones:

                    if var == 'tas':
                        historical_rcp_xr_sub = historical_rcp_xr.isel(year=slice(0,240))
                        historical_rcp_xr_sub_c = historical_rcp_xr_sub[var] - 273.15
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
                            historical_rcp_mean_values = historical_rcp_xr[var].isel(lon=slice(77, 135)).mean(dim=("plev","lat","lon")).values
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
                                                         'scenario':[rcp.split("_")[3] for i in range(n)],
                                                         'region' : [region for i in range(n)],
                                                         'variable':[var for i in range(n)],
                                                         'value':historical_rcp_mean_values})

                    agrupa_todo = pd.concat([agrupa_todo,agrupa_rcp])
        except :
            print("No hay archivos NetCFD para el RCP: {}".format(rcp))




    os.chdir("..")
agrupa_todo.to_csv("gcm_year_cmip6.csv",index=False)
