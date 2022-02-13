import xarray as xr
import pandas as pd
import re
import glob
import os
from termcolor import colored

# Cambiamos de directorio
os.chdir("../../data/cmip6")

# Lista de GCM
gcms = ["CESM2-WACCM","NorESM2-LM","GFDL-ESM4"]

# Lista de variables
var = "co2"

# Lista de macroregiones
macroregiones = ["america","eurafrica","asia"]
agrupa_todo = pd.DataFrame()
root_path = os.getcwd()

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
    rcps.remove("historical")
    gcm_path = os.getcwd()
    # Por cada directorio listamos los archivos

    for region in macroregiones:

        print(colored("**********************************************************************", 'red'))
        print(colored("GCM ----> {} REGION-----> {}".format(gcm,region), 'yellow'))
        print(colored("**********************************************************************", 'red'))

        # Nos cambiamos al directorio donde se encuetran los datos de CO2 historicos
        os.chdir("{}/{}/co2".format(gcm_path,"historical"))
        archivos_historicos = glob.glob("*.nc")
        archivos_historicos.sort()

        acumula_historical_rcp_mean_values = []
        acumula_anios = []

        for archivo in archivos_historicos:
            print("Abriendo el archivo nc historico ---> {}".format(archivo))
            try:
                historical = xr.open_mfdataset(archivo)
                historical_rcp_var = historical.groupby('time.year').mean()

                if region=="america":
                    historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >190) & (historical_rcp_var.lon <340)).mean(dim=("plev","lat","lon")).values
                elif region =="eurafrica":
                    historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >330) | (historical_rcp_var.lon <60)).mean(dim=("plev","lat","lon")).values
                else:
                    historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >60) & (historical_rcp_var.lon <190)).mean(dim=("plev","lat","lon")).values

                anios = list(historical_rcp_var[var].year.values)
                historical_rcp_mean_values = [v*1000000 for v in historical_rcp_mean_values]

                print("anio ---> {}  valor ---> {}".format(anios,historical_rcp_mean_values))
                acumula_historical_rcp_mean_values += historical_rcp_mean_values
                acumula_anios += anios
            except:
                print(colored("No está el archivo NetCFD {}".format(archivo), 'red'))

        # Salimos del directorio de archivos históricos
        os.chdir("../..")

        dict_rcp_values = {'america' : {},
                           'eurafrica' : {},
                           'asia' : {}}
        dict_rcp_anios = {'america' : {},
                           'eurafrica' : {},
                           'asia' : {}}

        for rcp in rcps:
            # Nos cambiamos al directorio donde se encuetran los datos de CO2 del rcp
            os.chdir("{}/co2".format(rcp))
            archivos_rcp = glob.glob("*.nc")
            archivos_rcp.sort()

            dict_rcp_values[region][rcp] = acumula_historical_rcp_mean_values
            dict_rcp_anios[region][rcp] = acumula_anios

            for archivo in archivos_rcp:
                print("Abriendo el archivo nc del rcp {} ---> {}".format(rcp,archivo))
                try:
                    historical = xr.open_mfdataset(archivo)
                    historical_rcp_var = historical.groupby('time.year').mean()

                    if region=="america":
                        historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >190) & (historical_rcp_var.lon <340)).mean(dim=("plev","lat","lon")).values
                    elif region =="eurafrica":
                        historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >330) | (historical_rcp_var.lon <60)).mean(dim=("plev","lat","lon")).values
                    else:
                        historical_rcp_mean_values = historical_rcp_var[var].where( (historical_rcp_var.lon >60) & (historical_rcp_var.lon <190)).mean(dim=("plev","lat","lon")).values

                    anios = list(historical_rcp_var[var].year.values)
                    historical_rcp_mean_values = [v*1000000 for v in historical_rcp_mean_values]

                    print("anio ---> {}  valor ---> {}".format(anios,historical_rcp_mean_values))

                    dict_rcp_values[region][rcp] += historical_rcp_mean_values
                    dict_rcp_anios[region][rcp] += anios
                except:
                    print(colored("No está el archivo NetCFD {}. PATH: ".format(archivo, os.getcwd()), 'red'))



            n = len(dict_rcp_anios[region][rcp])

            print(colored("**********************************************************************", 'red'))
            print(colored("Guardamos datos del GCM ---> {} RCP ---> {}".format(gcm,rcp), 'yellow'))
            print(colored("N ---> {}".format(n), 'yellow'))            
            print(colored("**********************************************************************", 'red'))

            agrupa_rcp = pd.DataFrame.from_dict({'year' : dict_rcp_anios[region][rcp],
                                                 'climate_model' : [gcm for i in range(n)],
                                                 'scenario' : [rcp for i in range(n)],
                                                 'region' : [region for i in range(n)],
                                                 'variable' : [var for i in range(n)],
                                                 'value': dict_rcp_values[region][rcp]})

            agrupa_todo = pd.concat([agrupa_todo,agrupa_rcp])

            # Salimos del directorio de archivos del rcp
            os.chdir("../..")

    # Nos cambiamos al directorio donde se encuetran los datos de los gcm
    os.chdir(root_path)
agrupa_todo.to_csv("co2_gcm_year_cmip6.csv",index=False)
