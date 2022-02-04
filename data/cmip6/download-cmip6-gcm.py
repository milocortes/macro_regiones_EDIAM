from pyesgf.search import SearchConnection
from tqdm import tqdm
import pandas as pd

df_result_total = pd.DataFrame()


gcm_models = ['BCC-CSM2-MR','BCC-EM1','CanESM5','CESM2','CESM2-WACCM','CNRM-CM6-1','CNRM-ESM2-1','EC-EARTH3-Veg','GFDL-ESM4','GFDL-CM4','IPSL-CM6A-LR','MRI-ESM2-0','NorESM2-LM','SAM0-UNICON','UKESM1-0-LL']

cmip_url = ["https://esgf-node.llnl.gov/esg-search",
        "https://esgf-node.ipsl.upmc.fr/esg-search",
        "https://esgf-data.dkrz.de/esg-search",
        "https://esgf-index1.ceda.ac.uk/esg-search"]

rcps = ['historical','ssp126','ssp245','ssp370','ssp585']

climate_variables = ['co2mass','co2','tas']

for url_node in cmip_url:
    print("############################################################3")
    print("Búsqueda en la liga {}".format(url_node))
    conn = SearchConnection(url_node, distrib=True)
    print("############################################################3")

    for model in tqdm(gcm_models):
        for rcp in rcps:
            for clim_var in climate_variables:
                ctx = conn.new_context(
                    project='CMIP6',
                    frequency='mon',
                    source_id=model,
                    experiment_id=rcp,
                    variant_label='r1i1p1f1',
                    variable=clim_var)

                if ctx.hit_count > 0:

                    df_result_parcial = pd.DataFrame().from_dict({"climate_model" : [model],
                                                                    "rcp" : [rcp],
                                                                    "climate_variable":[clim_var]})
                    print("Nodo: {}\nModelo: {}\nRCP: {}\nVariable climática: {}".format(url_node,model,rcp,clim_var))
                    for i in tqdm(range(len(ctx.search()))):
                        try:
                            result = ctx.search()[i]
                            files = result.file_context().search()

                            for file in files:
                                #print("Nodo: {}\nModelo: {}\nRCP: {}\nVariable climática: {}\nURL: {}".format(url_node,model,rcp,clim_var,file.opendap_url))
                                #df_result_parcial["url"] = [file.opendap_url]
                                df_result_parcial["url"] = [file.download_url]
                                df_result_total = pd.concat([df_result_total,df_result_parcial])
                        except Exception as e:
                            print("Hubo un error")

df_result_total.to_csv('variant_label_r1i1p1f1.csv')

# Hacemos un subset de los modelos que cumplen con todos los datos necesarios
final_models = ['CESM2-WACCM','NorESM2-LM','GFDL-ESM4']

df_final_models = df_result_total[df_result_total["climate_model"].isin(final_models)]

# Generamos una columna de fechas
df_final_models["date"] = pd.to_datetime(df_final_models["url"].apply(lambda x: x.split("/")[-2][1:]), format='%Y%m%d')
# Generamos una columna de nodo
df_final_models["nodo"] = df_final_models["url"].apply(lambda x: x.split("/")[2])

# Generamos una columna de la perioricidad del dato
df_final_models["temporalidad"] = df_final_models["url"].apply(lambda x: x.split("/")[-1].split("_")[1])

# Nos quedamos con las últimas versiones de cada variable, rcp y modelo
final_models_last_date = []

for model in final_models:
    for nodo in set(df_final_models["nodo"] ):
        for rcp in rcps:
            for cv in climate_variables:
                parcial_final_models = df_final_models.query("climate_model == '{}' and nodo =='{}' and climate_variable=='{}' and rcp == '{}' and temporalidad=='Amon'".format(model,nodo,cv,rcp))

                urls = set(list(parcial_final_models[parcial_final_models["date"]== parcial_final_models["date"].min()]["url"]))

                for url in urls:
                    final_models_last_date.append({'climate_model' : model,
                                                   'nodo' : nodo,
                                                   'rcp' : rcp,
                                                   'climate_variable' : cv,
                                                   'url' : url})

df_final_models_last_date = pd.DataFrame(final_models_last_date)

for model in final_models:
    for nodo in set(df_final_models_last_date["nodo"]):
        for cv in climate_variables:
            consulta = "climate_model=='{}' and nodo=='{}' and climate_variable=='{}'".format(model,nodo,cv)
            if set(df_final_models_last_date.query(consulta)["rcp"]):
                if len(set(df_final_models_last_date.query(consulta)["rcp"]))==5:
                    print("%------------------------------------------------------%")
                    print(consulta)
                    print(set(df_final_models_last_date.query(consulta)["rcp"]))
                    print(df_final_models_last_date.query(consulta).shape)


'''
    Usaremos los siguientes nodos para los modelos:
        * CESM2-WACCM ----> esgf-data.ucar.edu
        * NorESM2-LM  ----> noresg.nird.sigma2.no
        * GFDL-ESM4   ----> esgdata.gfdl.noaa.gov
'''

def get_url_final_models(model,nodo):
    parcial = df_final_models_last_date.query("climate_model=='{}' and nodo =='{}'".format(model,nodo))

    for rcp in rcps:
        f = open("{}_{}.txt".format(model,rcp),"x")
        for url in parcial.query("rcp=='{}'".format(rcp))["url"]:
            f.write("{}\n".format(url))

        f.close()

get_url_final_models("CESM2-WACCM","esgf-data.ucar.edu")
get_url_final_models("NorESM2-LM","noresg.nird.sigma2.no")
get_url_final_models("GFDL-ESM4","esgdata.gfdl.noaa.gov")


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import os



co2 = Dataset("co2_Amon_MRI-ESM1_rcp85_r1i1p1_200601-201512.nc",mode = "r")
np_co2 = np.array(co2["co2"])
np_co2[0].mean(axis=0)


tas = Dataset("tas_Amon_MRI-ESM1_rcp85_r1i1p1_200601-210012.nc",mode="r")
np_tas = np.array(tas["tas"])


for model in set(df_result_total["climate_model"]):
     for rcp in set(df_result_total["rcp"]):
         for clim_var in set(df_result_total["climate_variable"]):
             dimension = df_result_total.query("climate_model=='{}' and rcp=='{}' and climate_variable =='{}'".format(model,rcp,clim_var)).shape
    ...:             print("Modelo: {}. RCP: {}. ClimVar: {}. Lenght: {}".format(model,rcp,clim_var,dimension[0]))
