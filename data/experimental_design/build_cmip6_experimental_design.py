import pandas as pd

# Cargamos parámetros climáticos del cmip6
cmip6_climate = pd.read_csv("cmip6_climate_param.csv")
cmip6_climate_dict = cmip6_climate.to_dict('index')

# Cargamos los parámetros económicos

economic_params = pd.read_csv("parameters.csv")


regiones = ["america","eurafrica","asia","america","eurafrica","asia"]

n = economic_params.shape[0]

for i in range(6):
    df_region = economic_params
    df_region['β_T'] = [cmip6_climate_dict[i]['beta_delta_temp'] for _ in range(n)]
    df_region['CO2_base'] = [cmip6_climate_dict[i]['co2_base'] for _ in range(n)]
    df_region['CO2_Disaster'] = [cmip6_climate_dict[i]['co2_disaster'] for _ in range(n)]
    df_region['qsi'] = [cmip6_climate_dict[i]['qsi'] for _ in range(n)]
    df_region['δ_S'] = [cmip6_climate_dict[i]['delta_s'] for _ in range(n)]
    df_region['S_0'] = [cmip6_climate_dict[i]['s_0'] for _ in range(n)]
    df_region['Yre_N_0'] = [cmip6_climate_dict[i]['Yre_N_0'] for _ in range(n)]
    df_region['Yce_N_0'] = [cmip6_climate_dict[i]['Yce_N_0'] for _ in range(n)]
    df_region['Yre_S_0'] = [cmip6_climate_dict[i]['Yre_S_0'] for _ in range(n)]
    df_region['Yce_S_0'] = [cmip6_climate_dict[i]['Yce_S_0'] for _ in range(n)]
    df_region['region'] = [regiones[i] for _ in range(n)]
    df_region.to_csv("cmpi6_{}_{}_experimental_design.csv".format(cmip6_climate_dict[i]['climate_model'].split("-")[0],regiones[i]))

    
