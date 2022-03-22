import pandas as pd

path = "/home/milo/PCIC/Maestría/2doSemestre/seminario/github/data/experimental_design/"
regiones = ["america","eurafrica","asia"]
policies = ["P"+str(p) for p in range(8) for _ in range(400)]
parameter_set = [str(i+1 )for _ in range(8) for i in range(400)]

gcms = ["CESM2" ,"GFDL"]

df_total = pd.DataFrame()

for gcm in gcms:
    for region in regiones:
        csv_file_name = "cmpi6_{}_{}_experimental_design.csv".format(gcm,region)
        df_gcm_region = pd.DataFrame()
        print("Load data by GCM: {} and Region: {}".format(gcm,region))
        for policy in range(8):
            print("Policy {}".format(policy))
            df_gcm_region = pd.concat([df_gcm_region,pd.read_csv(path+csv_file_name)])
        df_gcm_region["policy"] = policies
        df_gcm_region["parameter_set"] = parameter_set
        df_gcm_region["region"] = [region for _ in range(400*8)]
        df_gcm_region["gcm"] = [gcm for _ in range(400*8)]

        # Create key
        df_gcm_region["key"] = df_gcm_region["gcm"] + "-" + df_gcm_region["region"] + "-" + df_gcm_region["parameter_set"] + "-" +  df_gcm_region["policy"]


        df_total = pd.concat([df_total,df_gcm_region])

# Load results data
results_df = pd.read_csv(path+"ediam_regions_results.csv")
results_df["key"] = results_df["gcm"].apply(lambda x : str(x)) + "-" + results_df["region"].apply(lambda x : str(x).replace("_oceania","")) + "-" + results_df["parameter_set"].apply(lambda x : str(x)) + "-" + results_df["policy"].apply(lambda x : str(x))
results_df = results_df[["key","resultados_2_c","resultados_3_c"]]

df_total_results = pd.merge(df_total,results_df,on = "key")
df_total_results.query("resultados_2_c !='RNF'",inplace = True)
df_total_results.query("resultados_3_c !='RNF'",inplace = True)

df_total_results["goal_2_c"] = 0
df_total_results["goal_2_c"][df_total_results["resultados_2_c"] == "true"] = 1

df_total_results["goal_3_c"] = 0
df_total_results["goal_3_c"][df_total_results["resultados_3_c"] == "true"] = 1
df_total_results.drop(columns=["Unnamed: 0"],inplace=True)
df_total_results.to_csv("data_comp_exp_ML.csv",index=False)
