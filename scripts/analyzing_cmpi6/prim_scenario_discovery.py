from __future__ import print_function

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import prim

# Import dataset
df_all = pd.read_csv("data_comp_exp_ML.csv")

# Check for missing values
df_all.isnull().sum()

region = "asia"
df = df_all.query("region =='{}'".format(region))
# Declare feature vector and target variable
X = df.drop(columns=["resultados_2_c","resultados_3_c","region","policy","parameter_set","gcm","key","goal_2_c","goal_3_c"])
response = df["goal_2_c"]

p = prim.Prim(X,
              response,
              threshold=1,
              threshold_type="=")
box = p.find_box()
box.show_tradeoff()
plt.show()
