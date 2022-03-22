import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import warnings

warnings.filterwarnings('ignore')

# Import dataset
df_all = pd.read_csv("data_comp_exp_ML.csv")

# Check for missing values
df_all.isnull().sum()

region = "asia"
df = df_all.query("region =='{}'".format(region))
# Declare feature vector and target variable
X = df.drop(columns=["resultados_2_c","resultados_3_c","region","policy","parameter_set","gcm","key","goal_2_c","goal_3_c"])
y = df["goal_2_c"]

# import XGBoost
import xgboost as xgb

# define data_dmatrix
data_dmatrix = xgb.DMatrix(data=X,label=y)

# split X and y into training and testing sets
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)

# Train the XGBoost classifier

# import XGBClassifier
from xgboost import XGBClassifier
# declare parameters
params = {
            'objective':'binary:logistic',
            'max_depth': 4,
            'alpha': 10,
            'learning_rate': 1.0,
            'n_estimators':100
        }
# instantiate the classifier
xgb_clf = XGBClassifier(**params)

# fit the classifier to the training data
xgb_clf.fit(X_train, y_train)

# Make predictions with XGBoost Classifier
# make predictions on test data
y_pred = xgb_clf.predict(X_test)

# check accuracy score
from sklearn.metrics import accuracy_score
print('XGBoost model accuracy score: {0:0.4f}'. format(accuracy_score(y_test, y_pred)))

# k-fold Cross Validation using XGBoost
from xgboost import cv

params = {"objective":"binary:logistic",'colsample_bytree': 0.3,'learning_rate': 0.1,
                'max_depth': 5, 'alpha': 10}

xgb_cv = cv(dtrain=data_dmatrix, params=params, nfold=3,
                    num_boost_round=50, early_stopping_rounds=10, metrics="auc", as_pandas=True, seed=123)

# Feature importance with XGBoost

xgb.plot_importance(xgb_clf)
plt.rcParams['figure.figsize'] = [6, 4]
plt.title(region)
plt.show()
