"""
To predict with all module and make the results
"""
import warnings
warnings.filterwarnings('ignore')

from volume_calculator import calculate_volumes
from formula_utils import log_transform
from amount_calculator import calculate_amounts
from sdec_fp_generator import calculate_sdec_fp
from collections import defaultdict
from xgboost import XGBRegressor as xgb
from catboost import CatBoostRegressor
import pandas as pd
import numpy as np
import joblib
import json
import math
import re
import os

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import json
import ast

# nano particle ready
nanoparticle = {'Core':'CdSe',
               'Shell':'',
               'Doping': '',
               'Doping Rate(%)': '',
               'Coating': '',
               'Diameter(nm)': 500}

model_folder = 'model'

data = pd.DataFrame([nanoparticle])

# Tox Model
cell_line_best_catboost = CatBoostRegressor()
cell_line_best_catboost.load_model(os.path.join(model_folder, 'best_tox_catboost.cbm'))

# Needed Data
df_atom = pd.read_excel('degenerated_electronic_configuration_without_spin.xlsx')

# Only Cell lines
cell_type = pd.read_csv('cell_type_test_data.csv')

# All Cell info
cell_info = pd.read_csv('cell_all_info_test.csv')


# To generate volumes
volumes = calculate_volumes(data)

# To generate amounts
amounts = calculate_amounts(volumes)

# To generate SDEC FP
sdec_fp = calculate_sdec_fp(amounts, df_atom=df_atom)

df_sdec_log = sdec_fp.apply(lambda x: x.map(log_transform))

cell_id = [col for col in cell_info.columns if 'Cell-identification' in col]

re_sdec = pd.concat([df_sdec_log] * len(cell_type), ignore_index=True)
x_data = pd.concat([re_sdec, cell_type], axis = 1)


### Save prediction Result ###
prediction = cell_line_best_catboost.predict(x_data)
result = {cell: float(pred) for cell, pred in zip(cell_id, prediction)}
cell_tissue = [col for col in cell_info.columns if 'Cell-tissue' in col]

cell_id_with_tissue = []

for index, row in cell_info.iterrows():
    id_columns_with_1 = [col for col in cell_id if row[col] == 1]
    tissue_columns_with_1 = [col for col in cell_tissue if row[col] == 1]
    combined_columns = id_columns_with_1 + tissue_columns_with_1

    cell_id_with_tissue.append(combined_columns)

df = pd.DataFrame(cell_id_with_tissue, columns=['Cell-identification', 'Cell-tissue'])

### output: JSON ###
nested_result = defaultdict(dict)
for _, row in df.iterrows():
    cell_id = row['Cell-identification']
    cell_tissue = row['Cell-tissue'].replace('Cell-tissue-organ-origin_', '')
    
    if cell_id in result:
        nested_result[cell_tissue][cell_id] = result[cell_id]

final_result = dict(nested_result)

import pprint
print(f"final result has been saved successfully")
pprint.pprint(json.dumps(final_result))

### output: csv file ###
data_for_df = []
for _, row in df.iterrows():
    cell_id = row['Cell-identification']
    cell_tissue = row['Cell-tissue'].replace('Cell-tissue-organ-origin_', '')

    if cell_id in result:
        data_for_df.append({
            'Cell-tissue': cell_tissue.lower(),
            'Cell-identification': cell_id.split('_')[-1],
            'Prediction': round(result[cell_id],3)
        })

# Convert the list of dictionaries into a DataFrame
df_pred = pd.DataFrame(data_for_df)

# Save the DataFrame to a CSV file
df_pred.sort_values(by='Cell-tissue').to_csv('result_from_model.csv', index=False)
print(f"final result has been saved successfully")