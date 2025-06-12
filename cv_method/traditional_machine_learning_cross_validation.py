"""
The modified Cross Validation for classical machine learning model

must_x_fold: the data is included mandatorily in training dataset

This cross validation worked with only training datasets

Created by Lukas Park (Jaehyeon Park)
"""

import clone
from sklearn.model_selection import KFold
from sklearn.metrics import (mean_absolute_error,
                             mean_absolute_percentage_error,
                             mean_squared_error,
                             root_mean_squared_error,
                             r2_score)
import pandas as pd
import numpy as np



def modified_cross_validation(model, x_train, y_train, must_x_fold, must_y_fold, fold=3):
        cv_results = {
            'Fold': [],
            'CV_MSE': [],
            'CV_RMSE': [],
            'CV_R2': [],
            'CV_MAPE': [],
            'CV_MAE':[]
        }
        
        final_results = {
            'CV_MSE_mean': None,
            'CV_MSE_std': None,
            'CV_RMSE_mean': None,
            'CV_RMSE_std': None,
            'CV_R2_mean': None,
            'CV_R2_std': None,
            'CV_MAPE_mean': None,
            'CV_MAPE_std': None,
            'CV_MAE_mean': None,
            'CV_MAE_std': None
        }
        
        kf = KFold(n_splits=fold, shuffle=True, random_state=123)
        
        fold_mse = []
        fold_rmse = []
        fold_r2 = []
        fold_mape = []
        fold_mae = []
        
        for fold_idx, (train_idx, val_idx) in enumerate(kf.split(x_train)):
            model_clone = clone(model)
            
            x_fold_train = x_train.iloc[train_idx]
            y_fold_train = y_train.iloc[train_idx]
            x_fold_val = x_train.iloc[val_idx]
            y_fold_val = y_train.iloc[val_idx]
            
            x_fold_train = pd.concat([x_fold_train, must_x_fold])
            y_fold_train = pd.concat([y_fold_train, must_y_fold])
            
            model_clone.fit(x_fold_train, y_fold_train)
            
            val_pred = model_clone.predict(x_fold_val)
            
            mse = mean_squared_error(y_fold_val, val_pred)
            rmse = np.sqrt(mse)
            r2 = r2_score(y_fold_val, val_pred)
            mape = mean_absolute_percentage_error(y_fold_val, val_pred)
            mae = mean_absolute_error(y_fold_val, val_pred)
            
            fold_mse.append(mse)
            fold_rmse.append(rmse)
            fold_r2.append(r2)
            fold_mape.append(mape)
            fold_mae.append(mae)
            
            cv_results['Fold'].append(fold_idx + 1)
            cv_results['CV_MSE'].append(mse)
            cv_results['CV_RMSE'].append(rmse)
            cv_results['CV_R2'].append(r2)
            cv_results['CV_MAPE'].append(mape)
            cv_results['CV_MAE'].append(mae)
        
        # Calculate final results
        final_results['CV_MSE_mean'] = np.mean(fold_mse)
        final_results['CV_MSE_std'] = np.std(fold_mse)
        final_results['CV_RMSE_mean'] = np.mean(fold_rmse)
        final_results['CV_RMSE_std'] = np.std(fold_rmse)
        final_results['CV_R2_mean'] = np.mean(fold_r2)
        final_results['CV_R2_std'] = np.std(fold_r2)
        final_results['CV_MAPE_mean'] = np.mean(fold_mape)
        final_results['CV_MAPE_std'] = np.std(fold_mape)
        final_results['CV_MAE_mean'] = np.mean(fold_mae)
        final_results['CV_MAE_std'] = np.std(fold_mae)
        
        return cv_results, final_results