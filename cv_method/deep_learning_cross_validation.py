"""
The modified Cross Validation for deep learning model

must_x_fold: the data is included mandatorily in training dataset

This cross validation worked with only training datasets

Created by Lukas Park (Jaehyeon Park)
"""

import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np

from sklearn.model_selection import KFold
from sklearn.metrics import (mean_absolute_error,
                             mean_absolute_percentage_error,
                             mean_squared_error,
                             root_mean_squared_error,
                             r2_score)

def modified_dl_cross_validation(model_class, 
                              model_params, 
                              x_train, y_train, 
                              must_x_fold,
                              must_y_fold, 
                              train_params, 
                              fold=3, 
                              device=None):
    
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
    
    
    # Convert data to PyTorch tensors
    x_train_tensor = torch.FloatTensor(x_train.astype(float).values).to(device)
    y_train_tensor = torch.FloatTensor(y_train.values).to(device)
    must_x_tensor = torch.FloatTensor(must_x_fold.astype(float).values).to(device)
    must_y_tensor = torch.FloatTensor(must_y_fold.values).to(device)
    
    kf = KFold(n_splits=fold, shuffle=True, random_state=123)
    
    fold_metrics = {
        'mse': [], 'rmse': [], 'r2': [], 'mape': [], 'mae': []
    }
    
    for fold_idx, (train_idx, val_idx) in enumerate(kf.split(x_train)):
        print(f"Processing Fold {fold_idx + 1}/{fold}")
        
        model = model_class(**model_params).to(device)
        optimizer = optim.Adam(model.parameters(), lr=train_params['learning_rate'])
        criterion = nn.MSELoss()
        
        x_fold_train = x_train_tensor[train_idx]
        y_fold_train = y_train_tensor[train_idx]
        x_fold_val = x_train_tensor[val_idx]
        y_fold_val = y_train_tensor[val_idx]
        
        x_fold_train = torch.cat([x_fold_train, must_x_tensor])
        y_fold_train = torch.cat([y_fold_train, must_y_tensor])
        
        train_dataset = torch.utils.data.TensorDataset(x_fold_train, y_fold_train)
        train_loader = torch.utils.data.DataLoader(
            train_dataset,
            batch_size=train_params['batch_size'],
            shuffle=True
        )
        
        val_dataset = torch.utils.data.TensorDataset(x_fold_val, y_fold_val)
        val_loader = torch.utils.data.DataLoader(
            val_dataset,
            batch_size=train_params['batch_size'],
            shuffle=False
        )
        
        # Training
        for epoch in range(train_params['epochs']):
            model.train()
            for batch_x, batch_y in train_loader:
                optimizer.zero_grad()
                outputs = model(batch_x)
                loss = criterion(outputs, batch_y)
                loss.backward()
                optimizer.step()
            
            if (epoch + 1) % 100 == 0:
                print(f"    Epoch {epoch + 1}/{train_params['epochs']}")
        
        # Calculate metrics
        model.eval()
        val_predictions = []
        val_targets = []

        with torch.no_grad():
            for batch_x, batch_y in val_loader:
                batch_pred = model(batch_x)
                val_predictions.append(batch_pred.cpu())
                val_targets.append(batch_y.cpu())
            
            val_pred = torch.cat(val_predictions).numpy()
            y_fold_val_np = torch.cat(val_targets).numpy()
            
            mse = mean_squared_error(y_fold_val_np, val_pred)
            rmse = np.sqrt(mse)
            r2 = r2_score(y_fold_val_np, val_pred)
            mape = mean_absolute_percentage_error(y_fold_val_np, val_pred)
            mae = mean_absolute_error(y_fold_val_np, val_pred)
            
            fold_metrics['mse'].append(mse)
            fold_metrics['rmse'].append(rmse)
            fold_metrics['r2'].append(r2)
            fold_metrics['mape'].append(mape)
            fold_metrics['mae'].append(mae)
            
            cv_results['Fold'].append(fold_idx + 1)
            cv_results['CV_MSE'].append(mse)
            cv_results['CV_RMSE'].append(rmse)
            cv_results['CV_R2'].append(r2)
            cv_results['CV_MAPE'].append(mape)
            cv_results['CV_MAE'].append(mae)
    
    # Calculate final cross-validation results
    for metric in ['MSE', 'RMSE', 'R2', 'MAPE', 'MAE']:
        metric_lower = metric.lower()
        final_results[f'CV_{metric}_mean'] = np.mean(fold_metrics[metric_lower])
        final_results[f'CV_{metric}_std'] = np.std(fold_metrics[metric_lower])
    
    cv_results_df = pd.DataFrame(cv_results)
    final_results_df = pd.DataFrame([final_results])
    
    return cv_results_df, final_results_df