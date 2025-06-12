import pandas as pd
from collections import defaultdict
import re
from formula_utils import parse_molecular_formula

def initialize_amount_columns(data):
    """Initialize amount columns in the dataframe"""
    components = ['Core', 'Doping', 'Shell', 'Coating']
    for component in components:
        col_name = f'Amounts of {component}'
        data[col_name] = 0
        data[col_name] = data[col_name].astype(object)
    return data

def calculate_coating_amount(particle_sa, coating_vol):
    """Calculate coating amount"""
    if coating_vol != 0:
        return float(particle_sa / coating_vol)
    return 0.0

def calculate_shell_amount(particle_sa, shell_vol):
    """Calculate shell amount"""
    if shell_vol != 0:
        return float(particle_sa / shell_vol)
    return 0

def calculate_doping_amounts(particle_vol, doping_ratio, doping_vol):
    """Calculate doping amounts for single or multiple dopings"""
    if '/' in doping_ratio:
        # Multiple dopings
        doping_vols = [float(x) for x in doping_vol.split('/')]
        doping_ratios = [(float(x) / 100) for x in doping_ratio.split('/')]
        doping_amount = [(ratio * particle_vol) / volume 
                        for ratio, volume in zip(doping_ratios, doping_vols)]
        return '/'.join(f"{float(x)}" for x in doping_amount), sum(doping_ratios)
    elif '/' not in doping_ratio:
        # Single doping
        doping_ratios = float(doping_ratio) / 100
        doping_amount = (doping_ratios * particle_vol) / float(doping_vol)
        return str(doping_amount), doping_ratios
    else:
        return 0

def calculate_core_amount(particle_vol, core_vol, total_doping_ratio=0):
    """Calculate core amount"""
    return float(((1 - total_doping_ratio) * particle_vol) / core_vol)

def calculate_amounts(data):
    """Calculate amounts for all components"""
    for idx, row in data.iterrows():
        particle_vol = row['Particle Volume (nm^3)']
        particle_sa = row['Particle Surface Area (nm^2)']
        
        # Calculate coating amount
        data.loc[idx, 'Amounts of Coating'] = calculate_coating_amount(
            particle_sa, row['Coating Volume (nm^3)'])
        
        # Calculate shell amount
        data.loc[idx, 'Amounts of Shell'] = calculate_shell_amount(
            particle_sa, row['Shell Volume (nm^3)'])
        
        # Calculate doping and core amounts
        if row['Doping'] != '' and row['Doping Rate(%)'] != '':
            doping_amount, total_doping_ratio = calculate_doping_amounts(
                particle_vol, row['Doping Rate(%)'], row['Doping Volume (nm^3)'])
            data.loc[idx, 'Amounts of Doping'] = doping_amount
            data.loc[idx, 'Amounts of Core'] = calculate_core_amount(
                particle_vol, row['Core Volume (nm^3)'], total_doping_ratio)
        else:
            data.loc[idx, 'Amounts of Core'] = calculate_core_amount(
                particle_vol, row['Core Volume (nm^3)'])
            data.loc[idx, 'Amounts of Doping'] = 0
            
    return data

def get_component_amounts(formula, amount):
    """Calculate amounts for each element in a component"""
    amount_components = {}
    if formula and amount:
        components = parse_molecular_formula(formula)
        for elem, count in components.items():
            if isinstance(amount, str) and '/' in amount:
                # Handle multiple components
                amount_list = [float(x) for x in amount.split('/')]
                if len(components) >= 2:
                    for amt in amount_list:
                        amount_components[elem] = count * amt
                else:
                    amount_components[elem] = count * float(amount_list[0])
            else:
                # Handle single component
                amount_components[elem] = count * float(amount)
    return amount_components