"""
To represent the structure of nanomaterials
SDEC FP has been employed and supplemented calculation of multi-components in this code.

This code is to generate the SDEC FP for a nanomaterial with function get_component_amounts.

Created by Jaehyeon Park, source from Ph.D Shin
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import json
from formula_utils import log_transform
from amount_calculator import get_component_amounts
import re

# Load volume data
shell_volume_data = pd.read_csv('shell_volume_list.csv')
doping_volume_data = pd.read_csv('doping_volume_list.csv')
core_volume_data = pd.read_csv('core_volume_list.csv')
coating_volume_data = pd.read_csv('coating_volume_list.csv')

def calculate_electronic_configuration(amount_components, df_atom_map):
    """Calculate electronic configuration for components"""
    ec = {}
    for atom, value in amount_components.items():
        if atom in df_atom_map.index:
            row = df_atom_map.loc[atom]
            calculated_row = row * value
            ec[atom] = calculated_row
    return ec

def calculate_sdec_fp(data, df_atom):
    """Calculate SDEC fingerprint for the data"""
    # Prepare atom mapping
    df_atom_map = df_atom.set_index('atom')
    try:
        df_atom_map = df_atom_map.drop(['AN'], axis=1)
    except:
        pass
        
    sdec_fp = []
    
    for idx, row in data.iterrows():
        core = row['Core']
        doping = row['Doping']
        shell = row['Shell']
        coating = row['Coating']
        num_core = row['Amounts of Core']
        num_doping = row['Amounts of Doping']
        num_shell = row['Amounts of Shell']
        num_coating = row['Amounts of Coating']

        # CORE COMPONENT NUMBER
        amount_component_core = {}
        match_core = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', core)
        core_component = {elem: float(count) if count else 1.0 for elem, count in match_core}
        for elem, count in core_component.items():
            amount_in_core = count * num_core
            amount_component_core[elem] = float(amount_in_core) if amount_in_core else 0

        # DOPING COMPONENT NUMBER
        amount_component_doping = {}
        if doping != '':
            match_doping = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(doping))
            match_doping_count = defaultdict(float)
            for elem, count in match_doping:
                count = float(count) if count else 1
                match_doping_count[elem] += count
            doping_component = {elem: count for elem, count in match_doping_count.items()}
            if len(doping_component) >= 2:
                doping_amount_list = [float(x) for x in num_doping.split('/')]
                for (elem, count), amount_each in zip(doping_component.items(), doping_amount_list):
                    amount_in_doping = count * amount_each
                    amount_component_doping[elem] = amount_in_doping
            else:
                doping_amount_list = float(num_doping)
                for elem, count in doping_component.items():
                    amount_in_doping = count * doping_amount_list
                    amount_component_doping[elem] = amount_in_doping
            pass
        # SHELL COMPONENT NUMBER
        amount_component_shell = {}
        if shell != '':
            match_shell = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(shell))
            match_shell_count = defaultdict(float)
            for elem, count in match_shell:
                count = float(count) if count else 1
                match_shell_count[elem] += count
            shell_component = {elem: count for elem, count in match_shell_count.items()}
            if shell_component:
                for elem, count in shell_component.items():
                    amount_in_shell = count * num_shell
                    amount_component_shell[elem] = amount_in_shell
            pass

        # COATING COMPONENT NUMBER
        amount_component_coating = {}
        if coating != '':
            if coating in coating_volume_data['Coating name'].values:
                coating_mf = coating_volume_data.loc[coating_volume_data['Coating name'] == coating, 'mf'].values[0]
                match_coating = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(coating_mf))
                match_coating_count = defaultdict(float)
                for elem, count in match_coating:
                    count = float(count) if count else 1
                    match_coating_count[elem] += count
                coating_component = {elem: count for elem, count in match_coating_count.items()}
                if coating_component:
                    for elem, count in coating_component.items():
                        amount_in_coating = count * num_coating
                        amount_component_coating[elem] = amount_in_coating
            elif '/' in coating:
                coating_list = coating.split('/')
                coating_total_list = []
                for coating_sub in coating_list:
                    if coating_sub in coating_volume_data['Coating name'].values:
                        coating_mf = coating_volume_data.loc[coating_volume_data['Coating name'] == coating_sub, 'mf'].values[0]
                        coating_total_list.append(coating_mf)
                    else:
                        coating_total_list.append(coating_sub)
                match_coating = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(coating_total_list))
                match_coating_count = defaultdict(float)
                for elem, count in match_coating:
                    count = float(count) if count else 1
                    match_coating_count[elem] += count
                coating_component = {elem: count for elem, count in match_coating_count.items()}
                if coating_component:
                    for elem, count in coating_component.items():
                        amount_in_coating = count * num_coating
                        amount_component_coating[elem] = amount_in_coating                    
            pass
    
    core_ec = {}
    for atom, value in amount_component_core.items():
        if atom in df_atom_map.index:
            row = df_atom_map.loc[atom]
            calculated_row = row * value
            core_ec[atom] = calculated_row

    doping_ec = {}
    for atom, value in amount_component_doping.items():
        if atom in df_atom_map.index:
            row = df_atom_map.loc[atom]
            caculated_row = row * value
            doping_ec[atom] = calculated_row

    shell_ec = {}
    for atom, value in amount_component_shell.items():
        if atom in df_atom_map.index:
            row = df_atom_map.loc[atom]
            calculated_row = row * value
            shell_ec[atom] = calculated_row

    coating_ec = {}
    for atom, value in amount_component_coating.items():
        if atom in df_atom_map.index:
            row = df_atom_map.loc[atom]
            calculated_row = row * value
            coating_ec[atom] = calculated_row

    df_atom_map_config = df_atom_map.transpose()
    combined_ec = pd.Series(0, index=df_atom_map_config.index)
    components = {**core_ec, **doping_ec, **coating_ec, **shell_ec}

    for sub, config in components.items():
        combined_ec += config
        
    sdec_data = {}
    for ec, value in combined_ec.items():
        sdec_data[ec] = value
    sdec_fp.append(sdec_data)
    
    return pd.DataFrame(sdec_fp)

