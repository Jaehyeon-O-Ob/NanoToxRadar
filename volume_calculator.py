"""
Volume Calculator

To calculate the volume of each components in nano particle
Core / Shell / Coating / Doping
"""
import warnings
warnings.filterwarnings('ignore')

import sys
import numpy as np
import re
from collections import Counter, defaultdict
import pandas as pd
import itertools
from itertools import product
import statistics
from radii_collection import metallic_radii, effective_ionic_radii, neutral_radii
from formula_utils import (get_possible_charges, 
                           find_valid_combinations,
                           calculate_stability,
                           calculate_stability_multiple,
                           calculate_stability_single,
                           parse_molecular_formula)
from rdkit import Chem

"""
Initialize constants and base data
"""
pt = Chem.GetPeriodicTable()
all_symbols = [pt.GetElementSymbol(i) for i in range(1, 119)]
valid_elements_regex = '|'.join(sorted(all_symbols, key=len, reverse=True))

# Load volume data
shell_volume_data = pd.read_csv('shell_volume_list.csv')
doping_volume_data = pd.read_csv('doping_volume_list.csv')
core_volume_data = pd.read_csv('core_volume_list.csv')
coating_volume_data = pd.read_csv('coating_volume_list.csv')

class FormulaError(Exception):
    """Custom exception for formula validation errors"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

"""
Base utility functions
"""
def formula_error_check(formula):
    """Validate chemical formula."""
    pattern = rf'({valid_elements_regex})(\d*\.?\d*)'
    elements = re.findall(pattern, formula)
    # print(f"Found elements: {elements}") - Debugging test
    
    if not elements:
        # WARNING MESSAGE
        raise FormulaError(f"No valid elements found in the formula: {formula}")
        # WARNING MESSAGE
    
    parsed_formula = ''.join(elem + count for elem, count in elements)
    # print(f"Parsed formula: {parsed_formula}") - Debugging test
    
    if parsed_formula != formula:
        # WARNING MESSAGE
        raise FormulaError(f"Invalid element symbol found in the formula: {formula}")
        # WARNING MESSAGE
        
    return formula

def sphere_volume(r):
    """Calculate sphere volume."""
    return (4/3)*np.pi*(r**3)

def sphere_surface(r):
    """Calculate sphere surface area."""
    return 4*np.pi*(r**2)

def parse_molecular_formula(formula):
    """Parse molecular formula into composition dictionary."""
    match = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(formula))
    match_count = defaultdict(float)
    for elem, count in match:
        count = float(count) if count else 1
        match_count[elem] += count
    return {elem: float(count) if count else 1.0 for elem, count in match_count.items()}

"""
Charge calculation functions
"""
def calculate_stability_multiple(combo, metals, total_required_charge):
    """Calculate stability score for multiple metal combinations."""
    stability_score = 0
    for (metal, count), charges in zip(metals.items(), combo):
        ideal_charge_per_atom = total_required_charge / sum(metals.values())
        metal_score = count * sum(abs(charge - ideal_charge_per_atom) for charge in charges)
        stability_score += metal_score
    return stability_score

def calculate_stability_single(combo, metals, total_required_charge):
    """Calculate stability score for single metal."""
    stability_score = 0
    for (metal, count), charge in zip(metals.items(), combo):
        ideal_charge_per_atom = total_required_charge / sum(metals.values())
        metal_score = count * abs(charge - ideal_charge_per_atom)
        stability_score += metal_score
    return stability_score

def calculate_stability(combo, metals, total_required_charge):
    """Calculate stability score based on charge combination type."""
    
    if not combo:
        raise FormulaError("Empty charge combination")
    try:
        if isinstance(combo[0], (list, tuple)):
            return calculate_stability_multiple(combo, metals, total_required_charge)
        else:
            return calculate_stability_single(combo, metals, total_required_charge)
    except IndexError:
        raise FormulaError("Invalid charge combination structure")
    
def mc_np_vol_surface(data):
    for idx, row in data.iterrows():
        if row['Diameter(nm)']:
            data.loc[idx, 'Particle Volume (nm^3)'] = sphere_volume(row['Diameter(nm)'])
            data.loc[idx, 'Particle Surface Area (nm^2)'] = sphere_surface(row['Diameter(nm)'])
            
    return data

"""
CORE VOLUME (STRING)
"""

def core_volume_process(data):
    for idx, row in data.iterrows():
        core = row['Core']
        if core in core_volume_data['Core'].values:
            core_volume = core_volume_data.loc[core_volume_data['Core'] == core, 'Core Volume (nm^3)'].values[0]
            data.loc[idx, 'Core Volume (nm^3)'] = core_volume
        else:
            try:
                formula_check = formula_error_check(core)
                
                elements = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', core)
                composition = {elem: float(count) if count else 1.0 for elem, count in elements}
                stable_core = []
                # IF CORE HAS OXYGEN - EFFECTIVE RADII
                if 'O' in composition:
                    oxygen_count = composition['O']
                    oxygen_charge = -2*oxygen_count
                    has_float = any(isinstance(count, float) and not count.is_integer() for count in composition.values())
                    # IF CORE HAS FLOAT NUMBER
                    if has_float:
                        metals = {elem: float(count) if count else 1.0 for elem, count in elements if elem != 'O'}
                        possible_charges = {}
                        for metal in metals.keys():
                            possible_charges[metal] = []
                            for charge_key in effective_ionic_radii.keys():
                                match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                if match:
                                    charge_element, charge = match.groups()
                                    if charge_element == metal:
                                        possible_charges[metal].append(int(charge))
                                        
                        valid_combinations = []
                        exact_combinations = []
                        approx_combinations = []
                        for combo in product(*(possible_charges[elem] for elem in metals)):
                            total = sum(metals[elem] * state for elem, state in zip(metals, combo))
                            if round(total + oxygen_charge, 2) == 0:
                                exact_combinations.append(combo)
                            elif abs(round(total + oxygen_charge, 2)) <= 2:
                                approx_combinations.append(combo)
                        if exact_combinations:
                            valid_combinations = exact_combinations
                        elif approx_combinations:
                            valid_combinations = approx_combinations
                
                        if valid_combinations:
                            most_stable = min(valid_combinations, key=lambda x: (
                                abs(round(sum(metals[elem] * charge for elem, charge in zip(metals, x)) + oxygen_charge, 2)),
                                statistics.stdev(x) if len(set(x)) > 1 else 0))
                            for subs, charge in zip(list(metals.keys()), most_stable):
                                stable_core.append(f"{subs}+{charge}")
                        else:
                            ## WARNING MESSAGE
                            print(f"Core material error: {core} is incorrect.")
                            sys.exit(1)
                            ## WARNING MESSAGE
                            
                        stable_core.extend(['O-2'] * int(oxygen_count))
                    # CORE HAS NO FLOAT
                    else:
                        metals = [elem for elem in composition if elem != 'O']
                        metals_comp = {elem: int(count) if count else 1 for elem, count in elements if elem !='O'}
                        if len(metals) == 1:
                            element = metals[0]
                            possible_charges = []
                            for charge_key in effective_ionic_radii.keys():
                                match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                if match:
                                    charge_element, charge = match.groups()
                                    if charge_element == element:
                                        possible_charges.append(int(charge))
                                        
                                        
                            charge_combinations = itertools.combinations_with_replacement(possible_charges, int(composition[element]))
                
                            valid_combinations = []
                            approx_combinations = []
                            exact_combinations = []
                            for combo in charge_combinations:
                                total_charge = sum(combo)
                                if total_charge + oxygen_charge == 0:
                                    exact_combinations.append(combo)
                                elif abs(total_charge + oxygen_charge) <= 2:
                                    approx_combinations.append(combo)
                                    
                            if exact_combinations:
                                valid_combinations = exact_combinations
                            elif approx_combinations:
                                valid_combinations = approx_combinations
                            else:
                                ## WARNING MEESAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MEESAGE
                            if valid_combinations:
                                most_stable = min(valid_combinations, 
                                                key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                for charge in most_stable:
                                    stable_core.append(f"{element}+{charge}")
                            else:
                                ## WARNING MEESAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MESSAGE
                            stable_core.extend(['O-2'] * int(oxygen_count))
                        # CORE HAS UPPER 2 METALS
                        elif len(metals) >= 2:
                            possible_charges = {}
                            for metal in metals:
                                possible_charges[metal] = []
                                for charge_key in effective_ionic_radii.keys():
                                    match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                    if match:
                                        charge_element, charge = match.groups()
                                        if charge_element == metal:
                                            possible_charges[metal].append(int(charge))

                
                            charge_combinations = itertools.product(*(itertools.combinations_with_replacement(possible_charges[metal], int(composition[metal])) for metal in metals))
                            
                            valid_combinations = []
                            approx_combinations = []
                            exact_combinations = []
                            for combo in charge_combinations:
                                total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                if total_charge + oxygen_charge == 0:
                                    exact_combinations.append(combo)
                                elif abs(total_charge + oxygen_charge) <= 2:
                                    approx_combinations.append(combo)
                            if exact_combinations:
                                valid_combinations = exact_combinations
                            elif approx_combinations:
                                valid_combinations = approx_combinations
                            else:
                                ## WARNING MESSAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MESSAGE
                            if valid_combinations:
                                most_stable = min(valid_combinations, 
                                                key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                for metal, charges in zip(metals_comp.keys(), most_stable):
                                    for charge in charges:
                                        stable_core.append(f"{metal}+{charge}")
                            else:
                                ## WARNING MESSAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MESSAGE
                            stable_core.extend(['O-2'] * int(oxygen_count))                    
            
            
                # IF CORE HAS NO OXYGEN
                else:
                    has_float = any(isinstance(count, float) and not count.is_integer() for count in composition.values())
                    if has_float:
                        metals = {elem: float(count) if count else 1.0 for elem, count in elements if elem != 'O'}
                        print(metals)
                        possible_charges = {}
                        for metal in metals.keys():
                            possible_charges[metal] = []
                            for charge_key in effective_ionic_radii.keys():
                                match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                if match:
                                    charge_element, charge = match.groups()
                                    if charge_element == metal:
                                        possible_charges[metal].append(int(charge))
                                        
                                        
                        valid_combinations = []
                        exact_combinations = []
                        approx_combinations = []
                        for combo in product(*(possible_charges[elem] for elem in metals)):
                            total = sum(metals[elem] * state for elem, state in zip(metals, combo))
                            if total == 0:
                                exact_combinations.append(combo)
                            elif abs(total) <= 2:
                                approx_combinations.append(combo)
                        if exact_combinations:
                            valid_combinations = exact_combinations
                        elif approx_combinations:
                            valid_combinations = approx_combinations
                
                        if valid_combinations:
                            most_stable = min(valid_combinations, key=lambda x: (
                                abs(round(sum(metals[elem] * charge for elem, charge in zip(metals, x)), 2)),
                                statistics.stdev(x) if len(set(x)) > 1 else 0))
                            for subs, charge in zip(list(metals.keys()), most_stable):
                                if charge > 0:
                                    stable_core.append(f"{subs}+{charge}")
                                else:
                                    stable_core.append(f"{subs}{charge}")
                        else:
                            ## WARNING MESSAGE
                            print(f"Core material error: {core} is incorrect.")
                            sys.exit(1)
                            ## WARNING MESSAGE
                    else:
                        metals = [elem for elem in composition]
                        metals_comp = {elem: int(count) if count else 1 for elem, count in elements if elem !='O'}
                        if len(metals) == 1:
                            element = metals[0]
                            stable_core.append(element)
                        
                        else:
                            possible_charges = {}
                            for metal in metals:
                                possible_charges[metal] = []
                                for charge_key in effective_ionic_radii.keys():
                                    match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                    if match:
                                        charge_element, charge = match.groups()
                                        if charge_element == metal:
                                            possible_charges[metal].append(int(charge))

                
                            charge_combinations = itertools.product(*(itertools.combinations_with_replacement(possible_charges[metal], int(composition[metal])) for metal in metals))
                            
                            valid_combinations = []
                            approx_combinations = []
                            exact_combinations = []
                            for combo in charge_combinations:
                                total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                if total_charge == 0:
                                    exact_combinations.append(combo)
                                elif abs(total_charge) <= 2:
                                    approx_combinations.append(combo)
                            if exact_combinations:
                                valid_combinations = exact_combinations
                            elif approx_combinations:
                                valid_combinations = approx_combinations
                            else:
                                ## WARNING MESSAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MESSAGE
                            if valid_combinations:
                                most_stable = min(valid_combinations, key=lambda x: 
                                                (abs(round(sum(metals_comp[elem] * charge[0] for elem, charge in zip(metals_comp.keys(), x)), 2)),
                                                statistics.stdev([charge[0] for charge in x]) if len(set(charge[0] for charge in x)) > 1 else 0))
                                for subs, charge in zip(metals_comp.keys(), most_stable):
                                    if charge > 0:
                                        stable_core.append(f"{subs}+{charge[0]}")
                                    else:
                                        stable_core.append(f"{subs}{charge[0]}")
                            else:
                                ## WARNING MESSAGE
                                print(f"Core material error: {core} is incorrect.")
                                sys.exit(1)
                                ## WARNING MESSAGE
        
                core_data = dict(Counter(stable_core))
                # VOLUME CALCULATOR
                if len(core_data) == 1:
                    for subs, count in core_data.items():
                        radius = metallic_radii.get(subs)
                        if radius:
                            data.loc[idx, 'Core Volume (nm^3)'] = float(count*sphere_volume(radius/1000))
                        else:
                            ## WARNING MESSAGE
                            print(f"Sorry, {subs} is out of domain")
                            sys.exit(1)
                            ## WARNING MESSAGE
                else:
                    total_volume = 0
                    for subs, count in core_data.items():
                        radius = effective_ionic_radii.get(subs)
                        if radius:
                            volume = sphere_volume(radius/1000)
                            has_float = any(isinstance(count, float) and not count.is_integer() for count in composition.values())
                            if has_float:
                                element = subs.split('+')[0].split('-')[0]
                                element = element.strip()
                                ratio = composition.get(element, 1.0)
                                total_volume += volume * ratio
                            else:
                                total_volume += volume * count
                        else:
                            ## WARNING MESSAGE
                            print(f"Sorry, {subs} is out of domain")
                            sys.exit(1)
                            ## WARNING MESSAGE
                    data.loc[idx, 'Core Volume (nm^3)'] = float(total_volume)
                    
            except FormulaError as e:
                print(f"Core material error: {core} is incorrect. {str(e)}")
                sys.exit(1)
                
    return data
                
"""
DOPING VOLUME (STRING)
"""
def doping_volume_process(data):
    for idx, row in data.iterrows():
        doping = row['Doping']
        # IF DOPING ITEMS ARE MULTIPLE
        if '/' in doping:
            doping_list = doping.split('/')
            doping_volume = []
            for elem in doping_list:
                volume = doping_volume_data.loc[doping_volume_data['Doping'] == elem, 'Doping Volume (nm^3)'].values[0]
                doping_volume.append(f"{volume}")
            data.loc[idx, 'Doping Volume (nm^3)'] = '/'.join(f"{float(x)}" for x in doping_volume)
        
        # IF DOPING ITEMS ARE SINGLE
        else:
            if doping in doping_volume_data['Doping'].values:
                doping_volume = doping_volume_data.loc[doping_volume_data['Doping'] == doping, 'Doping Volume (nm^3)'].values[0]
                data.loc[idx, 'Doping Volume (nm^3)'] = str(doping_volume)
                
            # IF DOPING ITEMS ARE NOTHING
            elif doping == '':
                data.loc[idx, 'Doping Volume (nm^3)'] = str(0)
                
    return data

"""
SHELL VOLUME (NOT STRING)
"""

def shell_volume_process(data):
    for idx, row in data.iterrows():
        shell = row['Shell']
        # IF SHELL ITEMS ARE MULTIPLE
        if '/' in shell:
            shell_list = shell.split('/')
            shell_volume = 0
            for elem in shell_list:
                volume = shell_volume_data.loc[shell_volume_data['Shell'] == elem, 'Shell Volume (nm^3)'].values[0]
                shell_volume += volume
            data.loc[idx, 'Shell Volume (nm^3)'] = shell_volume
        
        # IF SHELL ITEMS ARE SINGLE
        else:
            if shell in shell_volume_data['Shell'].values:
                shell_volume = shell_volume_data.loc[shell_volume_data['Shell'] == shell, 'Shell Volume (nm^3)'].values[0]
                data.loc[idx, 'Shell Volume (nm^3)'] = shell_volume
                
            # IF SHELL ITEMS ARE NOTHING
            elif shell == '':
                data.loc[idx, 'Shell Volume (nm^3)'] = 0
                
    return data

"""
COATING VOLUME - USERS INPUT THE MOLECULAR FORMULA DIRECTLY (NOT STRING)
(pm) -> (nm)
"""

def coating_volume_process(data):
    for idx, row in data.iterrows():
        coating = row['Coating']
        # IF COATING ITEMS ARE SINGLE
        if '/' not in coating:
            if coating in coating_volume_data['Coating name'].values:
                coating_volume = coating_volume_data.loc[coating_volume_data['Coating name'] == coating, 'Coating Volume (nm^3)'].values[0]
                data.loc[idx, 'Coating Volume (nm^3)'] = coating_volume
                
            # IF COATING ITEMS ARE NOTHING
            elif coating == '':
                data.loc[idx, 'Coating Volume (nm^3)'] = 0
                
            elif coating not in coating_volume_data['Coating name'].values:
                try:
                    formula_check = formula_error_check(coating)
                    match = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(coating))
                    match_count = defaultdict(float)
                    for elem, count in match:
                            count = float(count) if count else 1
                            match_count[elem] += count
                    coating_count = {elem: float(count) if count else 1.0 for elem, count in match_count.items()}
                    # COATING HAS NO OXYGEN
                    if 'O' not in coating_count:
                        if len(coating_count) == 1:
                            total_volume = 0
                            for subs, count in coating_count.items():
                                radius = metallic_radii.get(subs)
                                if radius:
                                    subs_volume = sphere_volume(radius / 1000)
                                    total_volume += count * subs_volume
                                else:
                                    # WARNING MESSAGE
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    # WARNING MESSAGE
                            data.loc[idx, 'Coating Volume (nm^3)'] = float(total_volume)
                        else:
                            # IF 'C' IN COATING - NEUTRAL RADII
                            if 'C' in coating_count:
                                total_volume = 0
                                for subs, count in coating_count.items():
                                    radius = neutral_radii.get(subs)
                                    if radius:
                                        subs_volume = sphere_volume(radius / 1000)
                                        total_volume += count * subs_volume
                                data.loc[idx, 'Coating Volume (nm^3)'] = float(total_volume)
                            # IF 'C' NOT IN COATING - EFFECTIVE RADII
                            else:
                                stable_coating = []
                                metals = [elem for elem in coating_count]
                                metals_comp = coating_count
                                possible_charges = {}
                                for metal in metals:
                                    possible_charges[metal] = []
                                    for charge_key in effective_ionic_radii.keys():
                                        match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                        if match:
                                            charge_element, charge = match.groups()
                                            if charge_element == metal:
                                                possible_charges[metal].append(int(charge))

                                charge_combinations = itertools.product(*(
                                    itertools.combinations_with_replacement(possible_charges[metal], int(coating_count[metal])) for metal in metals))
                                
                                valid_combinations = []
                                approx_combinations = []
                                exact_combinations = []
                                for combo in charge_combinations:
                                    total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                    if total_charge == 0:
                                        exact_combinations.append(combo)
                                    elif abs(total_charge) <= 2:
                                        approx_combinations.append(combo)
                                if exact_combinations:
                                    valid_combinations = exact_combinations
                                elif approx_combinations:
                                    valid_combinations = approx_combinations
                                else:
                                    ## WARNING MESSAGE ##
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    ## WARNING MESSAGE ##
                                if valid_combinations:
                                    most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                    for subs, charges in zip(metals, most_stable):
                                        for charge in charges:
                                            if charge > 0:
                                                stable_coating.append(f"{subs}+{charge}")
                                            else:
                                                stable_coating.append(f"{subs}{charge}")
                                else:
                                    ## WARNING MESSAGE ##
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    ## WARNING MESSAGE ##
                                coating_data = dict(Counter(stable_coating))
                                if len(coating_data) >= 2:
                                    total_volume = 0
                                    for subs, count in coating_data.items():
                                        radius = effective_ionic_radii.get(subs)
                                        if radius:
                                            subs_volume = sphere_volume(radius / 1000)
                                            total_volume += count * subs_volume
                                        else:
                                            ## WARNING MESSAGE
                                            print(f"Sorry, {subs} is out of domain")
                                            sys.exit(1)
                                            ## WARNING MESSAGE
                                    data.loc[idx, 'Coating Volume (nm^3)'] = float(total_volume)
                                        
                    # COATING HAS OXYGEN
                    else:
                        # "C" IN COATING WITH OXYGEN
                        if 'C' in coating_count:
                            total_volume = 0
                            for subs, count in coating_count.items():
                                radius = neutral_radii.get(subs)
                                if radius:
                                    subs_volume = sphere_volume(radius / 1000)
                                    total_volume += count * subs_volume
                                else:
                                    ## WARNING MESSAGE
                                    print(f"Sorry, {subs} is out of domain")
                                    sys.exit(1)
                                    ## WARNING MESSAGE
                            data.loc[idx, 'Coating Volume (nm^3)'] = float(total_volume)
                        # "C" NOT IN COATING WITH OXYGEN -> METAL + OXYGEN
                        else:
                            metals = [elem for elem in coating_count if elem != 'O']
                            metals_comp = {elem: int(count) if count else 1 for elem, count in match_count.items() if elem !='O'}
                            oxygen_count = coating_count['O']
                            oxygen_charge = -2*oxygen_count
                            stable_coating = []
                            if len(metals) == 1:
                                element = metals[0]
                                possible_charges = []
                                for charge_key in effective_ionic_radii.keys():
                                    match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                    if match:
                                        charge_element, charge = match.groups()
                                        if charge_element == element:
                                            possible_charges.append(int(charge))


                                charge_combinations = itertools.combinations_with_replacement(possible_charges, int(coating_count[element]))
                
                                valid_combinations = []
                                approx_combinations = []
                                exact_combinations = []
                                for combo in charge_combinations:
                                    total_charge = sum(combo)
                                    if total_charge + oxygen_charge == 0:
                                        exact_combinations.append(combo)
                                    elif abs(total_charge + oxygen_charge) <= 2:
                                        approx_combinations.append(combo)
                                        
                                if exact_combinations:
                                    valid_combinations = exact_combinations
                                elif approx_combinations:
                                    valid_combinations = approx_combinations
                                else:
                                    ## WARNING MESSAGE ##
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    ## WARNING MESSAGE ##
                                if valid_combinations:
                                    most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                    for charge in most_stable:
                                        stable_coating.append(f"{element}+{charge}")
                                else:
                                    ## WARNING MEESAGE
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    ## WARNING MESSAGE
                                stable_coating.extend(['O-2'] * int(oxygen_count))
                            else:
                                possible_charges = {}
                                for metal in metals:
                                    possible_charges[metal] = []
                                    for charge_key in effective_ionic_radii.keys():
                                        match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                        if match:
                                            charge_element, charge = match.groups()
                                            if charge_element == metal:
                                                possible_charges[metal].append(int(charge))

                                                
                                charge_combinations = itertools.product(*(
                                    itertools.combinations_with_replacement(possible_charges[metal], int(coating_count[metal])) for metal in metals))
                            
                                valid_combinations = []
                                approx_combinations = []
                                exact_combinations = []
                                for combo in charge_combinations:
                                    total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                    if total_charge + oxygen_charge == 0:
                                        exact_combinations.append(combo)
                                    elif abs(total_charge + oxygen_charge) <= 2:
                                        approx_combinations.append(combo)
                                if exact_combinations:
                                    valid_combinations = exact_combinations
                                elif approx_combinations:
                                    valid_combinations = approx_combinations
                                else:
                                    ## WARNING MESSAGE
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    ## WARNING MESSAGE
                                if valid_combinations:
                                    most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                    for subs, charges in zip(metals, most_stable):
                                        for charge in charges:
                                            if charge > 0:
                                                stable_coating.append(f"{subs}+{charge}")
                                            else:
                                                stable_coating.append(f"{subs}{charge}")
                                else:
                                    ## WARNING MESSAGE
                                    print(f"Coating material error: {coating} is incorrect. Give More Specific Molecular formula.")
                                    sys.exit(1)
                                    continue
                                    ## WARNING MESSAGE
                                stable_coating.extend(['O-2'] * int(oxygen_count))
                            coating_data = dict(Counter(stable_coating))
                            # VOLUME CALCULATOR
                            if coating_data:
                                total_volume = 0
                                for subs, count in coating_data.items():
                                    radius = effective_ionic_radii.get(subs)
                                    if radius:
                                        subs_volume = sphere_volume(radius / 1000)
                                        total_volume += count * subs_volume
                                    else:
                                        ## WARNING MESSAGE
                                        print(f"Sorry, {subs} is out of domain")
                                        sys.exit(1)
                                        ## WARNING MESSAGE
                                data.loc[idx, 'Coating Volume (nm^3)'] = float(total_volume)
                except FormulaError as e:
                    print(f"Coating material error: {coating} is incorrect. {str(e)}")
                    sys.exit(1)
            
        # IF COATING ITEMS ARE MULTIPLE
        else:
            coating_list = coating.split('/')
            coating_volume = []
            for coating_sub in coating_list:
                # IF COATING IN COATING VOLUME LIST
                if coating_sub in coating_volume_data['Coating name'].values:
                    subs_volume = coating_volume_data.loc[coating_volume_data['Coating name'] == coating_sub, 'Coating Volume (nm^3)'].values[0]
                    coating_volume.append(float(subs_volume))
                    
                # IF COATING NOT IN COATING VOLUME LIST
                else:
                    try:
                        formula_check = formula_error_check(coating_sub)
                        
                        match = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(coating_sub))
                        match_count = defaultdict(float)
                        for elem, count in match:
                            count = float(count) if count else 1
                            match_count[elem] += count
                        coating_count = {elem: float(count) if count else 1.0 for elem, count in match_count.items()}
                        # COATING HAS NO OXYGEN
                        if 'O' not in coating_count:
                            if len(coating_count) == 1:
                                total_volume = 0
                                for subs, count in coating_count.items():
                                    radius = metallic_radii.get(subs)
                                    if radius:
                                        subs_volume = sphere_volume(radius / 1000)
                                        total_volume += count * subs_volume
                                    else:
                                        # WARNING MESSAGE
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        # WARNING MESSAGE
                                coating_volume.append(total_volume)
                            else:
                                # IF 'C' IN COATING - NEUTRAL RADII
                                if 'C' in coating_count:
                                    total_volume = 0
                                    for subs, count in coating_count.items():
                                        radius = neutral_radii.get(subs)
                                        if radius:
                                            subs_volume = sphere_volume(radius / 1000)
                                            total_volume += count * subs_volume
                                    coating_volume.append(total_volume)
                                # IF 'C' NOT IN COATING - EFFECTIVE RADII
                                else:
                                    stable_coating = []
                                    metals = [elem for elem in coating_count]
                                    metals_comp = coating_count
                                    possible_charges = {}
                                    for metal in metals:
                                        possible_charges[metal] = []
                                        for charge_key in effective_ionic_radii.keys():
                                            match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                            if match:
                                                charge_element, charge = match.groups()
                                                if charge_element == metal:
                                                    possible_charges[metal].append(int(charge))

                                    charge_combinations = itertools.product(*(
                                        itertools.combinations_with_replacement(possible_charges[metal], int(coating_count[metal])) for metal in metals))
                                    
                                    valid_combinations = []
                                    approx_combinations = []
                                    exact_combinations = []
                                    for combo in charge_combinations:
                                        total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                        if total_charge == 0:
                                            exact_combinations.append(combo)
                                        elif abs(total_charge) <= 2:
                                            approx_combinations.append(combo)
                                    if exact_combinations:
                                        valid_combinations = exact_combinations
                                    elif approx_combinations:
                                        valid_combinations = approx_combinations
                                    else:
                                        ## WARNING MESSAGE ##
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        ## WARNING MESSAGE ##
                                    if valid_combinations:
                                        most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                        for subs, charges in zip(metals, most_stable):
                                            for charge in charges:
                                                if charge > 0:
                                                    stable_coating.append(f"{subs}+{charge}")
                                                else:
                                                    stable_coating.append(f"{subs}{charge}")
                                    else:
                                        ## WARNING MESSAGE ##
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        ## WARNING MESSAGE ##
                                    coating_data = dict(Counter(stable_coating))
                                    if len(coating_data) >= 2:
                                        total_volume = 0
                                        for subs, count in coating_data.items():
                                            radius = effective_ionic_radii.get(subs)
                                            if radius:
                                                subs_volume = sphere_volume(radius / 1000)
                                                total_volume += count * subs_volume
                                            else:
                                                ## WARNING MESSAGE
                                                print(f"Sorry, {subs} is out of domain")
                                                sys.exit(1)
                                                ## WARNING MESSAGE
                                        coating_volume.append(float(total_volume))
                                        
                        # COATING HAS OXYGEN
                        else:
                            # "C" IN COATING WITH OXYGEN
                            if 'C' in coating_count:
                                total_volume = 0
                                for subs, count in coating_count.items():
                                    radius = neutral_radii.get(subs)
                                    if radius:
                                        subs_volume = sphere_volume(radius / 1000)
                                        total_volume += count * subs_volume
                                    else:
                                        ## WARNING MESSAGE
                                        print(f"Sorry, {subs} is out of domain")
                                        sys.exit(1)
                                        ## WARNING MESSAGE
                                coating_volume.append(total_volume)
                            # "C" NOT IN COATING WITH OXYGEN -> METAL + OXYGEN
                            else:
                                metals = [elem for elem in coating_count if elem != 'O']
                                metals_comp = {elem: int(count) if count else 1 for elem, count in match_count.items() if elem !='O'}
                                oxygen_count = coating_count['O']
                                oxygen_charge = -2*oxygen_count
                                stable_coating = []
                                if len(metals) == 1:
                                    element = metals[0]
                                    possible_charges = []
                                    for charge_key in effective_ionic_radii.keys():
                                        match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                        if match:
                                            charge_element, charge = match.groups()
                                            if charge_element == element:
                                                possible_charges.append(int(charge))

            
                                    charge_combinations = itertools.combinations_with_replacement(possible_charges, int(coating_count[element]))
                    
                                    valid_combinations = []
                                    approx_combinations = []
                                    exact_combinations = []
                                    for combo in charge_combinations:
                                        total_charge = sum(combo)
                                        if total_charge + oxygen_charge == 0:
                                            exact_combinations.append(combo)
                                        elif abs(total_charge + oxygen_charge) <= 2:
                                            approx_combinations.append(combo)
                                            
                                    if exact_combinations:
                                        valid_combinations = exact_combinations
                                    elif approx_combinations:
                                        valid_combinations = approx_combinations
                                    else:
                                        ## WARNING MESSAGE ##
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        ## WARNING MESSAGE ##
                                    if valid_combinations:
                                        most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                        for charge in most_stable:
                                            stable_coating.append(f"{element}+{charge}")
                                    else:
                                        ## WARNING MEESAGE
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        ## WARNING MESSAGE
                                    stable_coating.extend(['O-2'] * int(oxygen_count))
                                else:
                                    possible_charges = {}
                                    for metal in metals:
                                        possible_charges[metal] = []
                                        for charge_key in effective_ionic_radii.keys():
                                            match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
                                            if match:
                                                charge_element, charge = match.groups()
                                                if charge_element == metal:
                                                    possible_charges[metal].append(int(charge))

                                                    
                                    charge_combinations = itertools.product(*(
                                        itertools.combinations_with_replacement(possible_charges[metal], int(coating_count[metal])) for metal in metals))
                                
                                    valid_combinations = []
                                    approx_combinations = []
                                    exact_combinations = []
                                    for combo in charge_combinations:
                                        total_charge = sum(sum(metal_combo) for metal_combo in combo)
                                        if total_charge + oxygen_charge == 0:
                                            exact_combinations.append(combo)
                                        elif abs(total_charge + oxygen_charge) <= 2:
                                            approx_combinations.append(combo)
                                    if exact_combinations:
                                        valid_combinations = exact_combinations
                                    elif approx_combinations:
                                        valid_combinations = approx_combinations
                                    else:
                                        ## WARNING MESSAGE
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        ## WARNING MESSAGE
                                    if valid_combinations:
                                        most_stable = min(valid_combinations, 
                                                        key=lambda x: calculate_stability(x, metals_comp, abs(oxygen_charge)))
                                        for subs, charges in zip(metals, most_stable):
                                            for charge in charges:
                                                if charge > 0:
                                                    stable_coating.append(f"{subs}+{charge}")
                                                else:
                                                    stable_coating.append(f"{subs}{charge}")
                                    else:
                                        ## WARNING MESSAGE
                                        print(f"Coating material error: {coating_sub} is incorrect. Give More Specific Molecular formula.")
                                        sys.exit(1)
                                        continue
                                        ## WARNING MESSAGE
                                    stable_coating.extend(['O-2'] * int(oxygen_count))
                                coating_data = dict(Counter(stable_coating))
                                # VOLUME CALCULATOR
                                if coating_data:
                                    total_volume = 0
                                    for subs, count in coating_data.items():
                                        radius = effective_ionic_radii.get(subs)
                                        if radius:
                                            subs_volume = sphere_volume(radius / 1000)
                                            total_volume += count * subs_volume
                                        else:
                                            ## WARNING MESSAGE
                                            print(f"Sorry, {subs} is out of domain")
                                            sys.exit(1)
                                            ## WARNING MESSAGE
                                    coating_volume.append(total_volume)
                    except FormulaError as e:
                        print(f"Coating material error: {coating_sub} is incorrect. {str(e)}")
                        sys.exit(1)

            total_coating_volume = sum(coating_volume)
            data.loc[idx, 'Coating Volume (nm^3)'] = total_coating_volume
            
    return data


def calculate_volumes(data):
    data1 = mc_np_vol_surface(data)
    data2 = core_volume_process(data1)
    data3 = doping_volume_process(data2)
    data4 = shell_volume_process(data3)
    data5 = coating_volume_process(data4)
    return data5