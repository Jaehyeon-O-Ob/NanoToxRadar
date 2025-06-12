"""
Formula Utils module

formula_error_check - this function will check whether formula, you input, is correct or not, 
if the formula has problem or gets wrong atom, you will have the answer which the formula is wrong.

calculate_stability_multiple - when the nano particle component with variable atoms have to be checked the electronic stability.

single, mulitple, normal will check each components

get_valid_combination - valid combination means that if stability of component is perfect there are valid combination by electron of atoms in the component,
for example, Fe3O4 has combination of (Fe+2, Fe+3, Fe+3, O-2, O-2, O-2, O-2) and (Fe+2, Fe+4, Fe+2, O-2, O-2, O-2, O-2)
but the former is more stable, so we need to check the valid combination by result of stability.

Created by Jaehyeon Park
"""
import numpy as np
import pandas as pd
import re
from itertools import product, combinations_with_replacement
import statistics
from collections import Counter, defaultdict
from radii_collection import effective_ionic_radii
from rdkit import Chem

class FormulaError(Exception):
    """Custom exception for formula validation errors"""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def initialize_periodic_table():
    """Initialize periodic table and valid elements regex"""
    pt = Chem.GetPeriodicTable()
    all_symbols = [pt.GetElementSymbol(i) for i in range(1, 119)]
    valid_elements_regex = '|'.join(sorted(all_symbols, key=len, reverse=True))
    return valid_elements_regex

def formula_error_check(formula, valid_elements_regex):
    """Validate chemical formula."""
    pattern = rf'({valid_elements_regex})(\d*\.?\d*)'
    elements = re.findall(pattern, formula)
    
    if not elements:
        raise FormulaError(f"No valid elements found in the formula: {formula}")
    
    parsed_formula = ''.join(elem + count for elem, count in elements)
    if parsed_formula != formula:
        raise FormulaError(f"Invalid element symbol found in the formula: {formula}")
        
    return formula        

def parse_molecular_formula(formula):
    match = re.findall(r'([A-Z][a-z]*)(\d*\.?\d*)', str(formula))
    match_count = defaultdict(float)
    for elem, count in match:
        count = float(count) if count else 1
        match_count[elem] += count
    return {elem: float(count) if count else 1.0 for elem, count in match_count.items()}

def calculate_stability_multiple(combo, metals, total_required_charge):
    stability_score = 0
    for (metal, count), charges in zip(metals.items(), combo):
        ideal_charge_per_atom = total_required_charge / sum(metals.values())
        metal_score = count * sum(abs(charge - ideal_charge_per_atom) for charge in charges)
        stability_score += metal_score
    return stability_score

def calculate_stability_single(combo, metals, total_required_charge):
    stability_score = 0
    for (metal, count), charge in zip(metals.items(), combo):
        ideal_charge_per_atom = total_required_charge / sum(metals.values())
        metal_score = count * abs(charge - ideal_charge_per_atom)
        stability_score += metal_score
    return stability_score

def calculate_stability(combo, metals, total_required_charge):
    if isinstance(combo[0], (list, tuple)):
        return calculate_stability_multiple(combo, metals, total_required_charge)
    else:
        return calculate_stability_single(combo, metals, total_required_charge)

def get_possible_charges(element, effective_ionic_radii):
    possible_charges = []
    for charge_key in effective_ionic_radii.keys():
        match = re.match(r'([A-Z][a-z]?)([+-]\d+)', charge_key)
        if match:
            charge_element, charge = match.groups()
            if charge_element == element:
                possible_charges.append(int(charge))
    return possible_charges

def find_valid_combinations(charge_combinations, oxygen_charge=0):
    valid_combinations = []
    approx_combinations = []
    exact_combinations = []
    
    for combo in charge_combinations:
        if isinstance(combo[0], (tuple, list)):
            total_charge = sum(sum(metal_combo) for metal_combo in combo)
        else:
            total_charge = sum(combo)
            
        if total_charge + oxygen_charge == 0:
            exact_combinations.append(combo)
        elif abs(total_charge + oxygen_charge) <= 2:
            approx_combinations.append(combo)
            
    if exact_combinations:
        valid_combinations = exact_combinations
    elif approx_combinations:
        valid_combinations = approx_combinations
        
    return valid_combinations


def log_transform(x):
    """Apply log transformation with sign preservation."""
    if x == 0:
        return 0
    sign = 1 if x > 0 else -1
    return sign*np.log10(1 + abs(x))
