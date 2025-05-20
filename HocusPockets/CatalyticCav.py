#!/usr/bin/env python3

# Description
###############################################################################
'''
Descrição do que esse arquivo faz
'''

# Imports
###############################################################################
from Initialize import *


# License
###############################################################################
'''
HocusPockets: Holistic Computacional Utility for Scoring Pockets

Authors: Machado, L.F.; Pascutti, P.G.; Torres, P.H.M.

[Federal University of Rio de Janeiro]

Contact info:
Carlos Chagas Filho Institute of Biophysics
LMDM - Laboratory for Molecular Modeling and Dynamics
Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
Health Sciences Centre - block G1-19
E-mail address: luizafuruno@ufrj.br
'''

class CatalyticCav:
    def __init__(self, cavities_path, rslrd_file, verbose_handler):
        self.cavities_path = cavities_path
        self.verbose_handler = verbose_handler
        self.rslrd_file = rslrd_file
        self.fukushima_dict = self.load_fukushima_data()

    def load_fukushima_data(self):
        fukushima_data = {}
        with open(self.rslrd_file, 'r') as file:
            for line in file:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split()
                res_id = parts[0]  # Ex: "1MET"
                fukushima_value = float(parts[13])  # Fukushima está na 14ª coluna (índice 13)
                fukushima_data[res_id] = fukushima_value
        return fukushima_data
    
    def get_residues_from_pdb(self, pdb_path):
        residues = set()
        with open(pdb_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    res_seq = int(line[22:26].strip()) - 1
                    res_name = line[17:20].strip()
                    residue_id = f"{res_seq}{res_name}"
                    residues.add(residue_id)
        return residues

    def compute_fukushima_average(self):
        cavity_programs = ['fpocket', 'ligsite', 'surfnet', 'pocketfinder']
        results = {}

        for program in cavity_programs:
            program_dir = os.path.join(self.cavities_path, program)
            if not os.path.exists(program_dir):
                continue

            for pocket_file in os.listdir(program_dir):
                if pocket_file.startswith("residue") and pocket_file.endswith(".pdb"):
                    pocket_path = os.path.join(program_dir, pocket_file)
                    residues = self.get_residues_from_pdb(pocket_path)
                    
                    fuku_values = [
                        self.fukushima_dict[res] for res in residues if res in self.fukushima_dict
                    ]

                    if fuku_values:
                        avg_fukushima = sum(fuku_values) / len(fuku_values)
                    else:
                        avg_fukushima = None

                    results[(program, pocket_file)] = avg_fukushima

        return results

                        