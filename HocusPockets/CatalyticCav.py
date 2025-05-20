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
    def __init__(self, cavities_path, rslrd_file, fukui_file, verbose_handler):
        self.cavities_path = cavities_path
        self.verbose_handler = verbose_handler
        self.rslrd_file = rslrd_file
        self.fukui_file = fukui_file
        self.fukushima_dict = self.load_fukushima_data()

    import os

    def process_fukui_cav(self):
        """
        Identifica resíduos catalíticos com base nos maiores valores de Fukui
        extraídos da coluna B-factor do arquivo PDB. Os resíduos com pelo menos
        um átomo entre os 20% maiores valores são considerados catalíticos e são
        salvos em um novo arquivo PDB.
        """
        print("entrou")
        if not os.path.exists(self.fukui_file):
            self.verbose_handler.output(f"[ERROR] Fukui file {self.fukui_file} not found.")
            return

        atoms = []
        values = []

        # Leitura do arquivo PDB
        with open(self.fukui_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    try:
                        bfactor = float(line[60:66].strip())
                        values.append(bfactor)
                        atoms.append((line, bfactor))
                    except ValueError:
                        continue

        if not values:
            self.verbose_handler.output("[WARNING] Fukui value not found in B-factor.")
            return

        # Define o ponto de corte para os 20% maiores valores
        values.sort()
        cutoff_index = int(0.8 * len(values))
        cutoff_value = values[cutoff_index]

        self.verbose_handler.output(f"[INFO] Fukui cutoff: {cutoff_value:.4f} (top 20%)")

        # Seleciona linhas com valor acima do cutoff
        selected_lines = []
        selected_residues = set()

        for line, b in atoms:
            if b >= cutoff_value:
                res_id = line[22:26].strip() + line[17:20].strip()  # número + nome
                selected_residues.add(res_id)

        # Salva todos os átomos dos resíduos selecionados
        for line, _ in atoms:
            res_id = line[22:26].strip() + line[17:20].strip()
            if res_id in selected_residues:
                selected_lines.append(line)

        # Cria diretório e salva o arquivo
        fukui_dir = os.path.join(self.cavities_path, 'fukui')
        os.makedirs(fukui_dir, exist_ok=True)
        output_file = os.path.join(fukui_dir, 'fukui_catalytic.pdb')

        with open(output_file, 'w') as out:
            out.writelines(selected_lines)


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

                        