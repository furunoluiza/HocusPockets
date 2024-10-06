#!/usr/bin/env python3

# Imports
import subprocess
import shutil
import os
import pandas as pd
from Initialize import *

# WhichCav: A Tool for Druggability-Driven and Reactivity-Based Target Prioritization

class CavityFinder:
    def __init__(self, input_path, output_path, verbose_handler, druggability_score):
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler
        self.druggability_score = druggability_score
    
    def find_cavities(self):
        base_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, base_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        if os.path.isdir(self.input_path):
            self.find_cavities_in_directory()
        else:
            self.find_cavities_in_file(protein_output_dir)
    
    def find_cavities_in_directory(self):
        fpocket_finder = FpocketFinder(self.input_path, self.output_path, self.verbose_handler, self.druggability_score)
        fpocket_finder.process_directory()
        concavity_finder = ConCavityFinder(self.input_path, self.output_path, self.verbose_handler, self.druggability_score)
        concavity_finder.process_directory()


    def find_cavities_in_file(self, protein_output_dir):
        fpocket_finder = FpocketFinder(self.input_path, protein_output_dir, self.verbose_handler, self.druggability_score)
        fpocket_finder.process_file()
        concavity_finder = ConCavityFinder(self.input_path, protein_output_dir, self.verbose_handler, self.druggability_score)
        concavity_finder.run_concavity(self.input_path)

class ConCavityFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, druggability_score):
        super().__init__(input_path, output_path, verbose_handler, druggability_score)
    
    def access_concavity(self, input_file,method_output_dir, method):
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        pockets_info = f"{base_name}_{method}_pocket.pdb"
        self.process_cavities(pockets_info, method_output_dir)
    
    def cleanup_files(self, input_file, method):
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        files_to_delete = [
            f"{base_name}_{method}_pocket.pdb",
            f"{base_name}_{method}_residue.pdb",
            f"{base_name}_{method}.dx",
            f"{base_name}_{method}.pml"
        ]
        current_directory = os.getcwd()
        for file_name in os.listdir(current_directory):
            if file_name.startswith(f"{base_name}_") and file_name.endswith(f"_{method}.scores"):
                files_to_delete.append(file_name)
        for file_name in files_to_delete:
            file_path = os.path.join(current_directory, file_name)
            if os.path.exists(file_path):
                os.remove(file_path)

    def run_concavity(self, input_file):
        grid_methods = ['pocketfinder', 'ligsite', 'surfnet']
        for method in grid_methods:
            method_output_dir = os.path.join(self.output_path, os.path.splitext(os.path.basename(input_file))[0], method)
            os.makedirs(method_output_dir, exist_ok=True)
            command = [
                'concavity',
                '-grid_method', method,
                '-extraction_method', 'search',
                '-res_map_method', 'blur',
                input_file,
                '-print_grid_pdb', '1',
                method
            ]
            return_code = subprocess.call(command)
            assert return_code == 0, clrs['r'] + f'Problem running ConCavity with {method} method' + clrs['n']
            self.access_concavity(input_file, method_output_dir, method)
            self.cleanup_files(input_file, method)
    
    def process_cavities(self, pockets_info, output_path):
        cavities = {}
        with open(pockets_info, 'r') as file:
            for line in file:
                if line.startswith("HETATM"):
                    cavity_index = line[22:26].strip() 
                    if cavity_index not in cavities:
                        cavities[cavity_index] = []
                    cavities[cavity_index].append(line)
        for index, coordinates in cavities.items():
            output_file = os.path.join(output_path, f"cavity_{index}.pdb")
            with open(output_file, 'w') as out_file:
                out_file.writelines(coordinates)
 
    def process_directory(self):
        for root, _, files in os.walk(self.input_path):
            for file in files:
                if file.endswith('.pdb'):  # Verifique as extensões dos arquivos
                    input_file = os.path.join(root, file)
                    # Crie um diretório para a proteína, se não existir
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)

                    # Ajuste o output_path para o diretório da proteína
                    self.run_concavity(input_file)  # Chame o método run_concavity para cada arquivo

class FpocketFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, druggability_score):
       super().__init__(input_path, output_path, verbose_handler, druggability_score)

    def run_fpocket(self, input_file):
        with open(os.devnull, 'w') as devnull:
            return_code = subprocess.call(['fpocket', '-f' + input_file], stdout=devnull, stderr=devnull)
        assert return_code == 0, clrs['r'] + f'Problem running fpocket tool' + clrs['n']

    @staticmethod
    def parse_info(pockets_info):
        with open(pockets_info, 'r') as file:
            lines = file.readlines()
        pockets = []
        current_pocket = {}
        for line in lines:
            if line.startswith('Pocket'):
                if current_pocket:
                    pockets.append(current_pocket)
                current_pocket = {'pocket_number': line.split()[1].strip(':')}
            elif ':' in line:
                key, value = line.split(':')
                current_pocket[key.strip()] = value.strip()
        if current_pocket:
            pockets.append(current_pocket)
        return pockets
    
    def druggable_pockets(self, pockets_info, druggability_score):
        pockets = self.parse_info(pockets_info)
        drug_pockets = [pocket for pocket in pockets if float(pocket['Druggability Score']) > 0.5]
        return drug_pockets

    def access_pockets(self, input_file):
        dir_name = os.path.dirname(input_file)
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        extension = os.path.splitext(input_file)[1]
        pockets_info = os.path.join(dir_name, f"{base_name}_out", f"{base_name}_info.txt")
        filtered_pockets = self.druggable_pockets(pockets_info, self.druggability_score)
        pockets_paths = []
        for pocket in filtered_pockets:
            pocket_number = pocket['pocket_number']
            if extension == ".cif":
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_vert.pqr")
            else:
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_vert.pqr")
            pocket['new_path'] = pocket_path
            pockets_paths.append(pocket)
        return pockets_paths

    def pockets_df(self, drug_pockets, protein_name):
        df = pd.DataFrame(drug_pockets)
        df.insert(0, 'protein', protein_name)
        return df

    def copy_and_cleanup(self, pockets_paths, input_file):
        updated_pockets_paths = []
        protein_name = os.path.splitext(os.path.basename(input_file))[0]
        protein_output_dir = os.path.join(self.output_path, protein_name)  # Diretório específico da proteína
        cavities_dir = os.path.join(protein_output_dir, 'fpocket')  # Novo diretório para as cavidades
        os.makedirs(cavities_dir, exist_ok=True)  # Cria o diretório de cavidades se não existir

        # Limpa o diretório de cavidades se já existir
        if os.path.exists(cavities_dir):
            shutil.rmtree(cavities_dir)  
            os.makedirs(cavities_dir)  # Cria novamente após remoção

        for pocket in pockets_paths:
            pocket_number = pocket['pocket_number']
            original_path = pocket['new_path']
            if os.path.isfile(original_path):
                new_path = os.path.join(cavities_dir, os.path.basename(original_path))  # Salva no diretório de cavidades
                shutil.copy(original_path, new_path)
                updated_pockets_paths.append({'pocket_number': pocket_number, 'path': new_path})

        # Limpa arquivos temporários
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        out_dir = os.path.join(os.path.dirname(input_file), f"{base_name}_out")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)  # Remove o diretório de saída do fpocket, se existir

        return updated_pockets_paths
    
    def process_file(self):
        protein_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, protein_name)  # NOVO: Define diretório específico da proteína
        os.makedirs(protein_output_dir, exist_ok=True)  # NOVO: Cria diretório se não existir
        
        self.run_fpocket(self.input_path)
        pockets_info = self.access_pockets(self.input_path)
        update_pockets_paths = self.copy_and_cleanup(pockets_info, self.input_path)
        df = self.pockets_df(update_pockets_paths, os.path.basename(self.input_path))
        output_file = os.path.join(self.output_path, 'pockets_fpocket.csv')
        df.to_csv(output_file, sep='\t', index=False)
    
    def process_directory(self):
        all_dfs = []
        for root, _, files in os.walk(self.input_path):
            for file in files:
                if file.endswith(('.pdb', '.cif')):
                    input_file = os.path.join(root, file)
                    # Muda o diretório de saída para o diretório específico da proteína
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)  # NOVO: Define diretório específico da proteína
                    os.makedirs(protein_output_dir, exist_ok=True)  # NOVO: Cria diretório se não existir

                    self.run_fpocket(input_file)
                    filtered_pockets = self.access_pockets(input_file)
                    updated_pockets_paths = self.copy_and_cleanup(filtered_pockets, input_file)
                    pockets_df = self.pockets_df(updated_pockets_paths, os.path.basename(input_file))

                    # Muda o local do arquivo CSV para o diretório específico da proteína
                    pockets_df.to_csv(os.path.join(protein_output_dir, 'pockets_fpocket.csv'), sep='\t', index=False)  # NOVO: Salva o CSV no diretório da proteína
                    all_dfs.append(pockets_df)