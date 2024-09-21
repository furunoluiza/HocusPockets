#!/usr/bin/env python3

# Imports
import subprocess
import shutil
import os
import pandas as pd
from Initialize import *

# BINDiS: A Tool for Druggability-Driven and Reactivity-Based Target Prioritization

class CavityFinder:
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.output_path = output_path

class FpocketFinder(CavityFinder):
    def __init__(self, input_path, output_path):
       super().__init__(input_path, output_path)

    def run_fpocket(self, input_file):
        return_code = subprocess.call(['fpocket', '-f' + input_file])
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
    
    def druggable_pockets(self, pockets_info):
        pockets = self.parse_info(pockets_info)
        drug_pockets = [pocket for pocket in pockets if float(pocket['Druggability Score']) > 0.8]
        
        return drug_pockets

    def access_pockets(self, input_file):
        dir_name = os.path.dirname(input_file)
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        extension = os.path.splitext(input_file)[1]
        pockets_info = os.path.join(dir_name, f"{base_name}_out", f"{base_name}_info.txt")
        filtered_pockets = self.druggable_pockets(pockets_info)

        pockets_paths = []
        for pocket in filtered_pockets:
            pocket_number = pocket['pocket_number']
            if extension == ".cif":
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.cif")
            else:
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.pdb")

        pocket['new_path'] = pocket_path
        pockets_paths.append(pocket)

        return pockets_paths

    def pockets_df(self, drug_pockets, protein_name):
        df = pd.DataFrame(drug_pockets)
        df['protein'] = protein_name

        return df

    def copy_and_cleanup(self, pockets_paths, input_file):
        print("Pocket_path", pockets_paths)
        updated_pockets_paths = []
        for pocket in pockets_paths:
            pocket_number = pocket['pocket_number']
            original_path = pocket['new_path']
            if os.path.isfile(original_path):
                new_path = os.path.join(self.output_path, os.path.basename(original_path))
                shutil.copy(original_path, new_path)
                updated_pockets_paths.append({
                'pocket_number': pocket_number,
                'new_path': new_path
            })
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        out_dir = os.path.join(os.path.dirname(input_file), f"{base_name}_out")
        print("Output: " , out_dir)
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)

        return updated_pockets_paths
    
    def process_directory(self):
        # If the input is a directory, process all .pdb files in it
        all_dfs = []
        for root, _, files in os.walk(self.input_path):
            for file in files:
                if file.endswith(('.pdb', '.cif')):
                    input_file = os.path.join(root, file)
                    self.run_fpocket(input_file)
                    filtered_pockets = self.access_pockets(input_file)
                    updated_pockets_paths = self.copy_and_cleanup(filtered_pockets, input_file)
                    pockets_df = self.pockets_df(updated_pockets_paths, os.path.basename(input_file))
                    all_dfs.append(pockets_df)
                    for pocket in updated_pockets_paths:
                        pocket_number = pocket['pocket_number']
                        new_path = pocket['new_path']
                        pockets_df.loc[pockets_df['pocket_number'] == pocket_number, 'pocket_path'] = new_path

        # Concatenate all dataframes into one
        if all_dfs:
            final_df = pd.concat(all_dfs, ignore_index=True)
            final_df.to_csv(os.path.join(self.output_path, 'all_pockets.csv'), sep='\t', index=False)


def main():
    check_fpocket_installed()
    args = args_parse()
    verify_input(args)
    output_path = args.output if args.output else os.getcwd()
    cavity_finder = FpocketFinder(args.input_path, output_path)
    if os.path.isdir(args.input_path):
        cavity_finder.process_directory()
    else:
        cavity_finder.run_fpocket(args.input_path)
        pockets_info = cavity_finder.access_pockets(args.input_path)
        update_pockets_paths = cavity_finder.copy_and_cleanup(pockets_info, args.input_path)
        df = cavity_finder.pockets_df(update_pockets_paths, os.path.basename(args.input_path))
        df.to_csv('pockets.csv', sep='\t', index=False)

if __name__ == "__main__":
    main()