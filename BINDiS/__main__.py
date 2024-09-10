#!/usr/bin/env python3

# Imports
import subprocess
import shutil
import os
import pandas as pd
from Initialize import *

# BINDiS: Best Interaction Network Druggability Score
# A tool for prioritizing biological targets within the proteome 
# of microorganisms based on the druggability of interaction sites

class PocketProcess:
    def __init__(self, input_file, output_path):
        self.input_file = input_file
        self.output_path = output_path

    def run_fpocket(self):
        return_code = subprocess.call(['fpocket', '-f' + self.input_file])
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

    def access_pockets(self):
        dir_name = os.path.dirname(self.input_file)
        base_name = os.path.splitext(os.path.basename(self.input_file))[0]
        extension = os.path.splitext(self.input_file)[1]
        pockets_info = os.path.join(dir_name, f"{base_name}_out", f"{base_name}_info.txt")
        filtered_pockets = self.druggable_pockets(pockets_info)

        pockets_paths = []
        for pocket in filtered_pockets:
            pocket_number = pocket['pocket_number']
            if extension == ".cif":
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.cif")
            else:
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.pdb")
            pockets_paths.append((pocket_number, pocket_path))

        return pockets_paths

    def pockets_df(seff, drug_pockets):
        df = pd.DataFrame(drug_pockets)
        df.to_csv('pockets.csv', sep='\t', index=False)
    
        return df

    def copy_and_cleanup(self, pockets_paths):
        updated_pockets_paths = []
        for pocket_number, original_path in pockets_paths:
            if os.path.isfile(original_path):
                new_path = os.path.join(self.output_path, os.path.basename(original_path))
                shutil.copy(original_path, new_path)
                updated_pockets_paths.append({
                'pocket_number': pocket_number,
                'new_path': new_path
            })
        base_name = os.path.splitext(os.path.basename(self.input_file))[0]
        out_dir = os.path.join(os.path.dirname(self.input_file), f"{base_name}_out")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)

        return updated_pockets_paths


def main():
    check_fpocket_installed()
    args = args_parse()
    verify_input(args)
    output_path = args.output if args.output else os.getcwd()
    pocket_process = PocketProcess(args.input_file, output_path)
    pocket_process.run_fpocket()
    pockets_info = pocket_process.access_pockets()
    update_pockets_paths = pocket_process.copy_and_cleanup(pockets_info)
    df = pocket_process.pockets_df(update_pockets_paths)

if __name__ == "__main__":
    main()