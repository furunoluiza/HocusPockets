#!/usr/bin/env python3

# Imports
import pyfiglet
import subprocess
import shutil
import os
from init import *

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

    def max_drugscore(self, pockets_info):
        pockets = self.parse_info(pockets_info)
        max_pocket = max(pockets, key=lambda x: float(x['Druggability Score']))
    
        return max_pocket['pocket_number']

    def access_pockets(self):
        dir_name = os.path.dirname(self.input_file)
        base_name = os.path.splitext(os.path.basename(self.input_file))[0]
        extension = os.path.splitext(self.input_file)[1]
        pockets_info = os.path.join(dir_name, f"{base_name}_out", f"{base_name}_info.txt")
        max_pocket = self.max_drugscore(pockets_info)

        if extension == ".cif":
            max_pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{max_pocket}_atm.cif")
        else:
            max_pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{max_pocket}_atm.pdb")

        return max_pocket, max_pocket_path

    def copy_and_cleanup(self, max_pocket_path):
        if max_pocket_path:
            shutil.copy(max_pocket_path, self.output_path)
        base_name = os.path.splitext(os.path.basename(self.input_file))[0]
        out_dir = os.path.join(os.path.dirname(self.input_file), f"{base_name}_out")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)

def print_banner():
    banner = pyfiglet.figlet_format("BINDiS")
    print(clrs['p'] + banner + clrs['n'])
    print(clrs['p'] + "BINDiS: Best Interaction Network Druggability Score" + clrs['n'])
    print(clrs['y'] + "=" * 80 + clrs['n'])

def main():
    print_banner()
    args = args_parse()
    verify_input(args)
    output_path = args.output if args.output else os.getcwd()
    pocket_process = PocketProcess(args.input_file, output_path)
    pocket_process.run_fpocket()
    max_pocket, max_pocket_path = pocket_process.access_pockets()
    print(f"Pocket '{max_pocket}' identified with the highest druggability score!")
    pocket_process.copy_and_cleanup(max_pocket_path)

if __name__ == "__main__":
    main()