#!/usr/bin/env python3

# Imports
import pyfiglet
import subprocess
import shutil
import os
from Bio.PDB import *
from init import *

# BINDiS: Binding Interaction Druggability Scoring
# A tool for prioritize biological targets of the proteome of microorganisms
# based of the druggability of interaction sites

# Dictionary for colors
clrs = {
    "r": "\033[1;91m",
    "g": "\033[1;92m",
    "y": "\033[1;93m",
    "b": "\033[1;94m",
    "p": "\033[1;95m",
    "c": "\033[1;96m",
    "n": "\033[1;0m",
}

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
        best_pocket = self.max_drugscore(pockets_info)

        if extension == ".cif":
            best_pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{best_pocket}_atm.cif")
        else:
            best_pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{best_pocket}_atm.pdb")

        return best_pocket, best_pocket_path

    def copy_best_pocket(self, best_pocket_path):
        if best_pocket_path:
            shutil.copy(best_pocket_path, self.output_path)
            print(f'Arquivo {os.path.basename(best_pocket_path)} copiado para {os.getcwd()}')
    
    def cleanup(self):
        base_name = os.path.splitext(os.path.basename(self.input_file))[0]
        out_dir = os.path.join(os.path.dirname(self.input_file), f"{base_name}_out")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)

def print_banner():
    banner = pyfiglet.figlet_format("BINDiS")
    print(clrs['p'] + banner + clrs['n'])
    print(clrs['b'] + "BINDiS: Binding Interaction Druggability Scoring" + clrs['n'])
    print(clrs['n'] + "A tool to prioritize biological targets of the proteome of microorganisms" + clrs['n'])
    print(clrs['n'] + "based on the druggability of interaction sites." + clrs['n'])
    print(clrs['y'] + "=" * 80 + clrs['n'])


def main():
    print_banner()
    args = args_parse()
    verify_input(args)
    output_path = args.output if args.output else os.getcwd()
    pocket_process = PocketProcess(args.input_file, output_path)
    pocket_process.run_fpocket()
    best_pocket, best_pocket_path = pocket_process.access_pockets()
    print(f'Pocket com maior score é: {best_pocket}')
    pocket_process.copy_best_pocket(best_pocket_path)
    pocket_process.cleanup()

if __name__ == "__main__":
    main()