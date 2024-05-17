#!/usr/bin/env python

# Imports
import subprocess
import shutil
from utils import *
from init import *

# BINDiS: Binding Site Detector
# A tool for identify and prioritize bindig sites on proteins
# Developed by Luíza Furuno Machado

#Run fpocket
def run_fpocket(input_file, output_path):
    return_code = subprocess.call(['fpocket', '-f' + input_file])
    assert return_code == 0, clrs['r'] + f'Problem running fpocket tool' + clrs['n']
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    fpocket_output_dir = f"{base_name}_out"
    if os.path.isdir(fpocket_output_dir):
        shutil.move(fpocket_output_dir, output_path)

    #output n esta funcionando

def main():
    args = args_parse()
    print('Welcome to ' + clrs['p'] + 'BINDiS!' + clrs['n'])
    #install_dependeces()
    verify_input(args)
    run_fpocket(args.input_file, args.output)

if __name__ == "__main__":
    main()