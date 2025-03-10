#!/usr/bin/env python3

# Description
###############################################################################
'''
Descrição do que esse arquivo faz
'''

# Imports
###############################################################################
import os
import sys
import argparse
import textwrap as tw
import subprocess

# License
###############################################################################
'''
HocusPockets: Holistic Computacional Utility for Scoring Pockets

Authors: Machado, L.F.; Torres, P.H.M.

[Federal University of Rio de Janeiro]

Contact info:
Carlos Chagas Filho Institute of Biophysics
LMDM - Laboratory for Molecular Modeling and Dynamics
Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
Health Sciences Centre - block G1-19
E-mail address: luizafuruno@ufrj.br
'''

# Dictionary for colors
clrs = {
    "r": "\033[1;91m", # red
    "g": "\033[1;92m", # green
    "y": "\033[1;93m", # yellow
    "b": "\033[1;94m", # blue
    "p": "\033[1;95m", # purple
    "c": "\033[1;96m", # cyan
    "n": "\033[1;0m", # default
}

description = tw.dedent("""\033[1;95m
     _   _                      ____            _        _       
    | | | | ___   ___ _   _ ___|  _ \ ___   ___| | _____| |_ ___ 
    | |_| |/ _ \ / __| | | / __| |_) / _ \ / __| |/ / _ \ __/ __|
    |  _  | (_) | (__| |_| \__ \  __/ (_) | (__|   <  __/ |_\__ \\
    |_| |_|\___/ \___|\__,_|___/_|   \___/ \___|_|\_\___|\__|___/

         Holistic Computational Utility for Scoring Pockets
             
        \033[1;0m""")

epilogue = tw.dedent("""
    Hocus Pockus generates a comprehensive dataframe containing the identified cavities, 
    ranked by their druggability and reactivity scores, providing detailed information 
    on each site's potential as a biological target.
    """)

# Set arguments
def args_parse():
    parser = argparse.ArgumentParser(prog='HocusPockets',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description=description,
                                    epilog=epilogue)
    parser.add_argument(
        '-i',
        '--input',
        required=True, 
        action="store",
        dest="input_path", 
        metavar ="",
        help='Path to a file or directory contaning a protein structure in .pdb format (required!)'
    )
    parser.add_argument(
        '-o',
        '--output',
        metavar="",
        help='Defines the output path'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enables verbose mode'
    )
    parser.add_argument(
        '--distance',
        type=float,
        default=6.0,
        metavar="",
        help='Maximum distance between grid points and residues (default: 6.0)'
    )
    parser.add_argument(
        '--similarity',
        type=float,
        default=0.6,
        metavar="",
        choices=[x / 100 for x in range(1, 101)],
        help='Minimum residue similarity threshold between pockets (default: 0.6)'
    )
    parser.add_argument(
        '--fpocket-score',
        type=float,
        default=0.5,
        metavar="",
        help='Minimum druggability score for Fpocket (default: 0.5).'
    )
    args = parser.parse_args()
    print(description)
    return args

#Verify fpocket program
def check_fpocket_installed(verbose_handler):
    try:
        subprocess.run(['fpocket', '--help'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        verbose_handler.output("'fpocket' is installed and accessible.")
    except subprocess.CalledProcessError:
        print("Error: 'fpocket' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'fpocket' is not installed. Please install it and try again.")
        sys.exit(1)

def check_concavity_installed(verbose_handler):
    try:
        subprocess.run(['concavity'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        verbose_handler.output("'ConCavity' is installed and accessible.")
    except subprocess.CalledProcessError:
        print("Error: 'ConCavity' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'Concavity' is not installed. Please install it and try again.")
        sys.exit(1)

# Verify input
def verify_input(args):
    if os.path.isfile(args.input_path):
        assert args.input_path.endswith(('.pdb')), clrs['r'] + f'The input file "{args.input_file}" must be a .pdb extension.' + clrs['n']
    elif os.path.isdir(args.input_path):
        files = os.listdir(args.input_path)
        pdb_files = [f for f in files if f.endswith(('.pdb'))]
        other_files = [f for f in files if not f.endswith(('.pdb'))]
        assert len(pdb_files) > 0, clrs['r'] + f'The directory "{args.input_path}" does not contain any .pdb files.' + clrs['n']
        for file in other_files:
            print(clrs['y'] + f'Warning: File "{file}" is not a .pdb. Skipping...' + clrs['n'])
    else:
        assert False, clrs['r'] + f'The input path "{args.input_path}" is neither a file nor a directory.' + clrs['n']
    if args.output:
        assert os.path.isdir(args.output), clrs['r'] + f'The output path "{args.output}" does not exist or is not a directory.' + clrs['n']