#!/usr/bin/env python3

#Imports
import os
import sys
import argparse
import textwrap as tw
import subprocess

# Hocus Pockus: Holistic Computacional Utility for Scoring Pockets

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
         _   _                        ____            _
        | | | | ___   ___ _   _ ___  |  _ \ ___   ___| | ___   _ ___
        | |_| |/ _ \ / __| | | / __| | |_) / _ \ / __| |/ / | | / __|
        |  _  | (_) | (__| |_| \__ \ |  __/ (_) | (__|   <| |_| \__ \\
        |_| |_|\___/ \___|\__,_|___/ |_|   \___/ \___|_|\_\\__,__|___/
                        
              Holistic Computational Utility for Scoring Pockets
             
        \033[1;0m""")

epilogue = tw.dedent("""
    Hocus Pockus generates a comprehensive dataframe containing the identified cavities, 
    ranked by their druggability and reactivity scores, providing detailed information 
    on each site's potential as a biological target.
    """)

# Set arguments
def args_parse():
    parser = argparse.ArgumentParser(prog='WhichCav',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description=description,
                                    epilog=epilogue)
    parser.add_argument(
        '-i',
        '--input',
        required=True, 
        action="store",
        dest="input_path", 
        metavar = "",
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
        '--fpocket-score',
        type=float,
        default=0.0,
        metavar="",
        help='Minimum druggability score for Fpocket (default: 0.5). Only pockets with a score greater than this will be considered.'
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
        assert args.input_path.endswith(('.pdb', '.cif')), clrs['r'] + f'The input file "{args.input_file}" must be a .pdb or a .cif extension.' + clrs['n']
    elif os.path.isdir(args.input_path):
        files = os.listdir(args.input_path)
        pdb_cif_files = [f for f in files if f.endswith(('.pdb', '.cif'))]
        other_files = [f for f in files if not f.endswith(('.pdb', '.cif'))]
        assert len(pdb_cif_files) > 0, clrs['r'] + f'The directory "{args.input_path}" does not contain any .pdb or .cif files.' + clrs['n']
        for file in other_files:
            print(clrs['y'] + f'Warning: File "{file}" is not a .pdb or .cif. Skipping...' + clrs['n'])
    else:
        assert False, clrs['r'] + f'The input path "{args.input_path}" is neither a file nor a directory.' + clrs['n']
    if args.output:
        assert os.path.isdir(args.output), clrs['r'] + f'The output path "{args.output}" does not exist or is not a directory.' + clrs['n']