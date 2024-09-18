#!/usr/bin/env python3

#Imports
import os
import sys
import argparse
import textwrap as tw
import subprocess

# WhichCav: A Tool for Druggability-Driven and Reactivity-Based Target Prioritization

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
         __          ___     _      _      _____            
         \ \        / / |   (_)    | |    / ____|           
          \ \  /\  / /| |__  _  ___| |__ | |     __ ___   __
           \ \/  \/ / | '_ \| |/ __| '_ \| |    / _` \ \ / /
            \  /\  /  | | | | | (__| | | | |___| (_| |\ V / 
             \/  \/   |_| |_|_|\___|_| |_|\_____\__,_| \_/  
                        
A Tool for Druggability-Driven and Reactivity-Based Target Prioritization
             
        \033[1;0m""")

epilogue = tw.dedent("""
    WhichCav generates a comprehensive dataframe containing the identified cavities, 
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
        '-f',
        '--file',
        required=True, 
        action="store",
        dest="input_file", 
        metavar = "",
        help='Path to a file contaning a protein structure in .pdb format (required!)'
    )
    parser.add_argument(
        '-o',
        '--output',
        metavar="",
        help='Defines the output path'
    )
    args = parser.parse_args()
    print(description)
    return args

#Verify fpocket program
def check_fpocket_installed():
    try:
        subprocess.run(['fpocket', '--help'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("Error: 'fpocket' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'fpocket' is not installed. Please install it and try again.")
        sys.exit(1)

# Verify input
def verify_input(args):
    assert os.path.isfile(args.input_file), clrs['r'] + f'The input file "{args.input_file}" does not exist or is not a file.' + clrs['n']
    assert args.input_file.endswith(('.pdb', '.cif')), clrs['r'] + f'The input file "{args.input_file}" must be a .pdb or a .cif extension.' + clrs['n']
    if args.output:
        assert os.path.isdir(args.output), clrs['r'] + f'The output path "{args.output}" does not exist or is not a directory.' + clrs['n']