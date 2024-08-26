#!/usr/bin/env python3

#Imports
import os
import argparse
from main import clrs

# BINDiS: Best Interaction Network Druggability Score
# A tool for prioritizing biological targets within the proteome 
# of microorganisms based on the druggability of interaction sites

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

description = tw.dedent(clrs[''])

 ____ ___ _   _ ____  _ ____  
| __ )_ _| \ | |  _ \(_) ___| 
|  _ \| ||  \| | | | | \___ \ 
| |_) | || |\  | |_| | |___) |
|____/___|_| \_|____/|_|____/ 

# Set arguments
def args_parse():
    parser = argparse.ArgumentParser(description='BINDiS is a tool to prioritize biological targets within the proteome of microorganisms based on the druggability of interaction sites')
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
    return args

# Verify input
def verify_input(args):
    assert os.path.isfile(args.input_file), clrs['r'] + f'The input file "{args.input_file}" does not exist or is not a file.' + clrs['n']
    assert args.input_file.endswith(('.pdb', '.cif')), clrs['r'] + f'The input file "{args.input_file}" must be a .pdb or a .cif extension.' + clrs['n']
    if args.output:
        assert os.path.isdir(args.output), clrs['r'] + f'The output path "{args.output}" does not exist or is not a directory.' + clrs['n']