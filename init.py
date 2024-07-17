#!/usr/bin/env python3

#Imports
import os
import argparse
from main import clrs

# BINDiS: Binding Interaction Druggability Scoring
# A tool for prioritizing biological targets within the proteome 
# of microorganisms based on the druggability of interaction sites

# Set arguments
def args_parse():
    parser = argparse.ArgumentParser(description='BINDiS')
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