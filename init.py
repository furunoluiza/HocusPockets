#!/usr/bin/env python

#Imports
import os
import argparse
from utils import *

# BINDiS: Binding Site Detector
# A tool for identify and prioritize bindig sites on proteins
# Developed by Luíza Furuno Machado

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
        '-i',
        '--index',
        action="store",
        dest="index",
        type=int,
        metavar="",
        help='Index of this target protein'
    )
    parser.add_argument(
        '-o',
        '--output',
        default='./results/',
        metavar="",
        help='Defines the output path'
    )
    args = parser.parse_args()
    return args

# Verify input
def verify_input(args):
    assert os.path.isfile(args.input_file), clrs['r'] + f'The input file "{args.input_file}" does not exist or is not a file.' + clrs['n']
    assert args.input_file.endswith('.pdb'), clrs['r'] + f'The input file "{args.input_file}" must be a .pdb extension.' + clrs['n']
    if not os.path.exists(args.output):
        os.makedirs(args.output)
