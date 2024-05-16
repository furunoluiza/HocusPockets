#!/usr/bin/env python

#Imports
import os
import argparse

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
        default='results/',
        metavar="",
        help='Defines the output path'
    )
    args = parser.parse_args()
    return args

# Verify input
def verify_input(args):
    assert os.path.isfile(args.input_file), clrs['r'] + f'The input file "{args.input_file}" does not exist or is not a file.' + clrs['n']
    assert args.input_file.endswith('.pdb'), f'The input file "{args.input_file}" must be a .pdb extension.' + clrs['n']

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