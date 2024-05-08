#!/usr/bin/env python

import os
import argparse

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
        #Vai ser um indice de 0-1?
        #Esse arg vai ser required?
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
    if not os.path.isfile(args.input_file):
        print(f'Error: The input file "{args.input_file}" does not exist or is not a file.')
        quit()

def main():
    args = args_parse()
    print('Welcome to the BINDiS!')
    verify_input(args)

if __name__ == "__main__":
    main()