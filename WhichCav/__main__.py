#!/usr/bin/env python3

# Imports
import os
from Initialize import *
from CavityFinder import *

# BINDiS: A Tool for Druggability-Driven and Reactivity-Based Target Prioritization

def main():
    check_fpocket_installed()
    args = args_parse()
    verify_input(args)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    cavity_finder = CavityFinder(args.input_path, output_path)
    cavity_finder.find_cavities()

if __name__ == "__main__":
    main()