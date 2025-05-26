#!/usr/bin/env python3

# Description
###############################################################################
'''
    This script is part of the HocusPockets toolchain for identifying and analyzing 
    druggable pockets in protein structures. It performs the following key tasks:

    - Parses command-line arguments for input/output paths and analysis parameters.
    - Verifies the existence and format of the input (.pdb) files or directories.
    - Checks for the availability and correct installation of required external tools:
    fpocket, ConCavity, PRIMoRDIA, and MOPAC.
    - Logs detailed messages and warnings through a customizable verbose handler.

    This module ensures the environment is correctly configured before the main
    pocket analysis routines are executed.
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

Authors: Machado, L.F.; Pascutti, P.G.; Torres, P.H.M.

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

def args_parse():
    """
    Parses command-line arguments provided by the user.
    Input:
        Nothing.
    Returns:
        args [Namespace] - A namespace containing parsed argument values.
    """
    parser = argparse.ArgumentParser(prog='HocusPockets',
                                    formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description=description,
                                    epilog=epilogue)
    parser.add_argument(
        '-i', '--input',
        required=True, 
        action="store",
        dest="input_path", 
        metavar ="",
        help='Path to a file or directory contaning a protein structure in .pdb format (required!)'
    )
    parser.add_argument(
        '-o', '--output',
        metavar="",
        help='Defines the output path'
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
    parser.add_argument(
        '--eband-limit',
        type=int,
        default=3,
        metavar="",
        help='Sets the HOMO-LUMO energy gap limit (eV) for orbital selection in PRIMoRDIA input (default: 3)'
    )
    args = parser.parse_args()
    print(description)
    return args

def check_fpocket_installed(verbose_handler):
    """
    Checks whether 'fpocket' is installed and accessible from the command line.
    Input:
        verbose_handler - Object responsible for displaying log messages.
    Exit:
        Terminates the script if 'fpocket' is not found or fails to execute.
    """
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
    """
    Checks whether 'ConCavity' is installed and accessible from the command line.
    Input:
        verbose_handler - Object responsible for displaying log messages.
    Exit:
        Terminates the script if 'ConCavity' is not found or fails to execute.
    """
    try:
        subprocess.run(['concavity'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        verbose_handler.output("'ConCavity' is installed and accessible.")
    except subprocess.CalledProcessError:
        print("Error: 'ConCavity' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'Concavity' is not installed. Please install it and try again.")
        sys.exit(1)

def check_primordia_installed(verbose_handler):
    """
    Checks whether 'PRIMoRDIA' is installed and accessible from the command line.
    Input:
        verbose_handler - Object responsible for displaying log messages.
    Exit:
        Terminates the script if 'PRIMoRDIA' is not found.
        Logs a warning if it returns a non-zero exit code.
    """
    try:
        result = subprocess.run(
            ['PRIMoRDiA_1.25v_LINUX64'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if result.returncode == 0:
            verbose_handler.output("'PRIMoRDIA' is installed and accessible.")
        else:
            verbose_handler.output("'PRIMoRDIA' is installed but returned a non-zero exit code.")
            verbose_handler.output(f"STDERR: {result.stderr.decode()}")
    except FileNotFoundError:
        print("Error: 'PRIMoRDIA' is not installed. Please install it and try again.")
        sys.exit(1)

def check_mopac_installed(verbose_handler):
    """
    Checks whether 'MOPAC' is installed and accessible from the expected path.
    Input:
        verbose_handler - Object responsible for displaying log messages.
    Exit:
        Terminates the script if 'MOPAC' is not found or fails to run.
    """
    try:
        subprocess.run(
            ['/opt/mopac/MOPAC2016.exe'],
            input=b'\n',
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        verbose_handler.output("'MOPAC' is installed and accessible.")
    except FileNotFoundError:
        print("Error: 'MOPAC' is not installed. Please install it and try again.")
        sys.exit(1)

def verify_input(args, verbose_handler):
    """
    Verifies that input and output paths are valid and that required tools are installed.
    Input:
        args [Namespace] - Parsed command-line arguments containing input/output paths and options.
        verbose_handler - Object responsible for displaying log messages.
    Exit:
        Terminates the script if:
            - Any required tool is missing.
            - Input path is invalid or lacks .pdb files.
            - Output path is not a directory.
    """
    check_fpocket_installed(verbose_handler)
    check_concavity_installed(verbose_handler)
    check_primordia_installed(verbose_handler)
    check_mopac_installed(verbose_handler)
    if os.path.isfile(args.input_path):
        if not args.input_path.endswith(('.pdb')):
            verbose_handler.output(f'The input file "{args.input_file}" must be a .pdb extension.')
            sys.exit(1)
    elif os.path.isdir(args.input_path):
        files = os.listdir(args.input_path)
        pdb_files = [f for f in files if f.endswith(('.pdb'))]
        other_files = [f for f in files if not f.endswith(('.pdb'))]
        if len(pdb_files) == 0:
            verbose_handler.output(f'The directory "{args.input_path}" does not contain any .pdb files.')
            sys.exit(1)
        for file in other_files:
            verbose_handler.output(f'Warning: File "{file}" is not a .pdb. Skipping...')
            sys.exit(1)
    else:
        verbose_handler.output(f'The input path "{args.input_path}" is neither a file nor a directory.')
        sys.exit(1)
    if args.output:
        if not os.path.isdir(args.output):
            verbose_handler.output(f'The output path "{args.output}" does not exist or is not a directory.')
            sys.exit(1)