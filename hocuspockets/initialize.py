
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

import os
import sys
import logging
import argparse
import textwrap as tw
import subprocess
from hocuspockets.constants import PRIMORDIA_BINARY, DEFAULT_MOPAC_FLAGS

# License
###############################################################################
'''
    HocusPockets: Holistic Computacional Utility for Scoring Pockets

    Authors: Machado, L.F.; Torres, P.H.M.; Spelta, G.I.; Oliveira Filho, J.A., Pascutti, P.G.

    [Federal University of Rio de Janeiro]

    Contact info:
    Carlos Chagas Filho Institute of Biophysics
    LMDM - Laboratory for Molecular Modeling and Dynamics
    Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
    Health Sciences Centre - block G1-19
    E-mail address: luizafuruno@ufrj.br
'''

logger = logging.getLogger(__name__)

description = tw.dedent("""\033[1;95m
     _   _                      ____            _        _       
    | | | | ___   ___ _   _ ___|  _ \\ ___   ___| | _____| |_ ___ 
    | |_| |/ _ \\ / __| | | / __| |_) / _ \\ / __| |/ / _ \\ __/ __|
    |  _  | (_) | (__| |_| \\__ \\  __/ (_) | (__|   <  __/ |_\\__ \\\\
    |_| |_|\\___/ \\___|\\__,_|___/_|   \\___/ \\___|_|\\_\\___|\\__|___/

         Holistic Computational Utility for Scoring Pockets
            
        \033[1;0m""")

epilogue = tw.dedent("""
    Hocus Pockets is an automated pipeline that identifies and characterizes protein cavities
    by integrating geometric, physicochemical, and quantum reactivity analyses,
    enabling the ranking of potencial ligand-biding and druggable sites.
    """)

def args_parse() -> argparse.Namespace:
    """
    Parses command-line arguments provided by the user.
    
    Return
    ------
        argparse.Namespace
            A namespace containing parsed argument values.
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
        metavar ="PATH",
        help='Path to a file or directory contaning a protein structure in .pdb format (required!)'
    )
    parser.add_argument(
        '-o', '--output',
        metavar="PATH",
        help='Defines the output path'
    )
    parser.add_argument(
        '--distance',
        type=float,
        default=6.0,
        metavar="VALUE",
        help='Maximum distance between grid points and residues (default: 6.0)'
    )
    parser.add_argument(
        '--similarity',
        type=float,
        default=0.5,
        metavar="VALUE",
        choices=[x / 100 for x in range(1, 101)],
        help='Minimum residue similarity threshold between pockets (default: 0.5)'
    )
    parser.add_argument(
        '--fpocket-score',
        type=float,
        default=0.5,
        metavar="VALUE",
        help='Minimum druggability score for Fpocket (default: 0.5).'
    )
    parser.add_argument(
        '--eband-limit',
        type=int,
        default=3,
        metavar="VALUE",
        help='Sets the HOMO-LUMO energy gap limit (eV) for orbital selection in PRIMoRDIA input (default: 3)'
    )
    parser.add_argument(
        '--parallel',
        action="store_true",
        help='Enable parallel execution for reactivity calculations (default: sequential)'
        )
    parser.add_argument(
        '--n-jobs',
        type=int,
        default=None,
        metavar="VALUE",
        help='Number of processes to use when running in parallel (default: all CPU cores)'
    )
    parser.add_argument(
        '--verbose',
        action="store_true",
        help='Enable verbose logging to a file (verbose.log)'
    )
    parser.add_argument(
        '--mopac-flags',
        type=str,
        default=DEFAULT_MOPAC_FLAGS,
        metavar="FLAGS",
        help=f'MOPAC keyword flags for the .mop header (default: "{DEFAULT_MOPAC_FLAGS}")'
    )
    args = parser.parse_args()
    print(description)
    return args

def check_fpocket_installed(logger: logging.Logger) -> None:
    """
    Checks whether 'fpocket' is installed and accessible from the command line.

    Parameters
    ----------
        logger : logging.Logger
            Object responsible for displaying log messages.
    """
    try:
        subprocess.run(['fpocket', '--help'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info("'fpocket' is installed and accessible.")
    except subprocess.CalledProcessError:
        print("Error: 'fpocket' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'fpocket' is not installed. Please install it and try again.")
        sys.exit(1)

def check_concavity_installed(logger: logging.Logger) -> None:
    """
    Checks whether 'ConCavity' is installed and accessible from the command line.

    Parameters
    ----------
        logger : logging.Logger
            Object responsible for displaying log messages.
    """
    try:
        subprocess.run(['concavity'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info("'ConCavity' is installed and accessible.")
    except subprocess.CalledProcessError:
        print("Error: 'ConCavity' is installed but could not be executed.")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: 'Concavity' is not installed. Please install it and try again.")
        sys.exit(1)

def check_primordia_installed(logger: logging.Logger) -> None:
    """
    Checks whether 'PRIMoRDIA' is installed and accessible from the command line.

    Parameters
    ----------
        logger : logging.Logger
            Object responsible for displaying log messages.
    """
    try:
        result = subprocess.run(
            [PRIMORDIA_BINARY],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if result.returncode == 0:
            logger.info("'PRIMoRDIA' is installed and accessible.")
        else:
            logger.warning("'PRIMoRDIA' is installed but returned a non-zero exit code.")
            logger.warning(f"STDERR: {result.stderr.decode()}")
    except FileNotFoundError:
        print("Error: 'PRIMoRDIA' is not installed. Please install it and try again.")
        sys.exit(1)

def check_mopac_installed(logger: logging.Logger) -> None:
    """
    Checks whether 'MOPAC' is installed and accessible from the expected path.

    Parameters
    ----------
        logger : logging.Logger
            Object responsible for displaying log messages.
    """
    try:
        subprocess.run(
            ['mopac'],
            input=b'\n',
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        logger.info("'MOPAC' is installed and accessible.")
    except FileNotFoundError:
        print("Error: 'MOPAC' is not installed. Please install it and try again.")
        sys.exit(1)

def verify_input(args: argparse.Namespace, logger: logging.Logger) -> None:
    """
    Verifies that input and output paths are valid and that required tools are installed.
    
    Parameters
    ----------
        args : argparse.Namespace
            Parsed command-line arguments containing input/output paths and options.
        logger : logging.Logger
            Object responsible for displaying log messages.
    """
    check_fpocket_installed(logger)
    check_concavity_installed(logger)
    check_primordia_installed(logger)
    check_mopac_installed(logger)
    if os.path.isfile(args.input_path):
        if not args.input_path.endswith(('.pdb')):
            logger.error(f'The input file "{args.input_path}" must be a .pdb extension.')
            sys.exit(1)
    elif os.path.isdir(args.input_path):
        files = os.listdir(args.input_path)
        pdb_files = [f for f in files if f.endswith(('.pdb'))]
        other_files = [f for f in files if not f.endswith(('.pdb'))]
        if len(pdb_files) == 0:
            logger.error(f'The directory "{args.input_path}" does not contain any .pdb files.')
            sys.exit(1)
        for file in other_files:
            logger.warning(f'File "{file}" is not a .pdb. Skipping...')
    else:
        logger.error(f'The input path "{args.input_path}" is neither a file nor a directory.')
        sys.exit(1)
    if args.output:
        if not os.path.isdir(args.output):
            logger.error(f'The output path "{args.output}" does not exist or is not a directory.')
            sys.exit(1)