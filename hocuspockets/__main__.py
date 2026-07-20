
import logging
import os
import sys
from typing import Optional

from hocuspockets.initialize import args_parse, verify_input
from hocuspockets.cavity_finder import CavityFinder
from hocuspockets.reactivity import Reactivity
from hocuspockets.color import clrs

# Description
###############################################################################
'''
    This is the main execution script for the HocusPockets pipeline.

    It includes the VerboseHandler class for managing logging and output behavior,
    and defines the main() function, which serves as the entry point for the workflow.

    Main responsibilities:
    - Parse and validate command-line arguments.
    - Check for required software dependencies.
    - Initialize the logging system.
    - Run the reactivity descriptor and cavity detection modules.
    - Provide user feedback via verbose or log file output.

    This script is intended to be run directly and should not be imported as a module.
'''

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


def setup_logging(output_path: Optional[str] = None, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger('hocuspockets')
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if verbose and output_path:
        fh = logging.FileHandler(os.path.join(output_path, 'hocuspockets.log'))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(fh)

    return logger


def main() -> None:
    """
    Main function that orchestrates the HocusPockets pipeline.
    """
    args = args_parse()
    logger = setup_logging(args.output, verbose=args.verbose)

    print(clrs['g'] + "Starting HocusPockets!" + clrs['n'])
    verify_input(args, logger)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    reactivity_descriptors = Reactivity(args.input_path, output_path, logger, eband_limit=args.eband_limit)
    reactivity_descriptors.find_reactivity(parallel=args.parallel, n_jobs=args.n_jobs)
    cavity_finder = CavityFinder(args.input_path, output_path, logger, args.distance, args.similarity, args.fpocket_score)
    cavity_finder.find_cavities()
    print(clrs['g'] + "HocusPockets has completed its execution." + clrs['n'])


if __name__ == "__main__":
    main()
