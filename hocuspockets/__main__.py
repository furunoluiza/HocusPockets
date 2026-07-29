import os

from hocuspockets import setup_logging
from hocuspockets.initialize import args_parse, verify_input
from hocuspockets.cavity_finder import CavityFinder
from hocuspockets.reactivity import Reactivity
from hocuspockets.color import clrs

# Description
###############################################################################
"""
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
"""

# License
###############################################################################
"""
    HocusPockets: Holistic Computacional Utility for Scoring Pockets

    Authors: Machado, L.F.; Torres, P.H.M.; Spelta, G.I.; Oliveira Filho, J.A., Pascutti, P.G.

    [Federal University of Rio de Janeiro]

    Contact info:
    Carlos Chagas Filho Institute of Biophysics
    LMDM - Laboratory for Molecular Modeling and Dynamics
    Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
    Health Sciences Centre - block G1-19
    E-mail address: luizafuruno@ufrj.br
"""


def main() -> None:
    """
    Main function that orchestrates the HocusPockets pipeline.
    """
    args = args_parse()
    logger = setup_logging(args.output, verbose=args.verbose)

    logger.info(clrs["g"] + "Starting HocusPockets!" + clrs["n"])
    verify_input(args, logger)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    reactivity_descriptors = Reactivity(
        input_path = args.input_path,
        output_path = output_path,
        logger = logger,
        eband_limit=args.eband_limit,
        mopac_flags=args.mopac_flags,
    )
    reactivity_descriptors.find_reactivity(parallel=args.parallel, n_jobs=args.n_jobs)
    cavity_finder = CavityFinder(
        input_path = args.input_path,
        output_path = output_path,
        logger = logger,
        distance = args.distance,
        similarity = args.similarity,
        fpocket_score = args.fpocket_score,
        toggle_image = args.toggle_image
    )
    cavity_finder.find_cavities()
    logger.info(clrs["g"] + "HocusPockets has completed its execution." + clrs["n"])


if __name__ == "__main__":
    main()
