#!/usr/bin/env python3

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

# Imports
###############################################################################
import os
from Initialize import *
from CavityFinder import *
from Reactivity import *

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

class VerboseHandler:
    """
    Handles verbose output and logging to file based on user preference.

    Attributes:
        log_file_path [str]: Path to the log file.
        log_file [file object]: Opened log file for writing log messages.

    Methods:
        output(message): Outputs the message to console or log file depending on verbosity.
        close(): Closes the log file if it was opened.
    """
    def __init__(self, output_path):
        self.log_file_path = os.path.join(output_path, 'verbose.log') if output_path else 'verbose.log'
        self.log_file = open(self.log_file_path, 'w')

    def output(self, message):
        if self.log_file:
            self.log_file.write(message + '\n')

    def close(self):
        if self.log_file:
            self.log_file.close()

# Main Function
def main():
    """
    Main function that orchestrates the HocusPockets pipeline.
    """
    args = args_parse() 
    verbose_handler = VerboseHandler(args.output)
    print(clrs['g'] + "Starting HocusPockets!" + clrs['n'])
    verify_input(args, verbose_handler)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    reactivity_descriptors = Reactivity(args.input_path, output_path, verbose_handler, eband_limit=args.eband_limit)
    reactivity_descriptors.find_reactivity()
    cavity_finder = CavityFinder(args.input_path, output_path, verbose_handler, args.distance, args.similarity, args.fpocket_score)
    cavity_finder.find_cavities()
    print(clrs['g'] + "HocusPockets has completed it`s execution." + clrs['n'])
    verbose_handler.close()

# Execute
if __name__ == "__main__":
    main()