#!/usr/bin/env python3

# Description
###############################################################################
'''
Descrição do que esse arquivo faz
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
    def __init__(self, verbose, output_path):
        self.verbose = verbose
        self.log_file_path = os.path.join(output_path, 'verbose.log') if output_path else 'verbose.log'
        self.log_file = open(self.log_file_path, 'w') if not verbose else None

    def output(self, message):
        if self.verbose:
            print(message)
        else:
            if self.log_file:
                self.log_file.write(message + '\n')

    def close(self):
        if self.log_file:
            self.log_file.close()

# Main Function
def main():
    args = args_parse() 
    verbose_handler = VerboseHandler(args.verbose, args.output)
    print(clrs['g'] + "Starting HocusPockets!" + clrs['n'])
    #check_requirements(verbose_handler)
    verify_input(args)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    #reactivity_descriptors = Reactivity(args.input_path, output_path, verbose_handler)
    #reactivity_descriptors.find_reactivity()
    cavity_finder = CavityFinder(args.input_path, output_path, verbose_handler, args.distance, args.similarity, args.fpocket_score)
    cavity_finder.find_cavities()
    print(clrs['g'] + "HocusPockets has completed its execution." + clrs['n'])
    verbose_handler.close()


# Execute
if __name__ == "__main__":
    main()