#!/usr/bin/env python3

# Imports
import os
from Initialize import *
from CavityFinder import *

# Hocus Pockus: Holistic Computacional Utility for Scoring Pockets
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

def main():
    args = args_parse()
    verbose_handler = VerboseHandler(args.verbose, args.output)
    print(clrs['g'] + "Starting WhichCav!" + clrs['n'])
    check_fpocket_installed(verbose_handler)
    check_concavity_installed(verbose_handler)
    verify_input(args)
    output_path = os.path.abspath(args.output) if args.output else os.getcwd()
    cavity_finder = CavityFinder(args.input_path, output_path, verbose_handler, args.druggability_score)
    cavity_finder.find_cavities()
    print(clrs['g'] + "WhichCav has completed its execution." + clrs['n'])
    verbose_handler.close()

if __name__ == "__main__":
    main()