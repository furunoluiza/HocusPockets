#!/usr/bin/env python

# Imports
from utils import *

# BINDiS: Binding Site Detector
# A tool for identify and prioritize bindig sites on proteins
# Developed by Luíza Furuno Machado


def main():
    args = args_parse()
    print('Welcome to ' + clrs['p'] + 'BINDiS!' + clrs['n'])
    verify_input(args)

if __name__ == "__main__":
    main()