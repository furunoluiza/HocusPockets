#!/usr/bin/env python3

# Description
###############################################################################
'''
Descrição do que esse arquivo faz
'''

# Imports
###############################################################################
from Initialize import *

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

class Reactivity:
    def __init__(self, input_path, output_path, verbose_handler):
        """
            Initializes the CavityFinder object.
            
            Parameters:
            input_path      (str)  - Path to the input file or directory containing the protein structure.
            output_path     (str)  - Path where the output will be saved.
            verbose_handler (bool) - Object to handle verbose output.
            
            Returns:
            None
        """
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler

    import os

    def create_mopac_file(self, input_file, protein_output_dir):
        """
        Cria um arquivo .mop a partir de um arquivo PDB, incluindo um cabeçalho específico.
        
        Input:
            input_file [str] - Caminho do arquivo PDB de entrada.
            protein_output_dir [str] - Diretório onde o arquivo .mop será salvo.

        Return:
            [str] - Caminho do arquivo .mop gerado.
        """
        # Nome do arquivo de saída
        mopac_file = os.path.join(protein_output_dir, os.path.basename(input_file).replace(".pdb", ".mop"))
        
        # Cabeçalho do arquivo .mop
        header = """PM7 ALLVECS VECTOR AUX LARGE MOZYME eps=78.4 PDB 1SCF
        
        
"""
        
        # Lendo o PDB e extraindo as linhas "ATOM"
        atom_lines = []
        with open(input_file, "r") as pdb:
            for line in pdb:
                if line.startswith("ATOM"):  # Apenas linhas que começam com ATOM
                    atom_lines.append(line.strip())  # Mantém a linha completa

        # Escrevendo o arquivo .mop
        with open(mopac_file, "w") as mop:
            mop.write(header)  # Escreve o cabeçalho
            mop.write("\n".join(atom_lines))  # Escreve as linhas do PDB

        return mopac_file

    def process_directory(self):
        """
        Processes all protein files in a directory ...
        Input:
            Nothing.
        Return:
            
        """
        for root, _, files in os.walk(self.input_path): # Traverse through all files in the specified input directory
            for file in files:
                if file.endswith(('.pdb')):
                    input_file = os.path.join(root, file)
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)
                    reactivity_dir = os.path.join(protein_output_dir, "reactivity")
                    os.makedirs(reactivity_dir, exist_ok=True)
                    self.create_mopac_file(input_file, reactivity_dir)

    def process_file(self, protein_output_dir):
        """
        Processes a single protein file to extract ?????.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        reactivity_dir = os.path.join(protein_output_dir, "reactivity")
        os.makedirs(reactivity_dir, exist_ok=True)
        self.create_mopac_file(self.input_path, reactivity_dir)
        
    
    def find_reactivity(self):
        """
        Main method to find cavities in a file or a directory.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, base_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        if os.path.isdir(self.input_path):
            self.process_directory()
        else:
            self.process_file(protein_output_dir)