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

    @staticmethod
    def run_primordia(primordia_file, primordia_dir):
        """
        Runs PRIMoRDiA on the given input file.

        Input:
            primordia_file [str] - Path to the PRIMoRDiA input file.
            primordia_dir [str]  - Directory where PRIMoRDiA outputs will be saved.

        Return:
            Nothing.
        """
        # Obtém o nome base do arquivo sem extensão
        #protein_name = os.path.splitext(os.path.basename(primordia_file))[0]

        # Comando para executar o PRIMoRDiA
        command = ["/opt/PRIMoRDiA1.0v/PRIMoRDiA_1.25v", "-f", primordia_file, "-verbose"]

        # Garante que o diretório de saída existe
        os.makedirs(primordia_dir, exist_ok=True)

        # Redireciona a saída para um arquivo de log dentro do diretório de saída
        log_file = os.path.join(primordia_dir, f"primordia.log")

        with open(os.devnull, "w") as log:
            return_code = subprocess.call(command, cwd=primordia_dir, stdout=log, stderr=log)

        if return_code != 0:
            raise RuntimeError(f"Error running PRIMoRDiA on {primordia_file}. Check {log_file} for details.")

        print(f"PRIMoRDiA ran successfully. Log saved at {log_file}")
    

    @staticmethod
    def create_primordia_file(input_path, output_aux, primordia_dir):
        """
        Creates a PRIMoRDiA input file inside the 'primordia' directory within the reactivity directory.

        Input:
            input_path [str]     - Path to the input PDB file.
            output_aux [str]     - Path to the MOPAC .aux output file.
            reactivity_dir [str] - Directory where PRIMoRDiA outputs will be stored.

        Return:
            primordia_file [str] - Path to the created PRIMoRDiA input file.
        """
        # Obtém o nome da proteína sem extensão
        protein_name = os.path.splitext(os.path.basename(input_path))[0]

        # Garante que os caminhos são absolutos e formatados corretamente
        abs_input_path = os.path.abspath(input_path)
        
        # Caminho do arquivo de entrada do PRIMoRDiA
        primordia_file = os.path.join(primordia_dir, f"{protein_name}_primordia.input")

        # Escreve o conteúdo no arquivo
        with open(primordia_file, "w") as f:
            f.write("#RT normal\n")
            f.write("#PR eband 3 Rscript\n")
            f.write(f"3 {output_aux} none 0 0 {abs_input_path} mopac 0 0 0 0 EW\n")

        return primordia_file

    @staticmethod
    def run_mopac(mopac_file, reactivity_dir):
        """
        Runs MOPAC on the given .mop file and saves the output in the reactivity directory.

        Input:
            mopac_file [str]     - Path to the .mop input file.
            reactivity_dir [str] - Directory where MOPAC outputs will be saved.

        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(mopac_file))[0]
        output_aux = os.path.join(reactivity_dir, f"{base_name}.aux")
        output_out = os.path.join(reactivity_dir, f"{base_name}.out")
        output_err = os.path.join(reactivity_dir, f"{base_name}.err")

        with open(output_err, 'w') as f:
            process = subprocess.Popen(['/opt/mopac/MOPAC2016.exe', mopac_file], cwd=reactivity_dir, stdout=subprocess.DEVNULL, stderr=f)

        # Aguarda a criação dos arquivos enquanto o processo está rodando
        '''while process.poll() is None:  # process.poll() retorna None se o processo ainda estiver rodando
            if os.path.exists(output_aux) and os.path.exists(output_out):
                process.terminate()  # Encerra o processo se os arquivos foram criados
                return output_aux'''
        return output_aux
        
    @staticmethod
    def create_mopac_file(input_file, protein_output_dir):
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
        header = """PM7 1SCF ALLVECS VECTOR AUX LARGE MOZYME eps=78.4
        
        
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
                    mopac_dir = os.path.join(reactivity_dir, "mopac")
                    os.makedirs(mopac_dir, exist_ok=True)
                    mopac_file = self.create_mopac_file(input_file, mopac_dir)
                    output_aux = self.run_mopac(mopac_file, mopac_dir)
                    primordia_file = self.create_primordia_file(input_file, output_aux, reactivity_dir)
                    

    def process_file(self, protein_output_dir):
        """
        Processes a single protein file to extract ?????.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        reactivity_dir = os.path.join(protein_output_dir, "reactivity")
        mopac_dir = os.path.join(reactivity_dir, "mopac")
        primordia_dir = os.path.join(reactivity_dir, "primordia")
        os.makedirs(mopac_dir, exist_ok=True)
        os.makedirs(primordia_dir, exist_ok=True)
        mopac_file = self.create_mopac_file(self.input_path, mopac_dir)
        output_aux = self.run_mopac(mopac_file, mopac_dir)
        primordia_file = self.create_primordia_file(self.input_path, output_aux, primordia_dir)
        self.run_primordia(primordia_file, primordia_dir)
        
    
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