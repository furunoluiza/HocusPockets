#!/usr/bin/env python3

# Description
###############################################################################
'''
    This module defines the Reactivity class, which is responsible for calculating
    quantum chemical reactivity descriptors for protein structures using MOPAC
    (semi-empirical quantum chemistry) and PRIMoRDiA (electrostatic and orbital-based analysis).

    Main functionalities:
    - Converts .pdb files into MOPAC-compatible .mop files with appropriate headers.
    - Executes MOPAC and collects output (.aux files) for quantum calculations.
    - Generates PRIMoRDiA input files based on MOPAC results and the original structure.
    - Runs PRIMoRDiA to obtain residue-level electronic descriptors.
    - Supports processing of individual files or batch processing of entire directories.
    - Includes verbose logging support for integration with analysis pipelines.

    This module is typically used as part of the HocusPockets workflow.
'''

# Imports
###############################################################################
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from hocuspockets.initialize import *

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

class Reactivity:
    def __init__(self, input_path, output_path, verbose_handler, eband_limit=3):
        """
        Initializes the Reactivity object for computing reactivity descriptors using
        MOPAC and PRIMoRDiA.

        Parameters
        ----------
            input_path : str
                Path to a single .pdb file or directory containing multiple .pdb files.
            output_path : str
                Directory where output files will be saved.
            verbose_handler : object
                Instance of VerboseHandler for logging output messages.
            eband_limit : int, optional
                Energy range limit (in eV) for orbital selection. Default is 3.
        """
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler
        self.eband_limit = eband_limit
    
    def run_primordia(self, primordia_file, primordia_dir):
        """
        Runs PRIMoRDiA on the given input file.

        Parameters
        ----------
            primordia_file : str
                Path to the PRIMoRDiA input file.
            primordia_dir : str
                Directory where PRIMoRDiA outputs will be saved.
        """
        os.makedirs(primordia_dir, exist_ok=True)
        command = ["PRIMoRDiA_1.25v_LINUX64", "-f", primordia_file]
        self.verbose_handler.output(f"[INFO] Running PRIMoRDiA for: {primordia_file}")
        result = subprocess.run(command, cwd=primordia_dir, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result.returncode != 0:
            self.verbose_handler.output(f"[ERROR] PRIMoRDiA execution failed for: {primordia_file}")
            return
        else:
            self.verbose_handler.output(f"[SUCESS] PRIMoRDIA completed for: {primordia_file}")

    @staticmethod
    def create_primordia_file(input_path, output_aux, primordia_dir, eband_limit=3):
        """
        Creates a PRIMoRDiA input file using the .aux output from MOPAC
        and the original .pdb structure.

        Parameters
        ----------
            input_path : str
                Path to the input PDB file.
            output_aux : str
                Path to the MOPAC .aux output file.
            primordia_dir : str
                Directory where the input file will be saved.
            eband_limit : int, optional
                Energy range limit (in eV) for orbital selection. Default is 3.

        Return
        ------
            str
                Path to the generated PRIMoRDiA input file.
        """
        protein_name = os.path.splitext(os.path.basename(input_path))[0]
        abs_input_path = os.path.abspath(input_path)
        primordia_file = os.path.join(primordia_dir, f"{protein_name}_primordia.input")
        # Write PRIMoRDiA input with required keywords and paths
        with open(primordia_file, "w") as f:
            f.write("#RT normal\n")
            f.write(f"#PR eband {eband_limit} Rscript\n")
            f.write(f"3 {output_aux} none 0 0 {abs_input_path} mopac 0 0 0 0 EW\n")
        return primordia_file

    def run_mopac(self, mopac_file, reactivity_dir):
        """
        Runs MOPAC on the given .mop file and saves the output in the reactivity directory.

        Parameters
        ----------
            mopac_file : str
                Path to the .mop input file.
            reactivity_dir : str
                Directory where MOPAC outputs will be saved.

        Return
        ------
            str
                Path to the generated .aux output file.
        """
        base_name = os.path.splitext(os.path.basename(mopac_file))[0]
        # Define output file paths
        output_aux = os.path.join(reactivity_dir, f"{base_name}.aux")
        output_err = os.path.join(reactivity_dir, f"{base_name}.err")
        self.verbose_handler.output(f"[INFO] Running MOPAC for: {mopac_file}")
        # Run MOPAC and capture errors in a dedicated error file
        with open(output_err, 'w') as f:
            process = subprocess.run(['mopac', mopac_file], cwd=reactivity_dir, 
                                     stdout=subprocess.DEVNULL, stderr=f, check=True)
        return output_aux
        
    @staticmethod
    def create_mopac_file(input_file, protein_output_dir):
        """
        Converts a PDB file into a MOPAC-compatible .mop file, including a predefined header.

        Parameters
        ----------
            input_file : str
                Path to the input PDB file.
            protein_output_dir : str
                Directory where the resulting .mop file will be saved.

        Return
        ------
            str
                Path to the generated .mop file.
        """
        # Define the name of the output .mop file
        mopac_file = os.path.join(protein_output_dir, os.path.basename(input_file).replace(".pdb", ".mop"))
        # Define the MOPAC-specific header (includes PM7 method, solvent model, etc.)
        header = """PM7 1SCF ALLVECS VECTOR AUX LARGE MOZYME eps=78.4
        
        
"""     
        # Extract only ATOM lines from the PDB file
        atom_lines = []
        with open(input_file, "r") as pdb:
            for line in pdb:
                if line.startswith("ATOM"):
                    atom_lines.append(line.strip())
        # Write the .mop file with the header followed by ATOM entries
        with open(mopac_file, "w") as mop:
            mop.write(header)
            mop.write("\n".join(atom_lines))
        return mopac_file
    
    def compress_file(self, file_path):
        """Compress a file with gzip and remove the original."""
        if not os.path.isfile(file_path):
            return None
        gz_path = file_path + ".gz"
        with open(file_path, 'rb') as f_in:
            with gzip.open(gz_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(file_path)
        return gz_path

    def process_directory(self, n_jobs=None, parallel=False):
        """
        Processes all .pdb protein files in the input directory and applies
        MOPAC and PRIMoRDIA calculations on each.

        Parameters
        ----------
            n_jobs : int, optional
                Number of parallel threads to use (default = number of CPU cores).
            parallel : bool, optional
                Whether to run tasks in parallel using threading.
        """
        # Collect all PDB files first
        pdb_files = []
        for root, _, files in os.walk(self.input_path):
            for file in files:
                if file.endswith(".pdb"):
                    pdb_files.append(os.path.join(root, file))
        if not pdb_files:
            self.verbose_handler.output("[WARNING] No PDB files found.")
            return
        def worker(input_file):
            """Worker that runs MOPAC and PRIMoRDIA for one protein file."""
            protein_name = os.path.splitext(os.path.basename(input_file))[0]
            protein_output_dir = os.path.join(self.output_path, protein_name)
            reactivity_dir = os.path.join(protein_output_dir, "reactivity")
            mopac_dir = os.path.join(reactivity_dir, "mopac")
            primordia_dir = os.path.join(reactivity_dir, "primordia")
            try:
                os.makedirs(protein_output_dir, exist_ok=True)
                os.makedirs(mopac_dir, exist_ok=True)
                os.makedirs(primordia_dir, exist_ok=True)
                self.verbose_handler.output(f"[INFO] Processing reactivity for: {input_file}")
                mopac_file = self.create_mopac_file(input_file, mopac_dir)
                output_aux = self.run_mopac(mopac_file, mopac_dir)
                primordia_file = self.create_primordia_file(input_file, output_aux, primordia_dir)
                self.run_primordia(primordia_file, primordia_dir)
                if output_aux and os.path.isfile(output_aux):
                    gz_path = self.compress_file(output_aux)
                    self.verbose_handler.output(f"[INFO] Compressed AUX file -> {gz_path}")
                self.verbose_handler.output(f"[SUCCESS] Finished processing {os.path.basename(input_file)}")
            except Exception as e:
                self.verbose_handler.output(f"[ERROR] Failed to process {input_file}: {e}")
        # Run in parallel with threads
        if parallel:
            if n_jobs is None or n_jobs < 1:
                import multiprocessing as mp
                n_jobs = mp.cpu_count()
            self.verbose_handler.output(f"[INFO] Running in parallel with {n_jobs} threads.")
            with ThreadPoolExecutor(max_workers=n_jobs) as executor:
                futures = {executor.submit(worker, pdb): pdb for pdb in pdb_files}
                for future in as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        pdb = futures[future]
                        self.verbose_handler.output(f"[ERROR] Exception in thread for {pdb}: {e}")
        else:
            self.verbose_handler.output("[INFO] Running sequentially.")
            for pdb in pdb_files:
                worker(pdb)

    def process_file(self, protein_output_dir):
        """
        Processes a single .pdb file provided as input and runs MOPAC and PRIMoRDIA calculations.

        Parameters
        ----------
            protein_output_dir : str
                Path where output directories will be created.
        
        Return
        ------
            None

        Notes
        -----
            This method is used when the input is a single PDB file rather than a directory.
        """
        reactivity_dir = os.path.join(protein_output_dir, "reactivity")
        mopac_dir = os.path.join(reactivity_dir, "mopac")
        primordia_dir = os.path.join(reactivity_dir, "primordia")
        os.makedirs(mopac_dir, exist_ok=True)
        os.makedirs(primordia_dir, exist_ok=True)
        self.verbose_handler.output(f"[INFO] Processing single file for reactivity: {self.input_path}")
        mopac_file = self.create_mopac_file(self.input_path, mopac_dir)
        output_aux = self.run_mopac(mopac_file, mopac_dir)
        primordia_file = self.create_primordia_file(self.input_path, output_aux, primordia_dir)
        self.run_primordia(primordia_file, primordia_dir)
        if output_aux and os.path.isfile(output_aux):
            gz_path = self.compress_file(output_aux)
            self.verbose_handler.output(f"[INFO] Compressed AUX file -> {gz_path}")
        self.verbose_handler.output(f"[SUCCESS] Finished processing single file: {os.path.basename(self.input_path)}")
        
    def find_reactivity(self, parallel=False, n_jobs=None):
        """
        Main method to find cavities in a file or a directory.

        Parameters
        ----------
            parallel : bool, optional
                Whether to run tasks in parallel using threading.
            n_jobs : int, optional
                Number of parallel threads to use.
        
        Return
        ------
            None
        
        Notes
        -----
            This method orchestrates the entire reactivity calculation workflow, handling both
            single file and directory inputs. It creates necessary output directories and
            invokes the appropriate processing methods.

        Notes
        -----
            If the input path is a directory, it will process all .pdb files within that directory.
            If the input path is a single .pdb file, it will process that file directly.
        """
        base_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, base_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        if os.path.isdir(self.input_path):
            self.process_directory(n_jobs=n_jobs, parallel=parallel)
        else:
            self.process_file(protein_output_dir)