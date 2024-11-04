#!/usr/bin/env python3

# Imports
import subprocess
import shutil
import os
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull #volume
from Initialize import *

# Hocus Pockus: Holistic Computacional Utility for Scoring Pockets

class CavityFinder:
    def __init__(self, input_path, output_path, verbose_handler, fpocket_score):
        """
        Initializes the CavityFinder object.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities on fpocket.
        Return:
            Nothing.
        """
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler
        self.fpocket_score = fpocket_score
        self.cavity_count = {}
        self.hydrophobic_residues = {'ALA', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'VAL', 'PRO', 'CYS'} #verificar!
        self.polar_residues = {'SER', 'THR', 'CYS', 'ASN', 'GLN', 'ARG', 'LYS', 'HIS', 'ASP', 'GLU'}
    
    @staticmethod
    def extract_residues_from_pdb(cavity_file):
        """
        Extracts unique residues from a PDB file representing a cavity.
        Input:
        cavity_file [str] - Path to the PDB file from which residues are to be extracted.
        Return:
            [set] A set of unique residue identifiers, formatted as "RESNAMERESNUMBER".
        """
        residues = set()
        with open(cavity_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    res_number = line[22:26].strip()
                    residue_id = f"{res_name}{res_number}"
                    residues.add(residue_id)
        return residues

    def compare_residues(self, residues_a, residues_b):
        """
        Compares two sets of residues to determine their similarity.
        Input:
        residues_a [set] - A set of residue identifiers from the first cavity.
        residues_b [set] - A set of residue identifiers from the second cavity.
        Return:
            [bool] True if the similarity between the two sets of residues 
                is greater than or equal to 0.6, otherwise False.
        """
        if not residues_a or not residues_b:
            return False
        intersection = residues_a.intersection(residues_b) # Calculate the intersection of the two sets of residues.
        similarity = len(intersection) / max(len(residues_a), len(residues_b)) # Calculate the similarity as the size of the intersection divided by the size of the larger set.
        return similarity >= 0.8 #Ajustar?

    def calculate_polarity(self, residues):
        """
        Calculates the polarity score based on the presence of polar residues.
        Input:
        residues [set] - A set of residue identifiers associated with the cavity.
        Return:
        [float] The proportion of polar residues relative to the total number of residues.
        """
        if len(residues) == 0:
            return 0
        polar_count = sum(1 for res in residues if res[:3] in self.polar_residues)
        return polar_count / len(residues)

    def calculate_hydrophobicity(self, residues, volume):
        """
        Calculates the hydrophobicity score based on the density of hydrophobic residues.
        Input:
        residues [set]   - A set of residue identifiers associated with the cavity.
        volume   [float] - Volume of the cavity.
        Return:
        [float] The density of hydrophobic residues relative to the volume of the cavity.
        """
        if len(residues) == 0 or volume == 0:
            #print("Warning: No residues or zero volume provided for hydrophobicity calculation.")
            return 0
        hydrophobic_count = sum(1 for res in residues if res[:3] in self.hydrophobic_residues)
        hydrophobicity_density = hydrophobic_count / len(residues) if len(residues) > 0 else 0
        return hydrophobicity_density * volume

    def calculate_cavity_volume(self, lines):
        """
        Calculates the volume of a cavity defined by atomic coordinates from PDB lines.
        Input:
        lines [list] - A list of lines from a PDB file, where each line may represent atomic information.
        Return:
        [float] The calculated volume of the cavity based on its atomic coordinates,
                or 0 if insufficient points are provided or if an error occurs.
        """
        coordinates = []
        for line in lines:
            if line.startswith("ATOM"):  # Check if the line starts with "ATOM" to identify atomic entries
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append([x, y, z])
                except (ValueError, IndexError) as e:
                    print(f"Error processing line: {line.strip()} - {e}")
        if len(coordinates) < 3:  # Check if there are at least three coordinates to define a volume
            #print("Insufficient number of points to calculate the volume.")
            return 0
        hull = ConvexHull(coordinates) # Create a convex hull from the list of coordinates to define the cavity shape.
        return hull.volume if hull else 0 # Return the volume of the convex hull; if the hull could not be created, return 0.
 
    
    def score_cavity(self, pocket_file, residues):
        """
        Calculates a score for a cavity based on volume, hydrophobicity, and frequency of occurrence.
        Parameters:
        pocket_file [str] - Path to the file containing cavity information in PDB format.
        residues    [set] - A set of residues associated with the cavity.
        Return:
        [float] The calculated score for the cavity, which combines volume, hydrophobicity score, 
            and occurrence frequency.
        """
        with open(pocket_file, 'r') as f:
            lines = f.readlines()
        volume = self.calculate_cavity_volume(lines) #Faz sentido calcular o volume com os residuos?
        hydrophobicity_score = self.calculate_hydrophobicity(residues, volume)
        frequency = self.cavity_count.get(pocket_file, 0)
        polarity_score = self.calculate_polarity(residues)
        score = (volume * 0.3) + (hydrophobicity_score * 40) + (frequency * 5) + (polarity_score * 10)
        return score
    
    def save_cavity_scores(self, cavity_path, program_cavities):
        """
        Calculates and saves cavity scores for different programs into a text file.
        Input:
        cavity_path      [str]  - The path to the directory containing the cavity data.
        program_cavities [dict] - A dictionary with program names as keys and a list of cavities as values.
        Return:
            Nothing.
        """
        scores = []
        for program, cavities in program_cavities.items():
            for pocket_file, residues in cavities:
                pocket_file_path = os.path.join(cavity_path, program, pocket_file)
                score = self.score_cavity(pocket_file_path, residues)
                scores.append({
                    'program': program,
                    'pocket_file': pocket_file,
                    'path': pocket_file_path,
                    'score': score
                })
        protein_dir = os.path.join(self.output_path, os.path.basename(cavity_path))
        os.makedirs(protein_dir, exist_ok=True)
        output_file = os.path.join(protein_dir, "cavity_scores.txt")
        with open(output_file, 'w') as f:
            f.write("Program\tPocket\tFile Path\tScore\n")
            for entry in sorted(scores, key=lambda item: item['score'], reverse=True):
                f.write(f"{entry['program']}\t{entry['pocket_file']}\t{entry['path']}\t{entry['score']:.2f}\n")

    def compare_cavities(self, cavity_path):
        """
        Compares cavities across different programs (fpocket, ligsite, pocketfinder, surfnet).
        Input:
        cavity_path [str] - Path to the directory containing the subdirectories with cavities from different programs.
        Return:
            Nothing.
        """
        cavity_directories = ['fpocket', 'ligsite', 'pocketfinder', 'surfnet'] # List of directories corresponding to different cavity-detection programs
        program_cavities = {program: [] for program in cavity_directories}
        for program in cavity_directories: # Loop through each program directory to gather cavity information
            program_dir = os.path.join(cavity_path, program)
            for pocket_file in os.listdir(program_dir): # Iterate over all files in the program's directory
                if pocket_file.endswith('.pdb'):
                    pocket_path = os.path.join(program_dir, pocket_file)
                    residues = self.extract_residues_from_pdb(pocket_path) # Extract residues from the PDB file
                    program_cavities[program].append((pocket_path, residues)) # Store the file path and its residues
                    self.cavity_count[pocket_path] = self.cavity_count.get(pocket_path, 0) + 1
        # Compare cavities from different programs
        for i, prog_a in enumerate(cavity_directories):  # Outer loop for the first program, with index i
            cavities_a = program_cavities[prog_a]
            for prog_b in cavity_directories[i + 1:]:  # Inner loop for the second program, starting from the next program
                cavities_b = program_cavities[prog_b]
                for pocket_a, residues_a in cavities_a:  # Loop through cavities of the first program
                    for pocket_b, residues_b in cavities_b:  # Loop through cavities of the second program
                        if self.compare_residues(residues_a, residues_b):  # Compare the residues from the two cavity files
                            self.cavity_count[pocket_a] = self.cavity_count.get(pocket_a, 0) + 1
                            self.cavity_count[pocket_b] = self.cavity_count.get(pocket_b, 0) + 1
                            cav_key = f"{pocket_a}-{pocket_b}"  # Create a key for the cavity pair
                            print(f"Cavities {os.path.basename(pocket_a)} (from {prog_a}) and {os.path.basename(pocket_b)} (from {prog_b}) are considered equal.")
        self.save_cavity_scores(cavity_path, program_cavities)
    
    def find_cavities_in_directory(self):
        """
        Finds cavities when the input is a directory of files.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        fpocket_finder = FpocketFinder(self.input_path, self.output_path, self.verbose_handler, self.fpocket_score)
        fpocket_out = fpocket_finder.process_directory()
        concavity_finder = ConCavityFinder(self.input_path, self.output_path, self.verbose_handler, self.fpocket_score)
        concavity_finder.process_directory()
        for cavity in fpocket_out:
            self.compare_cavities(cavity)

    def find_cavities_in_file(self, protein_output_dir):
        """
        Finds cavities when the input is a single file.
        Input:
        protein_output_dir [str] - Directory where the protein output will be saved.
        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(self.input_path))[0]
        cavities_output_dir = os.path.join(protein_output_dir, base_name)
        fpocket_finder = FpocketFinder(self.input_path, protein_output_dir, self.verbose_handler, self.fpocket_score)
        fpocket_finder.process_file()
        concavity_finder = ConCavityFinder(self.input_path, protein_output_dir, self.verbose_handler, self.fpocket_score)
        concavity_finder.run_concavity(self.input_path)
        self.compare_cavities(cavities_output_dir)
    
    def find_cavities(self):
        """
        Main method to find cavities.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, base_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        if os.path.isdir(self.input_path):
            self.find_cavities_in_directory()
        else:
            self.find_cavities_in_file(protein_output_dir)

class ConCavityFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, fpocket_score):
        """
        Initializes the ConCavity object, which is a subclass of Cavity.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the Fpocket output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities on fpocket.
        Return:
            Nothing.
        """
        super().__init__(input_path, output_path, verbose_handler, fpocket_score)
    
    def calculate_geometric_center(self, lines):
        """
        Calculates the geometric center of a set of atomic coordinates.
        Input:
            lines [list of str] - A list of lines, where each line contains atomic coordinates.
        Return:
            np.ndarray or None - The geometric center as a numpy array, or None if no coordinates are provided.
        """
        coordinates = []
        for line in lines:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coordinates.append([x, y, z]) # Append the coordinates as a list to the coordinates list
        # Calculate and return the mean coordinates if any were collected
        # This computes the geometric center by averaging the coordinates
        return np.mean(coordinates, axis=0) if coordinates else None

    def process_cavities(self, pockets_info, residue_info, output_path):
        """
        Processes the pocket and residue information to identify and save residues close to each cavity.
        Input:
            pockets_info [str] - Path to the file containing pocket data (HETATM entries).
            residue_info [str] - Path to the file containing residue data (ATOM entries).
            output_path [str]  - Directory where the output files will be saved.
        Return:
            Nothing.
        """
        cavities = {}
        current_cavity = []
        with open(pockets_info, 'r') as file:  # Read the pocket information from the specified pocket file
            for line in file:
                if line.startswith("HETATM"): # Check if the line corresponds to a pocket (HETATM entry)
                    cavity_index = line[22:26].strip() 
                    if cavity_index not in cavities:
                        cavities[cavity_index] = []
                    cavities[cavity_index].append(line)
        residue_coordinates = []
        with open(residue_info, 'r') as file:  # Read the residue information from the specified residue file
            for line in file:
                if line.startswith("ATOM"): # Check if the line corresponds to an atom (ATOM entry)
                    residue_coordinates.append(line)
        # Process each cavity and find residues nearby
        for index, coordinates in cavities.items():
            geom_center = self.calculate_geometric_center(coordinates) # Calculate the geometric center of the cavity
            cavity_volume = self.calculate_cavity_volume(coordinates) # Calculate the volume of the cavity
            if geom_center is not None:
                close_residues = [] # List to hold residues close to the cavity
                for res_line in residue_coordinates:
                    x = float(res_line[30:38].strip())
                    y = float(res_line[38:46].strip())
                    z = float(res_line[46:54].strip())
                    distance = np.linalg.norm(np.array([x, y, z]) - geom_center) # Calculate the distance from the residue to the cavity's geometric center
                    # Define the distance limit based on the cavity's volume
                    distance_limit = max(cavity_volume ** (1/3), 8.0) #isso precisa ser ajustado? #tamanho dos aas
                    if distance < distance_limit: 
                        close_residues.append(res_line)
                output_file = os.path.join(output_path, f"cavity_{index}.pdb") # Create an output file to save residues close to the current cavity
                with open(output_file, 'w') as out_file:
                    out_file.writelines(close_residues)
                
    def access_concavity(self, input_file, method_output_dir, method):
        """
        Accesses and processes the output files generated by the ConCavity analysis for a specific method.
        Input:
            input_file [str]        - Path to the input PDB file used for analysis.
            method_output_dir [str] - Directory where the output files for the specific method are located.
            method [str]            - The grid method used during the ConCavity analysis.
        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        pockets_info = f"{base_name}_{method}_pocket.pdb"
        residue_info = f"{base_name}_{method}_residue.pdb"
        self.process_cavities(pockets_info, residue_info, method_output_dir)
    
    def cleanup_files(self, input_file, method):
        """
        Cleans up temporary files generated during the analysis.
        Input:
            input_file [str] - Path to the input file used.
            method [str]     - The grid method used, which helps identify specific files.
        Return:
            Nothing.
        """
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        files_to_delete = [
            f"{base_name}_{method}_pocket.pdb",
            f"{base_name}_{method}_residue.pdb",
            f"{base_name}_{method}.dx",
            f"{base_name}_{method}.pml"
        ]
        current_directory = os.getcwd()
        for file_name in os.listdir(current_directory):
            if file_name.startswith(f"{base_name}_") and file_name.endswith(f"_{method}.scores"):
                files_to_delete.append(file_name)
        for file_name in files_to_delete:
            file_path = os.path.join(current_directory, file_name)
            if os.path.exists(file_path):
                os.remove(file_path) # Remove the file

    def run_concavity(self, input_file):
        """
        Runs the ConCavity tool on the specified input PDB file using various grid methods.
        Input:
            input_file [str] - Path to the input PDB file for analysis.
        Return:
            Nothing.
        """
        grid_methods = ['pocketfinder', 'ligsite', 'surfnet'] # List of grid methods to be used in the ConCavity analysis
        for method in grid_methods: # Iterate through each grid method to process the input file
            method_output_dir = os.path.join(self.output_path, os.path.splitext(os.path.basename(input_file))[0], method)
            os.makedirs(method_output_dir, exist_ok=True)
            # Construct the command to run the ConCavity tool with specified parameters
            command = [
                'concavity',
                '-grid_method', method,
                '-extraction_method', 'search',
                '-res_map_method', 'blur',
                input_file,
                '-print_grid_pdb', '1',
                method
            ]
            return_code = subprocess.call(command)
            assert return_code == 0, clrs['r'] + f'Problem running ConCavity with {method} method' + clrs['n']
            self.access_concavity(input_file, method_output_dir, method)
            self.cleanup_files(input_file, method)

    def process_directory(self):
        """
        Processes all protein files in a directory for cavity detection using ConCavity.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        for root, _, files in os.walk(self.input_path): # Traverse through all files in the specified input directory
            for file in files:
                if file.endswith('.pdb'):
                    input_file = os.path.join(root, file)
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)
                    self.run_concavity(input_file)

class FpocketFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, fpocket_score):
       """
        Initializes the FpocketFinder object, which is a subclass of Cavity.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the Fpocket output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities based on the druggability score.
        Return:
            Nothing.
        """
       super().__init__(input_path, output_path, verbose_handler, fpocket_score)

    def run_fpocket(self, input_file):
        """
        Executes the Fpocket program on a given protein file for cavity detection.
        Input:
            input_file [str] - Path to the protein file on which Fpocket will be run.
        Return:
            Nothing.
        """
        with open(os.devnull, 'w') as devnull:
            return_code = subprocess.call(['fpocket', '-f' + input_file], stdout=devnull, stderr=devnull)
        assert return_code == 0, clrs['r'] + f'Problem running fpocket tool' + clrs['n']

    @staticmethod
    def parse_info(pockets_info):
        """
        Parses the information from a file that lists properties of various pockets identified in a protein.
        Input:
            pockets_info [str] - Path to the text file containing pocket data.
        Returns:
            [list of dicts] A list of dictionaries, where each dictionary represents an individual pocket and its associated properties.
        """
        # Open the file containing pocket information and read all lines
        with open(pockets_info, 'r') as file: 
            lines = file.readlines()
        pockets = [] 
        current_pocket = {}
        for line in lines:
            # If the line starts with 'Pocket', it signifies the start of a new pocket entry
            if line.startswith('Pocket'):
                if current_pocket:
                    pockets.append(current_pocket)
                current_pocket = {'pocket_number': line.split()[1].strip(':')} # Initialize a new pocket entry with its number
            # If the line contains a colon, it represents a key-value property for the current pocket
            elif ':' in line:
                key, value = line.split(':')
                current_pocket[key.strip()] = value.strip()
        if current_pocket:
            pockets.append(current_pocket)
        return pockets
    
    def druggable_pockets(self, pockets_info, fpocket_score):
        """
        Filters pockets based on druggability score.
        Input:
            pockets_info [str]    - Path to the Fpocket output file containing pocket information.
            fpocket_score [float] - The minimum score threshold to consider a pocket as druggable.
        Return:
            [list of dicts] A list of pockets that have a higher druggability score specified.
        """
        pockets = self.parse_info(pockets_info)
        drug_pockets = [pocket for pocket in pockets if float(pocket['Druggability Score']) > 0.5]
        return drug_pockets

    def access_pockets(self, input_file):
        """
        Retrieves paths to druggable pockets' files.
        Input:
            input_file [str] - Path to the input protein file.

        Return:
            [list of dicts] A list of dictionaries containing paths to each druggable pocket.
        """
        dir_name = os.path.dirname(input_file)
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        extension = os.path.splitext(input_file)[1]
        pockets_info = os.path.join(dir_name, f"{base_name}_out", f"{base_name}_info.txt")
        filtered_pockets = self.druggable_pockets(pockets_info, self.fpocket_score) # Get the filtered list of druggable pockets based on score threshold
        pockets_paths = []
        for pocket in filtered_pockets: # Iterate over each druggable pocket and construct its file path based on the file type
            pocket_number = pocket['pocket_number']
            if extension == ".cif": #Se eu transformar os arquivos pra .pdb não vai precisar
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.cif")
            else:
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_atm.pdb")
            pocket['new_path'] = pocket_path
            pockets_paths.append(pocket)
        return pockets_paths

    def copy_and_cleanup(self, pockets_paths, input_file):
        """
        Copies pocket files to a dedicated directory and cleans up temporary folders after processing.
        Input:
            pockets_paths [list of dicts] - List of dictionaries containing information about each pocket.
            input_file [str] - Path to the protein input file being processed.
            
        Return:
            [list of dicts] List of dictionaries for each pocket.
        """
        updated_pockets_paths = []
        protein_name = os.path.splitext(os.path.basename(input_file))[0]
        protein_output_dir = os.path.join(self.output_path, protein_name)
        cavities_dir = os.path.join(protein_output_dir, 'fpocket')
        os.makedirs(cavities_dir, exist_ok=True)
        if os.path.exists(cavities_dir):
            shutil.rmtree(cavities_dir)  
            os.makedirs(cavities_dir)
        for pocket in pockets_paths:
            pocket_number = pocket['pocket_number']
            original_path = pocket['new_path']
            if os.path.isfile(original_path):
                new_path = os.path.join(cavities_dir, os.path.basename(original_path))
                shutil.copy(original_path, new_path)
                updated_pockets_paths.append({'pocket_number': pocket_number, 'path': new_path})
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        out_dir = os.path.join(os.path.dirname(input_file), f"{base_name}_out")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)
        return updated_pockets_paths

    def process_file(self):
        """
        Processes a single protein file to extract cavities using Fpocket.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        protein_name = os.path.splitext(os.path.basename(self.input_path))[0]
        protein_output_dir = os.path.join(self.output_path, protein_name)
        os.makedirs(protein_output_dir, exist_ok=True)
        self.run_fpocket(self.input_path)
        pockets_info = self.access_pockets(self.input_path)
        update_pockets_paths = self.copy_and_cleanup(pockets_info, self.input_path)
    
    def process_directory(self):
        """
        Processes all protein files in a directory for cavity detection using Fpocket.
        Input:
            Nothing.
        Return:
            [str] Paths to each protein's output directory.
        """
        protein_outputs = []
        for root, _, files in os.walk(self.input_path): # Traverse through all files in the specified input directory
            for file in files:
                if file.endswith(('.pdb', '.cif')):
                    input_file = os.path.join(root, file)
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)
                    self.run_fpocket(input_file)
                    filtered_pockets = self.access_pockets(input_file)
                    updated_pockets_paths = self.copy_and_cleanup(filtered_pockets, input_file)
                    protein_outputs.append(protein_output_dir)
        return protein_outputs
