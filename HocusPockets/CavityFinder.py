#!/usr/bin/env python3

# Imports
import subprocess
import shutil
import os
import csv
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull #volume
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
        return similarity >= 0.6 #Ajustar? #salvar a equivalência para determinar esse valor

    def plot_hull_and_residues(self, hull_points, residue_points, cavity_index, method):
        """
        Plots the convex hull points of a cavity and the residues within 1Å.
        
        Inputs:
            hull_points [list]: List of points forming the convex hull (N x 3 array).
            residue_points [list]: List of residue coordinates close to the hull (N x 3 array).
            cavity_index [str]: Identifier for the cavity being plotted.
        """
        # Cria uma figura 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plota os pontos do Convex Hull
        hull = ConvexHull(hull_points)
        for simplex in hull.simplices:
            ax.plot(hull_points[simplex, 0], hull_points[simplex, 1], hull_points[simplex, 2], 'k-')

        # Plota todos os pontos da cavidade
        ax.scatter(hull_points[:, 0], hull_points[:, 1], hull_points[:, 2], color='blue', label='Hull Points')

        # Plota os resíduos próximos ao hull
        residue_points = np.array(residue_points)
        if len(residue_points) > 0:
            ax.scatter(residue_points[:, 0], residue_points[:, 1], residue_points[:, 2], color='red', label='Nearby Residues')

        # Configurações do gráfico
        ax.set_title(f"Cavity {cavity_index}_{method} - Convex Hull and Nearby Residues")
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")
        ax.set_zlabel("Z Coordinate")
        ax.legend()
        plt.show()

    def plot_cavity(self, coordinates, pocket_file):
        # Criação do Convex Hull
        hull = ConvexHull(coordinates)
        # Plotagem em 3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        # Plota os pontos atômicos
        coordinates = np.array(coordinates)
        ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], color="blue", alpha=0.5, label="Átomos")
        # Plota as faces do Convex Hull
        for simplex in hull.simplices:
            simplex = np.append(simplex, simplex[0])  # Fecha o loop para o desenho
            ax.plot(coordinates[simplex, 0], coordinates[simplex, 1], coordinates[simplex, 2], "k-")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.title(f"Visualização da Cavidade Convexa - {pocket_file}")
        plt.legend()
        plt.show()

    import numpy as np

    def calculate_centers_of_residues(self, residue_lines):
        """
        Calculates the geometric centers of residues based on ATOM entries.
        Input:
            residue_lines [list]: Lines of the PDB file containing ATOM entries.
        Return:
            [dict]: A dictionary where the keys are residue identifiers 
                and the values are their geometric centers.
        """
        residues = {}
        
        for line in residue_lines:
            # Extract residue information: amino acid name, chain ID, and residue number
            amino_acid = line[17:20].strip()  # Extract amino acid name (e.g., 'ARG')
            chain_id = line[21:22].strip()  # Extract chain id (e.g., 'A')
            residue_number = line[22:26].strip()  # Extract residue number (e.g., '248')
            
            # Form a unique residue identifier: amino acid + chain + residue number
            residue_id = f"{amino_acid} {chain_id} {residue_number}"
            
            # Extract atomic coordinates (x, y, z)
            coords = np.array([float(line[30:38].strip()), 
                            float(line[38:46].strip()), 
                            float(line[46:54].strip())])
            
            # Store the coordinates of each atom in the residue
            if residue_id not in residues:
                residues[residue_id] = []
            residues[residue_id].append(coords)
        
        # Calculate geometric centers
        centers = {}
        for residue_id, atom_coords in residues.items():
            # Calculate the mean of the atom coordinates to find the geometric center
            centers[residue_id] = np.mean(atom_coords, axis=0)
        return centers


    def calculate_cavity_volume(self, lines, pocket_file):
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
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coordinates.append([x, y, z])
                except (ValueError, IndexError) as e:
                    print(f"Error processing line: {line.strip()} - {e}")
        print(f"Number of coordinates for volume calculation: {len(coordinates)}")
        if len(coordinates) < 3:  # Check if there are at least three coordinates to define a volume
            #print("Insufficient number of points to calculate the volume.")
            return 0
        try:
            # Create a convex hull from the list of coordinates
            # The convex hull is the smallest convex shape that encloses all the points
            hull = ConvexHull(coordinates)
            #self.plot_cavity(coordinates, pocket_file) #Plotar cavidades
            return hull.volume  # Return the volume of the convex hull
        except Exception as e:
            print(f"Error creating convex hull: {e}")
            return 0
        
    def calculate_enclosure(self, grid_file, residue_file):
        return 0
    
    def calculate_hydrophilicity(self, grid_file, residue_file):
        hydrophilicity_scores = []
        return 0


    def score_cavity(self, volume, grid_file, residue_file):
        enclosure = self.calculate_enclosure(grid_file, residue_file)
        hydrophilicity = self.calculate_hydrophilicity(grid_file, residue_file)

        drug_score = 0.094 * volume + 0.60 * enclosure - 0.324 * hydrophilicity
        return drug_score, hydrophilicity, enclosure
    
    def save_cavity_scores(self, cavity_path, program_cavities, cavity_ids):
        """
        Calculates and saves cavity scores for different programs into a text file.
        Input:
            cavity_path      [str]  - The path to the directory containing the cavity data.
            program_cavities [dict] - A dictionary with program names as keys and a list of cavities as values.
            cavity_ids       [dict] - A dictionary mapping grid file paths to cavity IDs.
        Return:
            Nothing.
        """
        # Define output path early
        protein_dir = os.path.join(self.output_path, os.path.basename(cavity_path))
        os.makedirs(protein_dir, exist_ok=True)
        output_file = os.path.join(protein_dir, "cavity_scores.csv")
        scores = []
        for program, cavities in program_cavities.items():
            for cavity in cavities:
                grid_file_path = cavity['grid_path']
                residue_file_path = cavity['residue_path']
                cavity_id = cavity_ids.get(grid_file_path, None)

                # Read grid lines and calculate volume
                with open(grid_file_path, 'r') as f:
                    grid_lines = f.readlines()
                volume = self.calculate_cavity_volume(grid_lines, grid_file_path)
                
                # Calculate other scores (e.g., hydrophobicity, polarity)
                score, hydrophilicity_score, enclosure = self.score_cavity(volume, grid_file_path, residue_file_path)

                # Append the results to the scores list
                scores.append({
                    'cavity_id': cavity_id,
                    'program': program,
                    'grid_file': os.path.basename(grid_file_path),
                    'residue_file': os.path.basename(residue_file_path),
                    'grid_path': grid_file_path,
                    'residue_path': residue_file_path,
                    'volume': volume,
                    'hydrophilicity_score': hydrophilicity_score,
                    'enclosure': enclosure,
                    'score': score
                })

        # Write results to the output file
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Cavity_ID', 'Program', 'Grid_File', 'Residue_File', 'Grid_Path', 'Residue_Path',
                        'Volume', 'Hydrophilicity_Score', 'Enclosure', 'Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=';')
            writer.writeheader()
            for entry in sorted(scores, key=lambda item: item['score'], reverse=True):  # Ordena por 'score'
                writer.writerow({
                    'Cavity_ID': entry['cavity_id'],
                    'Program': entry['program'],
                    'Grid_File': entry['grid_file'],
                    'Residue_File': entry['residue_file'],
                    'Grid_Path': entry['grid_path'],
                    'Residue_Path': entry['residue_path'],
                    'Volume': f"{entry['volume']:.2f}",
                    'Hydrophilicity_Score': f"{entry['hydrophilicity_score']:.2f}",
                    'Enclosure': f"{entry['enclosure']:.2f}",
                    'Score': f"{entry['score']:.2f}"
                })


    def compare_cavities(self, cavity_path):
        """
        Compares cavities across different programs (fpocket, ligsite, pocketfinder, surfnet).
        Input:
        cavity_path [str] - Path to the directory containing the subdirectories with cavities from different programs.
        Return:
            Nothing.
        """
        cavity_ids = {}
        current_id = 1
        comparisons = []
        cavity_directories = ['fpocket', 'ligsite', 'pocketfinder', 'surfnet'] # List of directories corresponding to different cavity-detection programs
        program_cavities = {program: [] for program in cavity_directories}
        for program in cavity_directories: # Loop through each program directory to gather cavity information
            program_dir = os.path.join(cavity_path, program)
            for pocket_file in os.listdir(program_dir): # Iterate over all files in the program's directory
                 if pocket_file.startswith("grid_") and pocket_file.endswith(".pdb"):  # Only process grid files
                    pocket_path = os.path.join(program_dir, pocket_file)
                    residue_file = pocket_file.replace("grid_", "residue_")
                    residue_path = os.path.join(program_dir, residue_file)
                    residues = self.extract_residues_from_pdb(residue_path) # Extract residues from the PDB file
                    cavity_id = cavity_ids.get(pocket_path)
                    if cavity_id is None:  # If not, assign a new ID
                        cavity_id = current_id
                        cavity_ids[pocket_path] = current_id
                        current_id += 1
                    # Store the cavity information
                    program_cavities[program].append({
                        'grid_path': pocket_path,
                        'residue_path': residue_path,
                        'residues': residues,
                        'id': cavity_id
                    })
        # Compare cavities from different programs
        for i, prog_a in enumerate(cavity_directories):  # Outer loop for the first program, with index i
            cavities_a = program_cavities[prog_a]
            for prog_b in cavity_directories[i + 1:]:  # Inner loop for the second program, starting from the next program
                cavities_b = program_cavities[prog_b]
                for pocket_a in cavities_a:  # Loop through cavities of the first program
                    for pocket_b in cavities_b:  # Loop through cavities of the second program
                        if self.compare_residues(pocket_a['residues'], pocket_b['residues']):  # Compare the residues from the two cavity files
                            if pocket_a['id'] != pocket_b['id']:
                                print(f"Cavities {os.path.basename(pocket_a['residue_path'])} (from {prog_a}) and {os.path.basename(pocket_b['residue_path'])} (from {prog_b}) are considered equal.")
                            pocket_b['id'] = pocket_a['id']  # Ensure the second cavity has the same ID as the first one
                            cavity_ids[pocket_b['grid_path']] = pocket_a['id']  # Update the cavity ID dictionary
        self.save_cavity_scores(cavity_path, program_cavities, cavity_ids)
    
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

    def process_cavities(self, pockets_info, output_path, method):
        """
        Processes the pocket and residue information to identify and save residues close to each cavity.
        Input:
            pockets_info [str] - Path to the file containing pocket data (HETATM entries).
            output_path [str]  - Directory where the output files will be saved.
            method [str]       - The grid method used during the ConCavity analysis.
        Return:
            Nothing.
        """
        cavities = {}
        with open(pockets_info, 'r') as file:  # Read the pocket information from the specified pocket file
            for line in file:
                if line.startswith("HETATM"):  # Check if the line corresponds to a pocket (HETATM entry)
                    cavity_index = line[22:26].strip()
                    if cavity_index not in cavities:
                        cavities[cavity_index] = []
                    cavities[cavity_index].append(line)

        residue_lines = []
        with open(self.input_path, 'r') as file:  # Read the residue information from the .pdb input file
            for line in file:
                if line.startswith("ATOM"):  # Check if the line corresponds to an atom (ATOM entry)
                    residue_lines.append(line)
        # Calculate the geometric centers of residues
        residue_centers = self.calculate_centers_of_residues(residue_lines)

        for index, coordinates in cavities.items():
            # Save pocket-specific data to a separate file
            pocket_output_file = os.path.join(output_path, f"grid_{index}_{method}.pdb")
            with open(pocket_output_file, 'w') as pocket_file:
                pocket_file.writelines(coordinates)

            if len(coordinates) >= 3:  # Ensure there are enough points to calculate the convex hull
                hull_points = np.array([
                    [float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    for line in coordinates
                ])
                hull = ConvexHull(hull_points)

                # Identify residues close to the convex hull (<= 6.0 Å)
                close_residues = set()  # Set to store unique residues
                for residue_id, center in residue_centers.items():
                    for hull_point in hull.points:
                        if np.linalg.norm(center - hull_point) <= 6.0:  # Distance threshold
                            close_residues.add(residue_id)

                print(f"Residues close to pocket {index}: {close_residues}")

                # Collect ATOM lines of all atoms forming the identified residues
                close_residue_lines = [
                    line for line in residue_lines if f"{line[17:20]} {line[21]} {line[22:26].strip()}" in close_residues
                ]

                print(f"Lines of residues close to the cavity:\n{close_residue_lines}")
                output_file = os.path.join(output_path, f"residue_{index}_{method}.pdb")
                with open(output_file, 'w') as out_file:
                    out_file.writelines(close_residue_lines)
                # Plot the convex hull and nearby residue centers (optional)
                residue_coords = np.array([residue_centers[res_id] for res_id in close_residues])
                #self.plot_hull_and_residues(hull_points, residue_coords, index, method)


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
        self.process_cavities(pockets_info, method_output_dir, method)
    
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

    def process_fpocket_residues(self, pockets_paths, input_file):
        """
        Processes the pocket files generated by fpocket and identifies residues close to each pocket's grid points or geometric centers.
        Input:
            pockets_paths [list of dicts] - List of dictionaries containing pocket paths and information.
            input_file [str] - Path to the input protein file (.pdb format).
        Return:
            None. Residues near pockets are saved to files in the output directory.
        """
        # Read residue lines from the input protein file
        residue_lines = []
        with open(input_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):  # Process only ATOM lines
                    residue_lines.append(line)
        
        # Calculate geometric centers of residues
        residue_centers = self.calculate_centers_of_residues(residue_lines)

        for pocket in pockets_paths:
            pocket_number = pocket['pocket_number']
            pocket_path = pocket['path']  # Path to the pocket grid file

            # Load pocket grid coordinates
            pocket_coordinates = []
            with open(pocket_path, 'r') as pocket_file:
                for line in pocket_file:
                    if line.startswith("ATOM"):  # Consider only ATOM entries
                        x, y, z = map(float, [
                            line[30:38].strip(),
                            line[38:46].strip(),
                            line[46:54].strip()
                        ])
                        pocket_coordinates.append([x, y, z])

            if len(pocket_coordinates) < 3:
                continue  # Skip pockets with fewer than 3 points

            pocket_points = np.array(pocket_coordinates)
            hull = ConvexHull(pocket_points)

            # Identify residues close to the convex hull (<= 6.0 Å)
            close_residues = set()
            for residue_id, center in residue_centers.items():
                for hull_point in hull.points:
                    if np.linalg.norm(center - hull_point) <= 6.0:  # Distance threshold
                        close_residues.add(residue_id)

            # Also check residues near any pocket grid point (<= 4.0 Å)
            for residue_id, center in residue_centers.items():
                for grid_point in pocket_points:
                    if np.linalg.norm(center - grid_point) <= 4.0:
                        close_residues.add(residue_id)

            # Collect ATOM lines of all atoms forming the identified residues
            close_residue_lines = [
                line for line in residue_lines if f"{line[17:20]} {line[21]} {line[22:26].strip()}" in close_residues
            ]

            # Save residue information
            residue_output_file = os.path.join(
            os.path.dirname(pocket_path),  # Same directory as the grid file
            f"residue_{pocket_number}_fpocket.pdb"  # New naming format
        )
            with open(residue_output_file, 'w') as out_file:
                out_file.writelines(close_residue_lines)

            # Plot the convex hull and nearby residue centers (optional)
            residue_coords = np.array([residue_centers[res_id] for res_id in close_residues])
            #self.plot_hull_and_residues(pocket_points, residue_coords, pocket_number, "fpocket")


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
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_vert.pqr")
            else:
                pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_vert.pqr")
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
                new_path = os.path.join(cavities_dir, f"grid_{pocket_number}_fpocket.pdb")
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
        updated_pockets_paths = self.copy_and_cleanup(pockets_info, self.input_path)
        self.process_fpocket_residues(updated_pockets_paths, self.input_path)
    
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
                    self.process_fpocket_residues(updated_pockets_paths, input_file, self.output_path) 
                    protein_outputs.append(protein_output_dir)
        return protein_outputs
