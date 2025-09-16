#!/usr/bin/env python3

# Description
###############################################################################
'''
    This module implements the CavityFinder system, which is part of the HocusPockets framework
    for protein pocket detection, characterization, and scoring.

    It integrates multiple tools (Fpocket, ConCavity, Fukui-based analysis) to:
    - Identify potential binding cavities in protein structures.
    - Extract residues near cavity surfaces based on distance to convex hulls.
    - Compare cavity similarity across different methods using residue-based Jaccard index.
    - Generate CSV reports detailing overlap and scores.

    It supports both file and directory-based batch processing and is designed to be modular
    and extensible.
'''

# Imports
###############################################################################
import subprocess
import shutil
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull #volume
import matplotlib.pyplot as plt
from Initialize import *
from ScoringPockets import *
from CatalyticCav import *

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

class CavityFinder:
    def __init__(self, input_path, output_path, verbose_handler, distance, similarity, fpocket_score):
        """
        Initializes the CavityFinder object.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities on fpocket.
        distance        [float] DEFAULT: 6.0 - Maximum allowed distance between grid points and residues when analyzing cavity proximity.
        similarity      [float] DEFAULT: 0.5 - Threshold for residue similarity between cavities. Cavities with a similarity score equal to or above this value are considered similar.
        Return:
            Nothing.
        """
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler
        self.distance = distance
        self.similarity = similarity
        self.fpocket_score = fpocket_score
    
    @staticmethod
    def extract_residues_from_pdb(cavity_file):
        """
        Extracts unique residues from a PDB file representing a cavity.
        Input:
            cavity_file [str] - Path to the PDB file from which residues are to be extracted.
        Return:
            [set] A set of unique residue identifiers, formatted as "chain:resSeqresName".
        """
        residues = set()
        with open(cavity_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    res_name = line[17:20].strip()
                    residue_id = f"{chain_id}:{res_seq}{res_name}"
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
                is greater than or equal to 0.5, otherwise False.
        """
        if not residues_a or not residues_b:
            return False
        intersection = residues_a.intersection(residues_b) # Calculate the intersection of the two sets of residues.
        similarity = len(intersection) / max(len(residues_a), len(residues_b)) # Calculate the similarity as the size of the intersection divided by the size of the larger set.
        return similarity >= self.similarity

    def plot_hull_and_residues(self, hull_points, residue_points, cavity_index, method):
        """
        Plots the convex hull points of a cavity and the residues within a 1Å distance.
        Inputs:
            hull_points [list or np.ndarray]    - Array of 3D coordinates forming the convex hull of the cavity (Nx3).
            residue_points [list or np.ndarray] - Array of 3D coordinates of residues near the cavity surface (within 1Å).
            cavity_index [str]                  - Identifier for the cavity (e.g., index number).
            method [str]                        - Name of the method that generated the cavity (e.g., 'fpocket').
        """
        # Create a 3D figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Compute and plot the edges of the convex hull
        hull = ConvexHull(hull_points)
        for simplex in hull.simplices:
            ax.plot(
                hull_points[simplex, 0],
                hull_points[simplex, 1],
                hull_points[simplex, 2],
                'k-'
            )
        # Plot all cavity (hull) points in blue
        ax.scatter(hull_points[:, 0], hull_points[:, 1], hull_points[:, 2],
                color='blue', label='Hull Points')
        # Plot nearby residues (within 1Å) in red
        residue_points = np.array(residue_points)
        if len(residue_points) > 0:
            ax.scatter(residue_points[:, 0], residue_points[:, 1], residue_points[:, 2],
                    color='red', label='Nearby Residues')
        # Configure axis labels and title
        ax.set_title(f"Cavity {cavity_index}_{method} - Convex Hull and Nearby Residues")
        ax.set_xlabel("X Coordinate")
        ax.set_ylabel("Y Coordinate")
        ax.set_zlabel("Z Coordinate")
        ax.legend()
        # Display the 3D plot
        plt.show()


    def calculate_centers_of_residues(self, residue_lines):
        """
        Calculates the geometric centers of residues based on ATOM entries.
        Input:
            residue_lines [list] - Lines of the PDB file containing ATOM entries.
        Return:
            [dict] A dictionary where the keys are residue identifiers 
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

    def generate_overlap_report(self, cavity_path, program_cavities):
        """
        Generates a CSV file reporting the percentage of residue overlap between all cavities
        detected by different programs (e.g., fpocket, ligsite, surfnet, fukui).
        The report compares each pair of cavities (one from program A, one from program B),
        computes the Jaccard similarity between their residue sets, and outputs the result.
        Input:
            cavity_path [str]       - Path to the parent folder containing cavity subdirectories.
            program_cavities [dict] - Mapping of program name to a list of cavity dictionaries.
                                    Each cavity dict must contain 'residue_path' and 'residues'.

        Output:
            A CSV file named "<protein>_cavity_overlap.csv" saved inside `cavity_path`.
            Each row contains:
                - Cavity A filename
                - Program A
                - Cavity B filename
                - Program B
                - Shared residue percentage (Jaccard similarity × 100)
        """
        # Extract protein name based on the parent folder of the cavity_path
        protein_name = os.path.basename(os.path.dirname(cavity_path))
        # Construct full path to output CSV file
        csv_path = os.path.join(cavity_path, f"{protein_name}_cavity_overlap.csv")
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        # Open the CSV file for writing
        with open(csv_path, mode="w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            # Write the header row
            writer.writerow(["Cavity_A", "Program_A", "Cavity_B", "Program_B", "Shared_Residues_(%)"])
            # Get the list of all program names
            programs = list(program_cavities.keys())
            # Compare each program A with each program B (avoiding duplicates)
            for i, prog_a in enumerate(programs):
                for pocket_a in program_cavities[prog_a]:
                    res_a = set(pocket_a['residues'])  # Set of residues in cavity A
                    for prog_b in programs[i + 1:]:
                        for pocket_b in program_cavities[prog_b]:
                            res_b = set(pocket_b['residues'])  # Set of residues in cavity B
                            # Calculate Jaccard similarity (shared / union)
                            shared = res_a.intersection(res_b)
                            union = res_a.union(res_b)
                            percent = (len(shared) / len(union)) * 100 if union else 0.0
                            # Write one row to CSV
                            writer.writerow([
                                os.path.basename(pocket_a['residue_path']),
                                prog_a,
                                os.path.basename(pocket_b['residue_path']),
                                prog_b,
                                f"{percent:.2f}"
                            ])

    def compare_cavities(self, cavity_path):
        """
        Compares cavities across different programs (fpocket, ligsite, pocketfinder, surfnet, fukui).
        Input:
        cavity_path [str] - Path to the directory containing the subdirectories with cavities from different programs.
        Return:
            Nothing.
        """
        cavity_ids = {}
        current_id = 1
        cavity_directories = ['fpocket', 'ligsite', 'surfnet', 'pocketfinder', 'fukui'] # List of directories corresponding to different cavity-detection programs #retirei pocketfinder
        program_cavities = {program: [] for program in cavity_directories}
        # Step 1: Load all cavity residue sets from each program
        for program in cavity_directories: # Loop through each program directory to gather cavity information
            program_dir = os.path.join(cavity_path, program)
            if not os.path.exists(program_dir):
                continue
            for pocket_file in os.listdir(program_dir): # Iterate over all files in the program's directory
                is_grid = pocket_file.startswith("grid_") and pocket_file.endswith(".pdb")
                is_fukui_residue = program == 'fukui' and pocket_file.startswith("residue_") and pocket_file.endswith(".pdb")
                if is_grid or is_fukui_residue:
                    if program == 'fukui':
                        # For fukui, there's no grid file. Use None explicitly.
                        pocket_path = None
                        residue_path = os.path.join(program_dir, pocket_file)
                    else:
                        pocket_path = os.path.join(program_dir, pocket_file)
                        residue_file = pocket_file.replace("grid_", "residue_")
                        residue_path = os.path.join(program_dir, residue_file)
                    if os.path.exists(residue_path):
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
                #else: se for do metodo de fukui, ele nao tem o grid, mas mudei o nome do arquivo para residue_0_fukui
         # Step 2: Compare cavities between programs
        for i, prog_a in enumerate(cavity_directories):  # Outer loop for the first program, with index i
            cavities_a = program_cavities[prog_a]
            for prog_b in cavity_directories[i + 1:]:  # Inner loop for the second program, starting from the next program
                cavities_b = program_cavities[prog_b]
                for pocket_a in cavities_a:  # Loop through cavities of the first program
                    for pocket_b in cavities_b:  # Loop through cavities of the second program
                        if self.compare_residues(pocket_a['residues'], pocket_b['residues']):  # Compare the residues from the two cavity files
                            if pocket_a['id'] != pocket_b['id']:
                                self.verbose_handler.output(
                                    f"[MATCH] Cavities {os.path.basename(pocket_a['residue_path'])} (from {prog_a}) and {os.path.basename(pocket_b['residue_path'])} (from {prog_b}) share high similarity.")
                            pocket_b['id'] = pocket_a['id']  # Ensure the second cavity has the same ID as the first one
                            cavity_ids[pocket_b['grid_path']] = pocket_a['id']  # Update the cavity ID dictionary
        # Step 3: Generate report and score cavities
        self.generate_overlap_report(cavity_path, program_cavities)
        score_cavity = ScoringPockets(self.input_path, self.output_path, self.verbose_handler)
        score_cavity.score_routine(cavity_path, program_cavities, cavity_ids)
    

    def add_fukushima_scores_to_csv(self, cavity_path, fukushima_results, fukui_results):
        """
        Adds average Fukushima and Fukui_pot_zero scores to an existing cavity score CSV file.
        Input:
            cavity_path [str]        - Path to the directory where the CSV file is located.
            fukushima_results [dict] - Dictionary mapping (program, residue_file) -> average Fukushima score.
            fukui_results [dict]     - Dictionary mapping (program, residue_file) -> average Fukui_pot_zero score.
        Return:
            Nothing.
        """
        # Extract the protein name based on the parent directory name
        protein_name = os.path.basename(os.path.dirname(cavity_path))
        # Construct the full path to the CSV file
        csv_file = os.path.join(cavity_path, f"{protein_name}_cavity_scores.csv")
        updated_rows = []
        # Step 1: Read the existing CSV
        with open(csv_file, 'r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            fieldnames = reader.fieldnames
            # Step 2: Ensure new columns exist
            if 'Fukushima_Score' not in fieldnames:
                fieldnames.append('Fukushima_Score')
            if 'Fukui_Score' not in fieldnames:
                fieldnames.append('Fukui_Score')
            # Step 3: Update each row with scores if available
            for row in reader:
                key = (row['Program'], row['Residue_File'])
                # Get Fukushima score
                fuku_score = fukushima_results.get(key, None)
                row['Fukushima_Score'] = f"{fuku_score:.4f}" if fuku_score is not None else ""
                # Get Fukui score
                fukui_score = fukui_results.get(key, None)
                row['Fukui_Score'] = f"{fukui_score:.4f}" if fukui_score is not None else ""
                updated_rows.append(row)
        # Step 4: Write the updated CSV with the new fields
        with open(csv_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(updated_rows)


    def find_cavities_in_directory(self):
        """
        Finds cavities when the input is a directory of files.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        fpocket_finder = FpocketFinder(self.input_path, self.output_path, self.verbose_handler, self.distance, self.similarity, self.fpocket_score)
        fpocket_finder.process_directory()
        concavity_finder = ConCavityFinder(self.input_path, self.output_path, self.verbose_handler, self.distance, self.similarity, self.fpocket_score)
        protein_output_dirs = concavity_finder.process_directory()
        for protein_output_dir in protein_output_dirs:
            try:
                program_dir = os.path.join(protein_output_dir, 'cavities')  # Adiciona "cavities" ao final
                protein_name = os.path.splitext(os.path.basename(protein_output_dir))[0]
                input_pdb_file = os.path.join(self.input_path, protein_name + ".pdb")
                rslrd_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_residues.rslrd')
                fukui_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukui_pot_zero.pdb')
                fukushima_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukushima.pdb')
                catalytic_cavity = CatalyticCav(program_dir, input_pdb_file, fukui_file, fukushima_file, self.verbose_handler)
                results_fukushima, results_fukui = catalytic_cavity.compute_average()
                catalytic_cavity.process_fukui_cav()
                self.compare_cavities(program_dir)
                self.add_fukushima_scores_to_csv(program_dir, results_fukushima, results_fukui)
            except Exception as e:
                self.verbose_handler.output(f"[WARNING] Error processing: {protein_output_dir}: {e}")
                continue


    def find_cavities_in_file(self, protein_output_dir):
        """
        Finds cavities when the input is a single file.
        Input:
        protein_output_dir [str] - Directory where the protein output will be saved.
        Return:
            Nothing.
        """
        cavities_output_dir = os.path.join(protein_output_dir, 'cavities')
        protein_name = os.path.splitext(os.path.basename(self.input_path))[0]
        rslrd_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_residues.rslrd')
        fukui_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukui_pot_zero.pdb')
        fukushima_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukushima.pdb')
        fpocket_finder = FpocketFinder(self.input_path, protein_output_dir, self.verbose_handler, self.distance, self.similarity, self.fpocket_score)
        fpocket_finder.process_file()
        concavity_finder = ConCavityFinder(self.input_path, protein_output_dir, self.verbose_handler, self.distance, self.similarity, self.fpocket_score)
        concavity_finder.run_concavity(self.input_path, protein_output_dir)
        catalytic_cavity = CatalyticCav(cavities_output_dir, self.input_path, fukui_file, fukushima_file, self.verbose_handler)
        results_fukushima, results_fukui = catalytic_cavity.compute_average()
        catalytic_cavity.process_fukui_cav()
        self.compare_cavities(cavities_output_dir)
        self.add_fukushima_scores_to_csv(cavities_output_dir, results_fukushima, results_fukui)
        
    
    def find_cavities(self):
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
            self.find_cavities_in_directory()
        else:
            self.find_cavities_in_file(protein_output_dir)

class ConCavityFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, distance, similarity, fpocket_score):
        """
        Initializes the ConCavity object, which is a subclass of Cavity.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the Fpocket output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        distance        [float] DEFAULT: 6.0 - Maximum allowed distance between grid points and residues when analyzing cavity proximity.
        similarity      [float] DEFAULT: 0.5 - Threshold for residue similarity between cavities. Cavities with a similarity score equal to or above this value are considered similar.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities on fpocket.
        Return:
            Nothing.
        """
        super().__init__(input_path, output_path, verbose_handler, distance, similarity, fpocket_score)
    
    def process_cavities(self, input_file, pockets_info, output_path, method):
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
        with open(input_file, 'r') as file:  # Read the residue information from the .pdb input file
            for line in file:
                if line.startswith("ATOM"):  # Check if the line corresponds to an atom (ATOM entry)
                    residue_lines.append(line)
        # Calculate the geometric centers of residues
        residue_centers = self.calculate_centers_of_residues(residue_lines)
        for index, coordinates in cavities.items():
            # Save pocket-specific data to a separate file
            pocket_output_file = os.path.join(output_path, f"grid_{index}_{method}.pdb")
            with open(pocket_output_file, 'w') as pocket_file:
                for line in coordinates:
                    # Extract atom index and residue index
                    atom_index = int(line[6:11].strip())  # Extract atom index (e.g., 0) as integer
                    residue_index = int(line[22:26].strip())  # Extract residue index (e.g., 0) as integer
                    
                    # Combine atom and residue indices into a single value
                    combined_index = f"{atom_index}_{residue_index}"
                    
                    # Extract coordinates
                    x = float(line[30:38].strip())  # X coordinate
                    y = float(line[38:46].strip())  # Y coordinate
                    z = float(line[46:54].strip())  # Z coordinate
                    
                    # Extract the last two columns (occupancy and temperature factor)
                    occupancy = line[54:60].strip()
                    temperature_factor = line[60:66].strip()
                    
                    # Format the line in the same way as the ATOM entry
                    new_line = f"ATOM  {atom_index:5d}    H POC {residue_index:5d}     {x:7.3f} {y:7.3f} {z:7.3f}    {occupancy:6s}{temperature_factor:6s}\n"

                    # Write the new line to the file
                    pocket_file.write(new_line)
                pocket_file.write("TER\n")
                pocket_file.write("END\n")
            if len(coordinates) >= 3:  # Ensure there are enough points to calculate the convex hull
                hull_points = np.array([[float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())]
                    for line in coordinates])
                hull = ConvexHull(hull_points)
                # Identify residues close to the convex hull (<= 6.0 Å)
                close_residues = set()  # Set to store unique residues
                for residue_id, center in residue_centers.items():
                    for hull_point in hull.points:
                        if np.linalg.norm(center - hull_point) <= self.distance:  # Distance threshold
                            close_residues.add(residue_id)
                # Collect ATOM lines of all atoms forming the identified residues
                close_residue_lines = [line for line in residue_lines if f"{line[17:20]} {line[21]} {line[22:26].strip()}" in close_residues]
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
        self.process_cavities(input_file, pockets_info, method_output_dir, method)
    
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

    def run_concavity(self, input_file, output_path):
        """
        Runs the ConCavity tool on the specified input PDB file using various grid methods.
        Input:
            input_file [str] - Path to the input PDB file for analysis.
        Return:
            Nothing.
        """
        grid_methods = ['pocketfinder', 'ligsite', 'surfnet'] # List of grid methods to be used in the ConCavity analysis
        for method in grid_methods: # Iterate through each grid method to process the input file
            method_output_dir = os.path.join(output_path, 'cavities', method)
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
        output_dirs = []
        for root, _, files in os.walk(self.input_path): # Traverse through all files in the specified input directory
            for file in files:
                if file.endswith('.pdb'):
                    input_file = os.path.join(root, file)
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)
                    self.run_concavity(input_file, protein_output_dir)
                    output_dirs.append(protein_output_dir)
        return (output_dirs)

class FpocketFinder(CavityFinder):
    def __init__(self, input_path, output_path, verbose_handler, distance, similarity, fpocket_score):
       """
        Initializes the FpocketFinder object, which is a subclass of Cavity.
        Input:
        input_path      [str]                - Path to the input file or directory containing the protein structure.
        output_path     [str]                - Path where the Fpocket output will be saved.
        verbose_handler [bool]               - Object to handle verbose output.
        distance        [float] DEFAULT: 6.0 - Maximum allowed distance between grid points and residues when analyzing cavity proximity.
        similarity      [float] DEFAULT: 0.5 - Threshold for residue similarity between cavities. Cavities with a similarity score equal to or above this value are considered similar.
        fpocket_score   [float] DEFAULT: 0.5 - Threshold score for selecting cavities based on the druggability score.
        Return:
            Nothing.
        """
       super().__init__(input_path, output_path, verbose_handler, distance, similarity, fpocket_score)

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
        assert return_code == 0, clrs['r'] + f'Problem running fpocket on {input_file}' + clrs['n']

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
        drug_pockets = [pocket for pocket in pockets if float(pocket['Druggability Score']) > fpocket_score]
        return drug_pockets

    def process_fpocket_residues(self, pockets_paths, input_file):
        """
        Processes the pocket files generated by fpocket and identifies residues close to each pocket's grid points or geometric centers.
        Input:
            pockets_paths [list of dicts] - List of dictionaries containing pocket paths and information.
            input_file [str] - Path to the input protein file (.pdb format).
        Return:
            Nothing.
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
                        x, y, z = map(float, [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()])
                        pocket_coordinates.append([x, y, z])
            if len(pocket_coordinates) < 3:
                continue  # Skip pockets with fewer than 3 points
            pocket_points = np.array(pocket_coordinates)
            hull = ConvexHull(pocket_points)
            # Identify residues close to the convex hull (<= 6.0 Å)
            close_residues = set()
            for residue_id, center in residue_centers.items():
                for hull_point in hull.points:
                    if np.linalg.norm(center - hull_point) <= self.distance:  # Distance threshold
                        close_residues.add(residue_id)
            # Collect ATOM lines of all atoms forming the identified residues
            close_residue_lines = [line for line in residue_lines if f"{line[17:20]} {line[21]} {line[22:26].strip()}" in close_residues]
            # Save residue information
            residue_output_file = os.path.join(os.path.dirname(pocket_path), f"residue_{pocket_number}_fpocket.pdb")
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
            pocket_path = os.path.join(dir_name, f"{base_name}_out", "pockets", f"pocket{pocket_number}_vert.pqr")
            pocket['new_path'] = pocket_path
            pockets_paths.append(pocket)
        return pockets_paths

    def copy_and_cleanup(self, pockets_paths, input_file, output_path):
        """
        Copies pocket files to a dedicated directory and cleans up temporary folders after processing.
        Input:
            pockets_paths [list of dicts] - List of dictionaries containing information about each pocket.
            input_file [str] - Path to the protein input file being processed.
        Return:
            [list of dicts] List of dictionaries for each pocket.
        """
        updated_pockets_paths = []
        protein_output_dir = os.path.join(output_path, 'cavities')
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
        self.run_fpocket(self.input_path)
        pockets_info = self.access_pockets(self.input_path)
        updated_pockets_paths = self.copy_and_cleanup(pockets_info, self.input_path, self.output_path)
        self.process_fpocket_residues(updated_pockets_paths, self.input_path)
    
    def process_directory(self):
        """
        Processes all protein files in a directory for cavity detection using Fpocket.
        Input:
            Nothing.
        Return:
            [str] Paths to each protein's output directory.
        """
        #protein_outputs = []
        for root, _, files in os.walk(self.input_path): # Traverse through all files in the specified input directory
            for file in files:
                if file.endswith(('.pdb')):
                    input_file = os.path.join(root, file)
                    protein_name = os.path.splitext(os.path.basename(input_file))[0]
                    protein_output_dir = os.path.join(self.output_path, protein_name)
                    os.makedirs(protein_output_dir, exist_ok=True)
                    self.run_fpocket(input_file)
                    filtered_pockets = self.access_pockets(input_file)
                    updated_pockets_paths = self.copy_and_cleanup(filtered_pockets, input_file, protein_output_dir)
                    self.process_fpocket_residues(updated_pockets_paths, input_file) 
                    #protein_outputs.append(protein_output_dir)
        #return protein_outputs
