#!/usr/bin/env python3

# Description
###############################################################################
'''
Descrição do que esse arquivo faz
'''

# Imports
###############################################################################
import os
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull #volume

# License
###############################################################################
'''
HocusPockets: Holistic Computacional Utility for Scoring Pockets

Authors: Machado, L.F.; Torres, P.H.M.

[Federal University of Rio de Janeiro]

Contact info:
Carlos Chagas Filho Institute of Biophysics
LMDM - Laboratory for Molecular Modeling and Dynamics
Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
Health Sciences Centre - block G1-19
E-mail address: luizafuruno@ufrj.br
'''

class ScoringPockets:
    def __init__(self, input_path, output_path, verbose_handler):
        self.input_path = input_path
        self.output_path = output_path
        self.verbose_handler = verbose_handler
        self.hydro = {"ALA": -0.17, "ARG": -0.81, "ASN": -0.42, "ASP": -1.23, "CYS": 0.24,
                "GLN": -2.02, "GLU": -0.58, "GLY": -0.01, "HIS": -0.96, "ILE": 0.31,
                "LEU": 0.56, "LYS": -0.99, "MET": 0.23, "PHE": 1.13, "PRO": -0.45,
                "SER": -0.13, "THR": -0.14, "TRP": 1.85, "TYR": 0.94, "VAL": -0.07,
                "MME": 0.23, "MSE": 0.23, "UNK": -0.17}
        
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
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coordinates.append([x, y, z])
                    except (ValueError, IndexError) as e:
                        print(f"Error processing line: {line.strip()} - {e}")
            if len(coordinates) < 3:  # Check if there are at least three coordinates to define a volume
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
    
    def calculate_hydro(self, residue_file):
        """
        Calculates the total hydrophobicity of residues of the pocket.
        Input:
            residue_file [str]: Path to the residue file containing ATOM entries.
        Return:
            [float] Total hydrophobicity of residues in the residue file.
        """
        total_hydro = 0.0
        # Read the residue file and extract residue names
        with open(residue_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    # Extract the residue name (columns 17-20 in PDB format)
                    residue_name = line[17:20].strip()
                    # Add hydrophobicity value to the total (default to -0.17 for unknown residues)
                    total_hydro += self.hydro.get(residue_name, self.hydro["UNK"])
        return total_hydro

    def score_cavity(self, volume, grid_file, residue_file):
        enclosure = self.calculate_enclosure(grid_file, residue_file)
        hydro = self.calculate_hydro(residue_file)
        drug_score = 0.094 * volume + 0.60 * enclosure + 0.324 * hydro
        return drug_score, hydro, enclosure
    
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
                volume = self.calculate_cavity_volume(grid_lines)
                # Calculate other scores (e.g., hydrophobicity, polarity)
                score, hydrophobicity_score, enclosure = self.score_cavity(volume, grid_file_path, residue_file_path)
                # Append the results to the scores list
                scores.append({
                    'cavity_id': cavity_id,
                    'program': program,
                    'grid_file': os.path.basename(grid_file_path),
                    'residue_file': os.path.basename(residue_file_path),
                    'grid_path': grid_file_path,
                    'residue_path': residue_file_path,
                    'volume': volume,
                    'hydrophobicity_score': hydrophobicity_score,
                    'enclosure': enclosure,
                    'score': score
                })
        # Write results to the output file
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Cavity_ID', 'Program', 'Grid_File', 'Residue_File', 'Grid_Path', 'Residue_Path',
                        'Volume', 'Hydrophobicity_Score', 'Enclosure', 'Score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=',')
            writer.writeheader()
            for entry in sorted(scores, key=lambda item: item['score'], reverse=True):  # Sort by 'score'
                writer.writerow({
                    'Cavity_ID': entry['cavity_id'],
                    'Program': entry['program'],
                    'Grid_File': entry['grid_file'],
                    'Residue_File': entry['residue_file'],
                    'Grid_Path': entry['grid_path'],
                    'Residue_Path': entry['residue_path'],
                    'Volume': f"{entry['volume']:.2f}",
                    'Hydrophobicity_Score': f"{entry['hydrophobicity_score']:.2f}",
                    'Enclosure': f"{entry['enclosure']:.2f}",
                    'Score': f"{entry['score']:.2f}"
                })