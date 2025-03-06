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
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, distance_matrix

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
        """
        Initializes the CavityFinder object.
        
        Parameters:
        input_path      (str)  - Path to the input file or directory containing the protein structure.
        output_path     (str)  - Path where the output will be saved.
        verbose_handler (bool) - Object to handle verbose output.
        
        Returns:
        None
        
        Hydrophobicity values are defined based on:
        WIMLEY, W. C.; WHITE, S. H. Experimentally determined hydrophobicity scale for 
        proteins at membrane interfaces. Nature Structural Biology, 1996.
        
        * B (ASX) is assigned the same value as D (ASP).
        * Z (GLX) is assigned the same value as E (GLU).
        """
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

    @staticmethod
    def calculate_geometric_descriptors(lines):
        """
        Calculates the volume, surface area, radius, diameter, smallest size,
        and minimal enclosing cylinder for a cavity defined by atomic coordinates.

        Input:
        lines [list] - A list of lines from a PDB file, where each line may represent atomic information.

        Return:
        [tuple] (volume, surface, radius, diameter, smallest_size, cylinder_radius) - 
                Calculated geometric properties, or zeros if insufficient data is provided.
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

        if len(coordinates) < 4:  # At least four points are needed for a valid convex hull
            return 0, 0, 0, 0, 0, 0

        coordinates = np.array(coordinates)

        # Compute Convex Hull
        hull = ConvexHull(coordinates)
        volume = hull.volume
        surface = hull.area

        # Compute radius of the smallest enclosing sphere
        center = np.mean(coordinates, axis=0)
        distances = np.linalg.norm(coordinates - center, axis=1)
        radius = np.max(distances)

        # Compute diameter (longest pairwise distance)
        pairwise_distances = distance_matrix(coordinates, coordinates)
        diameter = np.max(pairwise_distances)

        # Compute smallest size (distance between closest parallel enclosing planes)
        smallest_size = np.inf
        best_planes = None

        for simplex in hull.simplices:
            # Get the three vertices forming the triangular face
            p1, p2, p3 = coordinates[simplex]

            # Compute normal vector of the triangular face
            normal = np.cross(p2 - p1, p3 - p1)
            normal /= np.linalg.norm(normal)  # Normalize the normal vector

            # Compute plane equation: ax + by + cz + d = 0
            d = -np.dot(normal, p1)

            # Compute distances of all points from this plane
            signed_distances = np.dot(coordinates, normal) + d

            # Identify the extreme points along the normal direction
            min_distance = np.min(signed_distances)
            max_distance = np.max(signed_distances)

            # Compute the plane separation distance
            plane_distance = max_distance - min_distance

            # Check if this is the smallest valid enclosing pair of planes
            if plane_distance < smallest_size:
                smallest_size = plane_distance
                best_planes = (simplex, normal, plane_distance)

        # Compute minimal enclosing cylinder radius
        cylinder_radius = 0
        for axis in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:  # Test cylinders aligned with x, y, z axes
            projection = np.dot(coordinates, axis)
            radial_distances = np.linalg.norm(coordinates - np.outer(projection, axis), axis=1)
            radius = np.max(radial_distances)
            if radius < cylinder_radius or cylinder_radius == 0:
                cylinder_radius = radius

        return volume, surface, radius, diameter, smallest_size, cylinder_radius

    def calculate_residue_descriptors(self, residue_file):
        """
        Calculates the average hydrophobicity of residues of the pocket, 
        along with frequencies of aromatic residues and Tyr(OH) atoms.
        
        Input:
            residue_file [str]: Path to the residue file containing ATOM entries.
        Return:
            [tuple]: A tuple with average hydrophobicity, aromatic residue frequency, 
                    and Tyr(OH) atom frequency.
        """
        total_hydro = 0.0
        residue_count = 0  # Counter for the number of residues
        aromatic_count = 0  # Counter for aromatic residues (F, Y, H, W)
        tyr_oh_count = 0  # Counter for Tyr(OH) atoms

        # Read the residue file and extract residue names and atoms
        with open(residue_file, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    # Extract residue name (columns 17-20 in PDB format)
                    residue_name = line[17:20].strip()
                    # Add hydrophobicity value to the total (default to -0.17 for unknown residues)
                    total_hydro += self.hydro.get(residue_name, self.hydro["UNK"])
                    residue_count += 1  # Increment the residue count

                    # Check for aromatic residues
                    if residue_name in ["PHE", "TYR", "HIS", "TRP"]:
                        aromatic_count += 1

                    # Check for Tyr(OH) atoms (atom name Otyr assumed to be "OH " in columns 13-16)
                    atom_name = line[12:16].strip()
                    if residue_name == "TYR" and atom_name == "OH":
                        tyr_oh_count += 1

        # Calculate average hydrophobicity
        avg_hydro = total_hydro / residue_count if residue_count > 0 else 0.0

        # Calculate frequencies
        aromatic_freq = aromatic_count / residue_count if residue_count > 0 else 0.0
        tyr_oh_freq = tyr_oh_count / residue_count if residue_count > 0 else 0.0

        return avg_hydro, aromatic_freq, tyr_oh_freq

    @staticmethod
    def score_cavity(volume, surface, radius, diameter, smallest_size,
                 cylinder_radius, hydrophobicity_score, aromatic_freq, otyr_freq):
        """
        Calcula o score da cavidade com pesos ajustados e normaliza no intervalo [0,1].
        """
        
        # Pesos ajustados (baseados em uma avaliação dos dados e impacto observado)
        weights = {
            'volume': 0.25,                # Volume teve destaque nos exemplos
            'surface': 0.2,
            'radius': 0.15,
            'diameter': 0.1,
            'smallest_size': 0.1,
            'cylinder_radius': 0.1,
            'hydrophobicity_score': 0.15,
            'aromatic_freq': 0.05,
            'otyr_freq': 0.05
        }

        # Score bruto calculado com base em descritores e seus pesos
        score_raw = (
            weights['volume'] * volume +
            weights['surface'] * surface +
            weights['radius'] * radius +
            weights['diameter'] * diameter +
            weights['smallest_size'] * smallest_size +
            weights['cylinder_radius'] * cylinder_radius +
            weights['hydrophobicity_score'] * hydrophobicity_score +
            weights['aromatic_freq'] * aromatic_freq +
            weights['otyr_freq'] * otyr_freq
        )
        
        # Normalizar com base no valor máximo observado nos dados
        # Podemos calcular dinamicamente o máximo esperado (ou usar uma constante fixada de acordo com o seu caso)
        # Vamos calcular normalizando pelo maior valor de score_raw nos dados

        # Exemplo estático, ajuste com seus próprios dados/experimentos se necessário
        normalization_factor = 700  # Um fator aproximado baseado no máximo dos seus dados de exemplo
        
        # Normalização para manter score entre [0,1]
        score_normalized = min(score_raw / normalization_factor, 1.0)

        return score_raw


    
    def save_cavity_scores(self, cavity_path, scores):
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
        # Write results to the output file
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Cavity_ID', 'Program', 'Grid_File', 'Residue_File', 'Grid_Path', 'Residue_Path',
                        'Volume', 'Radius', 'Surface', 'Diameter', 'Smallest_Size', 'Cylinder_radius', 'Hydrophobicity_Score', 'Aromatic_freq', 'Otyr_freq', 'Score']
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
                    'Volume': f"{entry['volume']:.2f}",  # Arredondado para 2 casas decimais
                    'Radius': f"{entry['radius']:.2f}",  # Arredondado para 2 casas decimais
                    'Surface': f"{entry['surface']:.2f}",  # Arredondado para 2 casas decimais
                    'Diameter': f"{entry['diameter']:.2f}",  # Arredondado para 2 casas decimais
                    'Smallest_Size': f"{entry['smallest_size']:.2f}",  # Arredondado para 2 casas decimais
                    'Cylinder_radius': f"{entry['cylinder_radius']:.2f}",  # Arredondado para 2 casas decimais
                    'Hydrophobicity_Score': f"{entry['hydrophobicity_score']:.2f}",  # Arredondado para 2 casas decimais
                    'Aromatic_freq': f"{entry['aromatic_freq']:.2f}",  # Arredondado para 2 casas decimais
                    'Otyr_freq': f"{entry['tyr_oh_freq']:.2f}",  # Arredondado para 2 casas decimais
                    'Score': f"{entry['score']:.2f}"  # Arredondado para 2 casas decimais
                })

    def score_routine(self, cavity_path, program_cavities, cavity_ids):
        scores = []
        for program, cavities in program_cavities.items():
            for cavity in cavities:
                grid_file_path = cavity['grid_path']
                residue_file_path = cavity['residue_path']
                cavity_id = cavity_ids.get(grid_file_path, None)
                # Read grid lines and calculate volume
                with open(grid_file_path, 'r') as f:
                    grid_lines = f.readlines()
                volume, surface, radius, diameter, smallest_size, cylinder_radius = self.calculate_geometric_descriptors(grid_lines)
                # Calculate other scores (e.g., hydrophobicity, polarity)
                hydrophobicity_score, aromatic_freq, tyr_oh_freq = self.calculate_residue_descriptors(residue_file_path)
                score = self.score_cavity(volume, surface, radius, diameter, smallest_size, cylinder_radius, hydrophobicity_score, aromatic_freq, tyr_oh_freq)
                # Append the results to the scores list
                scores.append({
                    'cavity_id': cavity_id,
                    'program': program,
                    'grid_file': os.path.basename(grid_file_path),
                    'residue_file': os.path.basename(residue_file_path),
                    'grid_path': grid_file_path,
                    'residue_path': residue_file_path,
                    'volume': volume,
                    'radius': radius,
                    'surface': surface,
                    'diameter': diameter,
                    'smallest_size': smallest_size,
                    'cylinder_radius': cylinder_radius,
                    'hydrophobicity_score': hydrophobicity_score,
                    'aromatic_freq': aromatic_freq,
                    'tyr_oh_freq': tyr_oh_freq,
                    'score': score
                })
        self.save_cavity_scores(cavity_path, scores)