# Description
###############################################################################
"""
This module defines the ScoringPockets class, which is responsible for computing
a variety of physicochemical and geometric descriptors for protein cavities
detected by multiple methods.

Key functionalities include:
- Computing geometric properties from convex hulls: volume, surface area, radius,
diameter, minimal thickness, and cylindrical fit.
- Calculating residue-based descriptors: average hydrophobicity, frequency of
aromatic residues, and the presence of catalytic Tyr(OH) groups.
- Aggregating these features into a composite cavity score (currently experimental).
- Saving scores to structured CSV reports for further analysis or visualization.

The scoring system is designed to assist in evaluating and ranking potential
ligand-binding pockets in protein structures.

Note: The scoring weights and normalization strategy are still subject to validation.
"""

import os
import csv
import logging
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull, distance_matrix
from hocuspockets.constants import HYDROPHOBICITY_SCALE, AROMATIC_RESIDUES
from hocuspockets import pdb_utils

logger = logging.getLogger(__name__)

# License
###############################################################################
"""
    HocusPockets: Holistic Computacional Utility for Scoring Pockets

    Authors: Machado, L.F.; Torres, P.H.M.; Spelta, G.I.; Oliveira Filho, J.A., Pascutti, P.G.

    [Federal University of Rio de Janeiro]

    Contact info:
    Carlos Chagas Filho Institute of Biophysics
    LMDM - Laboratory for Molecular Modeling and Dynamics
    Carlos Chagas Filho Av. 373 - Rio de Janeiro - RJ - Brazil
    Health Sciences Centre - block G1-19
    E-mail address: luizafuruno@ufrj.br
"""


class ScoringPockets:
    def __init__(
        self, input_path: str, output_path: str, logger: logging.Logger
    ) -> None:
        """
        Initializes the CavityFinder object.

        Hydrophobicity values are defined based on:
        WIMLEY, W. C.; WHITE, S. H. Experimentally determined hydrophobicity scale for
        proteins at membrane interfaces. Nature Structural Biology, 1996.

        Parameters
        ----------
            input_path : str
                Path to the input file or directory containing the protein structure.
            output_path : str
                Path where the output will be saved.
            logger : logging.Logger
                Logger object to handle verbose output.

        Notes
        -----
            B (ASX) is assigned the same value as D (ASP).
            Z (GLX) is assigned the same value as E (GLU).
        """
        self.input_path = input_path
        self.output_path = output_path
        self.logger = logger
        self.hydro = HYDROPHOBICITY_SCALE

    def plot_cavity(self, coordinates: List[List[float]], pocket_file: str) -> None:
        """
        Plots a 3D representation of a cavity including its convex hull and atoms.

        Parameters
        ----------
            coordinates : list
                List of 3D coordinates representing cavity atoms.
            pocket_file : str
                Name of the pocket file for the plot title.
        """
        hull = ConvexHull(coordinates)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        # Plot the atomic points
        coordinates = np.array(coordinates)
        ax.scatter(
            coordinates[:, 0],
            coordinates[:, 1],
            coordinates[:, 2],
            color="blue",
            alpha=0.5,
            label="Átomos",
        )
        # Plot the convex hull edges
        for simplex in hull.simplices:
            simplex = np.append(simplex, simplex[0])
            ax.plot(
                coordinates[simplex, 0],
                coordinates[simplex, 1],
                coordinates[simplex, 2],
                "k-",
            )
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.title(f"Convex Cavity Visualization - {pocket_file}")
        plt.legend()
        plt.show()

    @staticmethod
    def calculate_geometric_descriptors(
        lines: List[str], grid_file_path: str
    ) -> Tuple[float, float, float, float, float, float]:
        """
        Calculates the volume, surface area, radius, diameter, smallest size,
        and minimal enclosing cylinder for a cavity defined by atomic coordinates.

        Parameters
        ----------
            lines : list
                A list of lines from a PDB file, where each line may represent atomic information.
            grid_file_path : str
                Path to the grid file for error reporting.

        Return
        ------
            tuple
                (volume, surface, radius, diameter, smallest_size, cylinder_radius)
                Calculated geometric properties, or zeros if insufficient data is provided.
        """
        coordinates = pdb_utils.extract_atom_coordinates(lines)
        if len(coordinates) < 4:
            return 0, 0, 0, 0, 0, 0
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
        # Compute minimal enclosing cylinder radius
        cylinder_radius = 0
        for axis in [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
        ]:  # Test cylinders aligned with x, y, z axes
            projection = np.dot(coordinates, axis)
            radial_distances = np.linalg.norm(
                coordinates - np.outer(projection, axis), axis=1
            )
            radius = np.max(radial_distances)
            if radius < cylinder_radius or cylinder_radius == 0:
                cylinder_radius = radius
        return volume, surface, radius, diameter, smallest_size, cylinder_radius

    def calculate_residue_descriptors(
        self, residue_file: str
    ) -> Tuple[float, float, float, float]:
        """
        Calculates the average hydrophobicity of residues of the pocket,
        along with frequencies of aromatic residues and Tyr(OH) atoms.

        Parameters
        ----------
            residue_file : str
                Path to the residue file containing ATOM entries.

        Return
        ------
            tuple
                A tuple with average hydrophobicity, aromatic residue frequency,
                Tyr(OH) atom frequency, and carbon atom frequency.
        """
        total_hydro = 0.0
        residue_count = 0
        aromatic_count = 0
        tyr_oh_count = 0
        total_atom_count = 0
        carbon_atom_count = 0
        residues_seen = set()
        aromatic_residues = AROMATIC_RESIDUES
        with open(residue_file, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    residue_number = line[22:26].strip()
                    residue_id = (residue_name, residue_number, chain_id)
                    atom_name = line[12:16].strip()
                    total_atom_count += 1
                    if atom_name.startswith("C"):
                        carbon_atom_count += 1
                    if residue_id not in residues_seen:
                        residues_seen.add(residue_id)
                        residue_count += 1
                        total_hydro += self.hydro.get(residue_name, self.hydro["UNK"])
                        if residue_name in aromatic_residues:
                            aromatic_count += 1
                    if residue_name == "TYR" and atom_name == "OH":
                        tyr_oh_count += 1
        avg_hydro = total_hydro / residue_count if residue_count > 0 else 0.0
        aromatic_freq = aromatic_count / residue_count if residue_count > 0 else 0.0
        tyr_oh_freq = tyr_oh_count / residue_count if residue_count > 0 else 0.0
        carbon_freq = (
            carbon_atom_count / total_atom_count if total_atom_count > 0 else 0.0
        )
        return avg_hydro, aromatic_freq, tyr_oh_freq, carbon_freq

    @staticmethod
    def score_cavity(
        volume: float,
        surface: float,
        diameter: float,
        smallest_size: float,
        hydrophobicity_score: float,
        aromatic_freq: float,
        otyr_freq: float,
        radius: float,
        cylinder_radius: float,
    ) -> float:
        """
        Calculates the cavity score as a simple mean of selected descriptors.

        Parameters
        ----------
            volume : float
                Cavity volume.
            surface : float
                Cavity surface area.
            diameter : float
                Cavity diameter.
            smallest_size : float
                Smallest size of the cavity.
            hydrophobicity_score : float
                Average hydrophobicity score.
            aromatic_freq : float
                Aromatic residue frequency.
            otyr_freq : float
                Tyrosine OH frequency.
            radius : float
                Cavity radius.
            cylinder_radius : float
                Minimal enclosing cylinder radius.

        Return
        ------
            float
                Average score of the selected descriptors.
        """
        values = [
            volume,
            surface,
            diameter,
            smallest_size,
            hydrophobicity_score,
            aromatic_freq,
            otyr_freq,
            radius,
            cylinder_radius,
        ]
        return sum(values) / len(values) if values else 0.0

    def calculate_geometric_center(
        self, pdb_lines: List[str]
    ) -> Tuple[float, float, float]:
        center = pdb_utils.calculate_geometric_center(pdb_lines)
        return tuple(center)

    def save_cavity_scores(
        self, cavity_path: str, scores: List[Dict[str, Any]]
    ) -> None:
        """
        Saves the calculated cavity scores for different programs into a CSV file.

        Parameters
        ----------
            cavity_path : str
                Path to the directory containing the cavity data.
            scores : list
                List of dictionaries with cavity data and computed descriptors.
        """
        # Define output file path
        protein_code = os.path.basename(os.path.dirname(cavity_path))
        output_file = os.path.join(cavity_path, f"{protein_code}_cavity_scores.csv")
        # Define CSV field order, now including geometric center coordinates
        fieldnames = [
            "Protein_Code",
            "Cavity_ID",
            "Program",
            "Grid_File",
            "Residue_File",
            "Grid_Path",
            "Residue_Path",
            "Volume",
            "Radius",
            "Surface",
            "Diameter",
            "Smallest_Size",
            "Cylinder_Radius",
            "Hydrophobicity_Score",
            "Aromatic_Freq",
            "Otyr_Freq",
            "Carbon_Freq",
            "Geom_Center_X",
            "Geom_Center_Y",
            "Geom_Center_Z",
            "Geom_Center",
            "Score",
        ]
        # Write to CSV
        with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=",")
            writer.writeheader()
            # Sort cavities by descending score
            for entry in sorted(scores, key=lambda item: item["score"], reverse=True):
                writer.writerow(
                    {
                        "Protein_Code": entry["protein_code"],
                        "Cavity_ID": entry["cavity_id"],
                        "Program": entry["program"],
                        "Grid_File": entry["grid_file"],
                        "Residue_File": entry["residue_file"],
                        "Grid_Path": entry["grid_path"],
                        "Residue_Path": entry["residue_path"],
                        "Volume": f"{entry['volume']:.2f}",
                        "Radius": f"{entry['radius']:.2f}",
                        "Surface": f"{entry['surface']:.2f}",
                        "Diameter": f"{entry['diameter']:.2f}",
                        "Smallest_Size": f"{entry['smallest_size']:.2f}",
                        "Cylinder_Radius": f"{entry['cylinder_radius']:.2f}",
                        "Hydrophobicity_Score": f"{entry['hydrophobicity_score']:.2f}",
                        "Aromatic_Freq": f"{entry['aromatic_freq']:.2f}",
                        "Otyr_Freq": f"{entry['tyr_oh_freq']:.2f}",
                        "Carbon_Freq": f"{entry['carbon_freq']:.2f}",
                        "Geom_Center_X": f"{entry['geom_center_x']:.3f}",
                        "Geom_Center_Y": f"{entry['geom_center_y']:.3f}",
                        "Geom_Center_Z": f"{entry['geom_center_z']:.3f}",
                        "Geom_Center": entry[
                            "geom_center"
                        ],  # string formatada "x, y, z"
                        "Score": f"{entry['score']:.2f}",
                    }
                )
        self.logger.info(f"Cavity scores saved to {output_file}")

    def score_routine(
        self,
        cavity_path: str,
        program_cavities: Dict[str, List[Dict[str, Any]]],
        cavity_ids: Dict[Optional[str], int],
    ) -> None:
        """
        Main scoring routine that processes cavities from different programs,
        computes descriptors, and aggregates all scores.

        Parameters
        ----------
            cavity_path : str
                Path to the cavity directory.
            program_cavities : dict
                Dictionary mapping program names to lists of cavity dictionaries.
            cavity_ids : dict
                Dictionary mapping grid file paths to cavity identifiers.
        """
        protein_code = os.path.basename(os.path.dirname(cavity_path))
        scores = []
        total_programs = len(program_cavities)
        total_cavities = sum(len(c) for c in program_cavities.values())
        self.logger.debug(
            "score_routine: protein=%s, programs=%d, cavities=%d",
            protein_code,
            total_programs,
            total_cavities,
        )
        for program, cavities in program_cavities.items():
            self.logger.debug(
                "  Processing program: %s (%d cavities)", program, len(cavities)
            )
            for cavity in cavities:
                grid_file_path = cavity["grid_path"]
                residue_file_path = cavity["residue_path"]
                # Skip cavities without a grid file (e.g., Fukui-only cavities)
                if grid_file_path is None:
                    self.logger.warning(
                        f"Skipping cavity without grid: {residue_file_path}"
                    )
                    continue
                cavity_id = cavity_ids.get(grid_file_path, None)
                # Read grid lines and calculate volume
                with open(grid_file_path, "r") as f:
                    grid_lines = f.readlines()
                volume, surface, radius, diameter, smallest_size, cylinder_radius = (
                    self.calculate_geometric_descriptors(grid_lines, grid_file_path)
                )
                # Compute geometric center
                geom_center = self.calculate_geometric_center(grid_lines)
                geom_center_str = (
                    f"{geom_center[0]:.3f}, {geom_center[1]:.3f}, {geom_center[2]:.3f}"
                )
                # Calculate other scores (e.g., hydrophobicity, polarity)
                hydrophobicity_score, aromatic_freq, tyr_oh_freq, carbon_freq = (
                    self.calculate_residue_descriptors(residue_file_path)
                )
                score = self.score_cavity(
                    volume,
                    surface,
                    diameter,
                    smallest_size,
                    hydrophobicity_score,
                    aromatic_freq,
                    tyr_oh_freq,
                    radius,
                    cylinder_radius,
                )
                # Append the results to the scores list
                scores.append(
                    {
                        "protein_code": protein_code,
                        "cavity_id": cavity_id,
                        "program": program,
                        "grid_file": os.path.basename(grid_file_path),
                        "residue_file": os.path.basename(residue_file_path),
                        "grid_path": grid_file_path,
                        "residue_path": residue_file_path,
                        "volume": volume,
                        "radius": radius,
                        "surface": surface,
                        "diameter": diameter,
                        "smallest_size": smallest_size,
                        "cylinder_radius": cylinder_radius,
                        "hydrophobicity_score": hydrophobicity_score,
                        "aromatic_freq": aromatic_freq,
                        "tyr_oh_freq": tyr_oh_freq,
                        "carbon_freq": carbon_freq,
                        "geom_center_x": geom_center[0],
                        "geom_center_y": geom_center[1],
                        "geom_center_z": geom_center[2],
                        "geom_center": geom_center_str,
                        "score": score,
                    }
                )
        self.save_cavity_scores(cavity_path, scores)
