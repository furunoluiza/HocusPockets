#!/usr/bin/env python3

# Description
###############################################################################
'''
    This module defines the CatalyticCav class, which is responsible for analyzing
    cavity residues based on quantum reactivity descriptors calculated by PRIMoRDiA,
    specifically Fukui and Fukushima indices.

    Key functionalities:
    - Maps residues from the .rslrd descriptor file to full identifiers including chain ID,
    using the original input PDB structure as reference.
    - Computes average Fukui and Fukushima values for each predicted cavity from tools
    such as Fpocket and ConCavity.
    - Identifies catalytic residues by selecting those with atoms among the top 5% highest
    Fukui potentials, and outputs a PDB file with their full atomic representation.
    - Supports integration with other modules for comparative scoring and CSV reporting.

    This module is a core component of the HocusPockets workflow for identifying and prioritizing
    potentially druggable and catalytically relevant pockets in protein structures.
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

class CatalyticCav:
    def __init__(self, cavities_path, pdb_path, fukui_file, fukushima_file, verbose_handler):
        """
        Initializes the CatalyticCav object using original PDB residue order to match Fukui and Fukushima
        values (stored in the B-factor field) from corresponding PDBs that lack chain information.
        Input:
            cavities_path    [str] - Path to the directory containing cavity output files.
            pdb_path         [str] - Path to the original input PDB with full residue and chain info.
            fukui_file       [str] - PDB file where Fukui values are encoded in the B-factor field.
            fukushima_file   [str] - PDB file where Fukushima values are encoded in the B-factor field.
            verbose_handler  [VerboseHandler] - Logger for informative output.
        """
        self.cavities_path = cavities_path
        self.pdb_path = pdb_path
        self.fukui_file = fukui_file
        self.fukushima_file = fukushima_file
        self.verbose_handler = verbose_handler
        # Step 1: build ordered residue list from original PDB
        self.index_to_residue_id = self.extract_residue_index_mapping()
        # Step 2: parse values from fukui and fukushima PDBs
        self.fukui_dict = self.load_bfactor_values(self.fukui_file, descriptor_name="Fukui")
        self.fukushima_dict = self.load_bfactor_values(self.fukushima_file, descriptor_name="Fukushima")

    def extract_residue_index_mapping(self):
            """
            Extracts the list of residues from the original PDB, ordered by first appearance,
            and maps sequential index (starting at 1) → 'chain:resSeqresName'.
            Input:
                Nothing.
            Return:
                [dict]: Maps index (1-based) → 'chain:resSeqresName'
            """
            residue_map = {} # Final mapping from sequential index to full residue ID
            seen = set() # To track residues already encountered (avoid duplicates)
            index = 1 # Start indexing at 1 to match the .rslrd file
            with open(self.pdb_path, 'r') as pdb:
                for line in pdb:
                    if line.startswith("ATOM"):
                        chain_id = line[21].strip()
                        res_seq = line[22:26].strip()
                        res_name = line[17:20].strip()
                        res_id = f"{chain_id}:{res_seq}{res_name}" # Construct full residue ID

                        if res_id not in seen:
                            seen.add(res_id)
                            residue_map[index] = res_id # Map current index to the residue ID
                            index += 1
            self.verbose_handler.output(f"[INFO] Extracted {len(residue_map)} residues with chain IDs from PDB.")
            return residue_map

    def load_bfactor_values(self, pdb_file, descriptor_name="Descriptor"):
        """
        Parses a PDB file to extract average B-factor per residue, assuming residues are ordered numerically
        from 1 upward (no chain info), and maps these values to proper residue IDs using the index-to-ID map.
        Input:
            pdb_file [str]        - PDB file with descriptor values encoded in B-factor.
            descriptor_name [str] - Name for logging ("Fukui", "Fukushima", etc).
        Return:
            [dict] residue_id → average B-factor value
        """
        residue_values = {}
        current_residue_index = None
        current_values = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        res_index = int(line[22:26].strip())
                        bfactor = float(line[60:66].strip())
                    except ValueError:
                        continue
                    # New residue encountered
                    if current_residue_index != res_index:
                        if current_residue_index is not None and current_values:
                            residue_id = self.index_to_residue_id.get(current_residue_index)
                            if residue_id:
                                avg = sum(current_values) / len(current_values)
                                residue_values[residue_id] = avg
                        current_residue_index = res_index
                        current_values = [bfactor]
                    else:
                        current_values.append(bfactor)
            # Final residue
            if current_residue_index is not None and current_values:
                residue_id = self.index_to_residue_id.get(current_residue_index)
                if residue_id:
                    avg = sum(current_values) / len(current_values)
                    residue_values[residue_id] = avg
        self.verbose_handler.output(f"[INFO] Parsed {len(residue_values)} residue-level {descriptor_name} values from: {os.path.basename(pdb_file)}")
        return residue_values

    def get_residues_from_pdb(self, pdb_path):
        """
        Extracts a set of unique residues from a PDB file, using a composite ID that includes:
        chain ID, residue sequence number, and residue name.
        Input:
            pdb_path [str]: Path to the input PDB file.
        Return:
            residues[set]: A set of residue IDs in the format 'chain:resSeqresName' (e.g., 'A:42LYS').
        """
        residues = set()
        with open(pdb_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    res_name = line[17:20].strip()
                    res_id = f"{chain_id}:{res_seq}{res_name}" # Construct residue ID with chain:resSeqresName format
                    residues.add(res_id)
        return residues

    def compute_average(self):
        """
        Computes the average Fukushima and Fukui descriptor values for each predicted cavity,
        normalized by the number of atoms rather than residues.
        Each pocket file (residue*.pdb) is parsed to extract atoms and identify their corresponding
        residues. Descriptor values are then accumulated and averaged across atoms that belong to
        residues with known descriptors.
        Input:
            None
        Return:
            Tuple:
                - results_fukushima         [dict]: Maps (method, pocket_file) to average Fukushima value per atom.
                - results_fukui             [dict]: Maps (method, pocket_file) to average Fukui value per atom.
                - results_fukushima_no_norm [dict]: (program, pocket_file) -> total (non-normalized) Fukushima sum.
                - results_fukushima_sum_norm[dict]: (program, pocket_file) -> normalized Fukushima sum.
        """
        cavity_programs = ['fpocket', 'ligsite', 'surfnet', 'pocketfinder']
        results_fukushima = {}
        results_fukui = {}
        results_fukushima_no_norm = {}
        for program in cavity_programs:
            program_dir = os.path.join(self.cavities_path, program)
            if not os.path.exists(program_dir):
                self.verbose_handler.output(f"[WARNING] Skipping missing cavity directory: {program_dir}")
                continue
            for pocket_file in os.listdir(program_dir):
                if pocket_file.startswith("residue") and pocket_file.endswith(".pdb"):
                    pocket_path = os.path.join(program_dir, pocket_file)
                    total_fukui = 0.0
                    total_fukushima = 0.0
                    atom_count = 0
                    nonzero_atoms = 0
                    with open(pocket_path, 'r') as pdb:
                        for line in pdb:
                            if line.startswith("ATOM"):
                                chain_id = line[21].strip()
                                res_seq = line[22:26].strip()
                                res_name = line[17:20].strip()
                                residue_id = f"{chain_id}:{res_seq}{res_name}"
                                if residue_id in self.fukui_dict:
                                    total_fukui += self.fukui_dict[residue_id]
                                if residue_id in self.fukushima_dict:
                                    val_fuku = self.fukushima_dict[residue_id]
                                    total_fukushima += val_fuku
                                    if val_fuku != 0.0:
                                        nonzero_atoms += 1
                                atom_count += 1
                    avg_fukui = total_fukui / atom_count if atom_count else None
                    avg_fukushima = total_fukushima / atom_count if atom_count else None
                    fukushima_no_norm = total_fukushima
                    fukushima_sum_norm = (total_fukushima / nonzero_atoms if nonzero_atoms > 0 else None)
                    results_fukui[(program, pocket_file)] = avg_fukui
                    results_fukushima[(program, pocket_file)] = avg_fukushima
                    results_fukushima_no_norm[(program, pocket_file)] = fukushima_no_norm
                    results_fukushima_sum_norm[(program, pocket_file)] = fukushima_sum_norm
                    self.verbose_handler.output(
                    f"[INFO] {program}/{pocket_file} → Atoms: {atom_count}, "
                        f"Avg. Fukushima: {avg_fukushima}, Avg. Fukui: {avg_fukui}, "
                        f"Fukushima_No_Norm: {fukushima_no_norm}, "
                        f"Fukushima_Sum_Norm: {fukushima_sum_norm}"
                    )

        return results_fukushima, results_fukui, results_fukushima_no_norm, results_fukushima_sum_norm


    def process_fukui_cav(self):
        """
        Identifies catalytic residues based on the top 5% Fukui values (from B-factor column)
        in the Fukui PDB file (which lacks chain information), but maps back to residues with
        chain using the original input PDB.
        The output is a new PDB file containing all atoms of the residues considered catalytic.
        Input:
            Nothing.
        Return:
            Nothing.
        """
        if not os.path.exists(self.pdb_path):
            self.verbose_handler.output(f"[ERROR] Original PDB file not found: {self.pdb_path}")
            return
        fukui_bfactors = []
        # Step 1: Extract B-factors from the Fukui PDB (assumed to be in same atom order as original PDB)
        with open(self.fukui_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    try:
                        bfactor = float(line[60:66].strip())
                        fukui_bfactors.append(bfactor)
                    except ValueError:
                        fukui_bfactors.append(0.0)
        if not fukui_bfactors:
            self.verbose_handler.output("[WARNING] No B-factor values found in Fukui file.")
            return
        # Step 2: Compute cutoff
        sorted_bf = sorted(fukui_bfactors)
        cutoff_index = int(0.95 * len(sorted_bf))
        cutoff_value = sorted_bf[cutoff_index]
        self.verbose_handler.output(f"[INFO] Fukui cutoff: {cutoff_value:.4f} (top 5%)")
        # Step 3: Read original PDB and collect catalytic residues based on cutoff
        selected_residues = set()
        pdb_atoms = []
        with open(self.pdb_path, 'r') as file:
            for i, line in enumerate(file):
                if line.startswith("ATOM") and i < len(fukui_bfactors):
                    chain_id = line[21].strip()
                    res_seq = line[22:26].strip()
                    res_name = line[17:20].strip()
                    res_id = f"{chain_id}:{res_seq}{res_name}"
                    pdb_atoms.append((line, res_id, fukui_bfactors[i]))
                    if fukui_bfactors[i] >= cutoff_value:
                        selected_residues.add(res_id)
        # Step 4: Extract all atoms of selected residues
        selected_lines = [line for line, res_id, _ in pdb_atoms if res_id in selected_residues]
        # Step 5: Save to output file
        fukui_dir = os.path.join(self.cavities_path, 'fukui')
        os.makedirs(fukui_dir, exist_ok=True)
        output_file = os.path.join(fukui_dir, 'residue_0_fukui.pdb')
        with open(output_file, 'w') as out:
            out.writelines(selected_lines)
        self.verbose_handler.output(f"[SUCCESS] Catalytic residues saved to: {output_file}")