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
from hocuspockets.initialize import *

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
        Initialize a CatalyticCav object.

        The original PDB residue ordering is preserved to ensure
        correct mapping of Fukui and Fukushima values stored in
        the B-factor field of external PDB files that do not
        contain chain information.

        Parameters
        ----------
        cavities_path : str
            Path to the directory containing cavity output files.

        pdb_path : str
            Path to the original input PDB structure containing
            complete residue and chain information.

        fukui_file : str
            PDB file containing Fukui values encoded in the
            B-factor column.

        fukushima_file : str
            PDB file containing Fukushima values encoded in the
            B-factor column.

        verbose_handler : VerboseHandler
            Logger used to report execution information and warnings.

        Notes
        -----
        The residue order from the original structure is used as the
        reference for matching atoms and residues between the input
        PDB and the Fukui/Fukushima descriptor files.
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
            Parameters
            ----------
                None
            Return
            ------
                residue_map : dict
                    Maps index (1-based) → 'chain:resSeqresName'
            Notes
            -----
                This mapping is used to relate the sequential order of atoms in the Fukui/Fukushima
                PDB files (which lack chain information) back to the original residues with chain IDs.
            
            Example
            --------
                residue_map = self.extract_residue_index_mapping()
                # residue_map[1] might be 'A:42LYS', residue_map[2] might be 'A:43GLU', etc.
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
        Parses a PDB file to extract B-factor values for individual atoms (not residues),
        mapping them to the corresponding residue and chain using the reference PDB order.
        Parameters
        ----------
            pdb_file : str
                - PDB file with descriptor values encoded in B-factor.
            descriptor_name : str
                - Name for logging ("Fukui", "Fukushima", etc).
        Return
        ------
            atom_values : dict
                - Maps 'chain:resSeqresName:atomName' → B-factor value
                e.g. "A:42LYS:CA" → 0.054   
        Notes
        -----
            The method assumes that the order of atoms in the input PDB file matches
            the order of residues in the original PDB, allowing for correct mapping
            of B-factor values to residues with chain information.

        Example
        -------
            atom_values = self.load_bfactor_values("fukui.pdb", descriptor_name="Fukui")
        """
        atom_values = {}
        current_res_index = None
        res_id = None

        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue

                try:
                    res_index = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    bfactor = float(line[60:66].strip())
                except ValueError:
                    continue

                # Atualiza resíduo atual (a partir do PDB original)
                if res_index != current_res_index:
                    res_id = self.index_to_residue_id.get(res_index)
                    current_res_index = res_index

                if not res_id:
                    continue

                # Cria chave única: cadeia + resíduo + átomo
                atom_key = f"{res_id}:{atom_name}"
                atom_values[atom_key] = bfactor

        self.verbose_handler.output(
            f"[INFO] Parsed {len(atom_values)} atom-level {descriptor_name} values from: {os.path.basename(pdb_file)}"
        )
        return atom_values


    def get_residues_from_pdb(self, pdb_path):
        """
        Extracts a set of unique residues from a PDB file, using a composite ID that includes:
        chain ID, residue sequence number, and residue name.
        Parameters
        ----------
            pdb_path : str
                Path to the input PDB file.
        Return
        ------
            residues : str
                A set of residue IDs in the format 'chain:resSeqresName' (e.g., 'A:42LYS').
        Notes
        -----
            This method is useful for quickly identifying all residues present in a PDB file,
            regardless of atom-level details. It ignores duplicate residues and only considers
            lines starting with "ATOM".
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
        Computes the average Fukui and Fukushima descriptor values for each predicted cavity,
        using per-atom values rather than per-residue averages.
        Parameters
        ----------
            None
        Return
        ------
            results_fukushima : dict
                Maps (program, pocket_file) → average Fukushima value
            results_fukui : dict
                Maps (program, pocket_file) → average Fukui value
            results_fukushima_no_norm : dict
                Maps (program, pocket_file) → sum of Fukushima values (no normalization)
            results_fukushima_sum_norm : dict
                Maps (program, pocket_file) → sum of Fukushima values normalized by number of active atoms
        Notes
        -----
        Each cavity file (residue*.pdb) is parsed to extract atoms, whose chain, residue, and atom names
        are used to look up the corresponding Fukui/Fukushima reactivity values.
        
        Example:
            results_fukushima, results_fukui, results_fukushima_no_norm, results_fukushima_sum_norm = self.compute_average()
        """
        cavity_programs = ['fpocket', 'ligsite', 'surfnet', 'pocketfinder']
        results_fukushima = {}
        results_fukui = {}
        results_fukushima_no_norm = {}
        results_fukushima_sum_norm = {}

        for program in cavity_programs:
            program_dir = os.path.join(self.cavities_path, program)
            if not os.path.exists(program_dir):
                self.verbose_handler.output(f"[WARNING] Skipping missing cavity directory: {program_dir}")
                continue

            for pocket_file in os.listdir(program_dir):
                if not (pocket_file.startswith("residue") and pocket_file.endswith(".pdb")):
                    continue

                pocket_path = os.path.join(program_dir, pocket_file)
                total_fukui = 0.0
                total_fukushima = 0.0
                atom_count = 0
                active_atom_count = 0

                with open(pocket_path, 'r') as pdb:
                    for line in pdb:
                        if not line.startswith("ATOM"):
                            continue

                        # Extrai informações do átomo da cavidade
                        chain_id = line[21].strip()
                        res_seq = line[22:26].strip()
                        res_name = line[17:20].strip()
                        atom_name = line[12:16].strip()
                        atom_id = f"{chain_id}:{res_seq}{res_name}:{atom_name}"

                        # Recupera valor Fukui
                        if atom_id in self.fukui_dict:
                            total_fukui += self.fukui_dict[atom_id]

                        # Recupera valor Fukushima
                        if atom_id in self.fukushima_dict:
                            val_fuku = self.fukushima_dict[atom_id]
                            total_fukushima += val_fuku
                            if val_fuku != 0.0:
                                active_atom_count += 1

                        atom_count += 1

                # Calcula médias e normalizações
                avg_fukui = total_fukui / atom_count if atom_count > 0 else None
                avg_fukushima = total_fukushima / atom_count if atom_count > 0 else None
                fukushima_no_norm = total_fukushima
                fukushima_sum_norm = total_fukushima / active_atom_count if active_atom_count > 0 else None

                results_fukui[(program, pocket_file)] = avg_fukui
                results_fukushima[(program, pocket_file)] = avg_fukushima
                results_fukushima_no_norm[(program, pocket_file)] = fukushima_no_norm
                results_fukushima_sum_norm[(program, pocket_file)] = fukushima_sum_norm

                self.verbose_handler.output(
                    f"[INFO] {program}/{pocket_file} → Atoms: {atom_count}, Active atoms: {active_atom_count}, "
                    f"Avg_Fukushima={avg_fukushima}, Avg_Fukui={avg_fukui}, "
                    f"Sum_Fukushima={fukushima_no_norm}, Norm_Fukushima={fukushima_sum_norm}"
                )

        return results_fukushima, results_fukui, results_fukushima_no_norm, results_fukushima_sum_norm


    def process_fukui_cav(self):
        """
        Identifies catalytic residues based on the top 5% Fukui values (from B-factor column)
        in the Fukui PDB file (which lacks chain information), but maps back to residues with
        chain using the original input PDB.
        Parameters
        ----------
            None

        Return
        ------
            None

        Notes
        -----   
        The method reads the Fukui PDB file, computes the cutoff for the top 5% of B-factor values,
        and then extracts all atoms of residues that have at least one atom above this cutoff.
        The resulting atoms are saved to a new PDB file in the 'fukui' subdirectory of the cavities path.
        
        Example
        ------- 
            self.process_fukui_cav()
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