import sys
import os
import pandas as pd

def extract_residue_index_mapping(pdb_path):
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
    with open(pdb_path, 'r') as pdb:
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
    return residue_map

def load_bfactor_values(residue_map, pdb_file, descriptor_name="Descriptor"):
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
                        residue_id = residue_map.get(current_residue_index)
                        if residue_id:
                            avg = sum(current_values) / len(current_values)
                            residue_values[residue_id] = avg
                    current_residue_index = res_index
                    current_values = [bfactor]
                else:
                    current_values.append(bfactor)
        # Final residue
        if current_residue_index is not None and current_values:
            residue_id = residue_map.get(current_residue_index)
            if residue_id:
                avg = sum(current_values) / len(current_values)
                residue_values[residue_id] = avg
    return residue_values

def compute_values(fukui_dict, fukushima_dict, pocket_file):
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
    total_fukui = 0.0
    total_fukushima = 0.0
    atom_count = 0
    nonzero_atoms = 0
    with open(pocket_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                res_seq = line[22:26].strip()
                res_name = line[17:20].strip()
                residue_id = f"{chain_id}:{res_seq}{res_name}"
                if residue_id in fukui_dict:
                    total_fukui += fukui_dict[residue_id]
                if residue_id in fukushima_dict:
                    val_fuku = fukushima_dict[residue_id]
                    total_fukushima += val_fuku
                    if val_fuku != 0.0:
                        nonzero_atoms += 1
                atom_count += 1
    avg_fukui = total_fukui / atom_count if atom_count else None
    avg_fukushima = total_fukushima / atom_count if atom_count else None
    fukushima_no_norm = total_fukushima
    fukushima_sum_norm = (total_fukushima / nonzero_atoms if nonzero_atoms > 0 else None)
    return avg_fukui, avg_fukushima, fukushima_no_norm, fukushima_sum_norm

def compute_average(pdb_path, pocket_file, fukui_file, fukushima_file):
    residue_map = extract_residue_index_mapping(pdb_path)
    fukui_dict = load_bfactor_values(residue_map, fukui_file, descriptor_name="Fukui")
    fukushima_dict = load_bfactor_values(residue_map, fukushima_file, descriptor_name="Fukushima")
    avg_fukui, avg_fukushima, fukushima_no_norm, fukushima_sum_norm = compute_values(fukui_dict, fukushima_dict, pocket_file)
    return(avg_fukui, avg_fukushima, fukushima_no_norm, fukushima_sum_norm)

def main(pdb_dir, pockets_dir, output_hocus, pdbbind_csv):
    df = pd.read_csv(pdbbind_csv)
    if "ID" not in df.columns:
        print("[ERROR] CSV must contain a column named 'ID'.")
        sys.exit(1)
    new_cols = ["Avg_Fukui", "Avg_Fukushima", "Fukushima_NoNorm", "Fukushima_SumNorm"]
    for col in new_cols:
        if col not in df.columns:
            df[col] = None

    for file in os.listdir(pdb_dir):
        if not pdb_path.endswith(".pdb"):
            continue
        pdb_path = os.path.join(pdb_dir, file)
        protein_name = os.path.splitext(file)[0]
        protein_output_dir = os.path.join(output_hocus, protein_name)
        fukui_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukui_pot_zero.pdb')
        fukushima_file = os.path.join(protein_output_dir, 'reactivity', 'primordia', protein_name + '_PDB_RD', protein_name + '_fukushima.pdb')
        pocket_file = os.path.join(pockets_dir, f"{protein_name}_pocket.pdb")
        if not all(os.path.exists(p) for p in [fukui_file, fukushima_file, pocket_file]):
            print(f"[Warning] Missing files for {protein_name}, skipping.")
            continue

        avg_fukui, avg_fukushima, fukushima_no_norm, fukushima_sum_norm  = compute_average(pdb_path, pocket_file, fukui_file, fukushima_file)
        mask = df["ID"].astype(str).str.lower() == protein_name.lower()
        if mask.any():
            df.loc[mask, "Avg_Fukui"] = avg_fukui
            df.loc[mask, "Avg_Fukushima"] = avg_fukushima
            df.loc[mask, "Fukushima_NoNorm"] = fuku_no_norm
            df.loc[mask, "Fukushima_SumNorm"] = fuku_sum_norm
            print(f"[OK] Updated {protein_name}")
        else:
            print(f"[Info] ID {protein_name} not found in CSV, skipping.")
        
        df.to_csv(pdbbind_csv, index=False)
        print(f"\n[Done] Updated CSV saved at: {pdbbind_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python3 drug_pdbbind_csv.py <pdb_dir> <pockets_dir> <output_hocus> <pdbbind_csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])