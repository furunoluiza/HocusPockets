import os
import sys
import csv
import shutil
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# ========== Seleção de resíduos ==========
class ResidueSelect(Select):
    def __init__(self, selected_residues):
        self.selected_residues = selected_residues

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        res_id = residue.id
        return (chain_id, res_id) in self.selected_residues

# ========== Funções utilitárias ==========
def safe_copy(src, dst_dir, new_name=None):
    """Copia um arquivo para o diretório destino, opcionalmente renomeando."""
    if not src or not os.path.isfile(src):
        return None
    os.makedirs(dst_dir, exist_ok=True)
    dst = os.path.join(dst_dir, new_name if new_name else os.path.basename(src))
    shutil.copy2(src, dst)
    return dst

def extract_residues_from_pdb(cavity_file):
    residues = set()
    with open(cavity_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_number = line[22:26].strip()
                residue_id = f"{chain_id}:{res_number}{res_name}"
                residues.add(residue_id)
    return residues

def compute_similarity(cavity_path, pocket_path):
    """Similaridade entre cavidade prevista e pocket de referência."""
    if not (os.path.isfile(cavity_path) and os.path.isfile(pocket_path)):
        return None

    cav_residues = extract_residues_from_pdb(cavity_path)
    pocket_residues = extract_residues_from_pdb(pocket_path)

    if not cav_residues or not pocket_residues:
        return 0.0

    overlap = cav_residues & pocket_residues
    denom = min(len(cav_residues), len(pocket_residues))
    return len(overlap) / denom

def extract_nearby_residues_from_ligand(pdb_file, mol2_file, output_dir, name, cutoff_factor=2.0):
    """
    Extrai resíduos próximos ao ligante (dentro de raio = Rg*cutoff_factor).
    """
    mol = Chem.MolFromMol2File(mol2_file, sanitize=False)
    if mol is None:
        print(f"[ERROR] Falha ao carregar ligante: {mol2_file}")
        return None
    rg = rdMolDescriptors.CalcRadiusOfGyration(mol)
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    center = coords.mean(axis=0)
    cutoff = cutoff_factor * rg

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    selected_residues = set()
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if np.linalg.norm(atom.coord - center) <= cutoff:
                        selected_residues.add((chain.id, residue.id))
                        break

    os.makedirs(output_dir, exist_ok=True)
    pocket_path = os.path.join(output_dir, f"{name}_pocket.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(pocket_path, ResidueSelect(selected_residues))
    return pocket_path

def build_all_cavities_summary(struct_dir, pocket_gen_dir, output_csv, threshold=0.5):
    """Gera arquivo consolidado com todas as cavidades e flag IsReference."""
    rows_out = []
    for name in os.listdir(struct_dir):
        print(name)
        protein_dir = os.path.join(struct_dir, name)
        print(protein_dir)
        if not os.path.isdir(protein_dir):
            continue
        cav_dir = os.path.join(protein_dir, 'cavities')
        score_file = os.path.join(cav_dir, f"{name}_cavity_scores.csv")
        print(score_file)
        pocket_path = os.path.join(pocket_gen_dir, f"{name}_pocket.pdb")

        if not (os.path.isfile(score_file) and os.path.isfile(pocket_path)):
            continue

        with open(score_file, 'r') as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                cav_path = row.get("Residue_Path")
                is_ref = False
                print(cav_path)
                print(pocket_path)
                if cav_path and os.path.isfile(cav_path):
                    sim = compute_similarity(cav_path, pocket_path)
                    print(cav_path)
                    print(sim)
                    if sim is not None and sim >= threshold:
                        is_ref = True

                row_out = dict(row)
                row_out["Protein"] = name
                row_out["RowID"] = f"{name}_{row['Residue_File']}"
                row_out["Similarity"] = sim if sim is not None else 0.0
                row_out["IsReference"] = is_ref
                rows_out.append(row_out)

    if rows_out:
        fieldnames = list(rows_out[0].keys())
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows_out)
        print(f"[INFO] Summary CSV gerado: {output_csv}")
    else:
        print("[WARN] Nenhuma cavidade encontrada para gerar summary CSV.")

# ========== MAIN ==========
def main(pdbs_dir, struct_dir, ligants_dir):
    pocket_gen_dir = "pockets"
    os.makedirs(pocket_gen_dir, exist_ok=True)

    bundles_parent = "bundle"
    os.makedirs(bundles_parent, exist_ok=True)

    for name in os.listdir(struct_dir):
        protein_dir = os.path.join(struct_dir, name)
        if not os.path.isdir(protein_dir):
            continue

        base_id = name.split("_")[0]
        cav_dir = os.path.join(protein_dir, 'cavities')
        score_file = os.path.join(cav_dir, f"{name}_cavity_scores.csv")
        pdb_path = os.path.join(pdbs_dir, f"{base_id}.pdb")
        ligand_path = os.path.join(ligants_dir, f"{base_id}_ligand.mol2")

        if not (os.path.isfile(score_file) and os.path.isfile(pdb_path) and os.path.isfile(ligand_path)):
            print(f"[SKIP] Arquivos ausentes para {name}")
            continue

        pocket_path = extract_nearby_residues_from_ligand(pdb_path, ligand_path, pocket_gen_dir, name)
        if pocket_path is None:
            continue

        # Criar bundle
        """bundle_dir = os.path.join(bundles_parent, name)
        os.makedirs(bundle_dir, exist_ok=True)
        safe_copy(pdb_path, bundle_dir, new_name=f"{name}_protein_original.pdb")
        safe_copy(ligand_path, bundle_dir, new_name=f"{name}_ligand.mol2")
        safe_copy(pocket_path, bundle_dir, new_name=f"{name}_pocket_generated.pdb")
        safe_copy(score_file, bundle_dir, new_name=f"{name}_cavity_scores.csv")"""

    # Gera CSV consolidado
    
    build_all_cavities_summary(struct_dir, pocket_gen_dir, "all_cavities_summary.csv", threshold=0.5)

# ========== ENTRY POINT ==========
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python generate_summary.py <pdbs_dir> <struct_dir> <ligants_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])