#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import pandas as pd
import numpy as np

# ============================================================
# CLASSIFICAÇÃO BASEADA EM KD/KI
# ============================================================
def classificar_drogabilidade(kd_value):
    """Classifica o alvo como Drogável (d), Borderline (b) ou Não Drogável (n) com base em KD/KI (µM)."""
    if pd.isna(kd_value):
        return np.nan
    try:
        kd_value = float(kd_value)
    except ValueError:
        return np.nan
    # Converte nM → µM, se necessário
    if kd_value > 10000:
        kd_value /= 1000.0
    if kd_value < 0.3:
        return "d"  # Drogável
    elif kd_value <= 1.0:
        return "b"  # Borderline
    else:
        return "n"  # Não Drogável


# ============================================================
# MAPEAMENTO DE RESÍDUOS DO PDB ORIGINAL
# ============================================================
def extract_residue_index_mapping(pdb_path):
    """Gera mapeamento indexado de resíduos (1 → cadeia:resSeqresName) com base no PDB original."""
    residue_map = {}
    seen = set()
    index = 1
    with open(pdb_path, 'r') as pdb:
        for line in pdb:
            if line.startswith("ATOM"):
                chain_id = line[21].strip()
                res_seq = line[22:26].strip()
                res_name = line[17:20].strip()
                res_id = f"{chain_id}:{res_seq}{res_name}"
                if res_id not in seen:
                    seen.add(res_id)
                    residue_map[index] = res_id
                    index += 1
    return residue_map


# ============================================================
# LEITURA DOS VALORES B-FACTOR POR ÁTOMO (FUKUI / FUKUSHIMA)
# ============================================================
def load_bfactor_values(residue_map, pdb_file, descriptor_name="Descriptor"):
    """
    Lê os valores de B-factor por átomo, associando-os ao resíduo e cadeia corretos
    com base no índice sequencial do PDB original.
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

            if res_index != current_res_index:
                res_id = residue_map.get(res_index)
                current_res_index = res_index
            if not res_id:
                continue

            atom_key = f"{res_id}:{atom_name}"
            atom_values[atom_key] = bfactor

    print(f"[INFO] Parsed {len(atom_values)} atom-level {descriptor_name} values from: {os.path.basename(pdb_file)}")
    return atom_values


# ============================================================
# CÁLCULO DE MÉDIAS POR CAVIDADE
# ============================================================
def compute_values(fukui_dict, fukushima_dict, pocket_file):
    """
    Calcula médias e somas Fukushima/Fukui com base em valores por átomo.
    """
    total_fukui = 0.0
    total_fukushima = 0.0
    atom_count = 0
    active_atom_count = 0

    with open(pocket_file, 'r') as pdb:
        for line in pdb:
            if not line.startswith("ATOM"):
                continue
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            atom_id = f"{chain_id}:{res_seq}{res_name}:{atom_name}"

            if atom_id in fukui_dict:
                total_fukui += fukui_dict[atom_id]
            if atom_id in fukushima_dict:
                val_fuku = fukushima_dict[atom_id]
                total_fukushima += val_fuku
                if val_fuku != 0.0:
                    active_atom_count += 1
            atom_count += 1

    avg_fukui = total_fukui / atom_count if atom_count > 0 else None
    avg_fukushima = total_fukushima / atom_count if atom_count > 0 else None
    fukushima_no_norm = total_fukushima
    fukushima_sum_norm = total_fukushima / active_atom_count if active_atom_count > 0 else None

    return avg_fukui, avg_fukushima, fukushima_no_norm, fukushima_sum_norm


# ============================================================
# PIPELINE PRINCIPAL
# ============================================================
def compute_average(pdb_path, pocket_file, fukui_file, fukushima_file):
    residue_map = extract_residue_index_mapping(pdb_path)
    fukui_dict = load_bfactor_values(residue_map, fukui_file, descriptor_name="Fukui")
    fukushima_dict = load_bfactor_values(residue_map, fukushima_file, descriptor_name="Fukushima")
    return compute_values(fukui_dict, fukushima_dict, pocket_file)


def main(pdb_dir, pockets_dir, output_hocus, pdbbind_csv):
    print(f"[INICIANDO] Lendo arquivo de entrada: {pdbbind_csv}")
    df = pd.read_csv(pdbbind_csv)

    # ============================================================
    # 1️⃣ Classificação de Drogabilidade com base em KD/Ki
    # ============================================================
    colunas_possiveis = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)", "Kd/Ki (µM)", "Kd/Ki (uM)"]
    col_afinidade = next((c for c in df.columns if any(alt.lower() == c.lower() for alt in colunas_possiveis)), None)
    if col_afinidade is None:
        print("[ERRO] Coluna de afinidade KD/Ki não encontrada no CSV.")
        print("Colunas disponíveis:", df.columns.tolist())
        sys.exit(1)

    df["Categoria"] = df[col_afinidade].apply(classificar_drogabilidade)
    print("[INFO] Classificação de drogabilidade concluída:")
    print(df["Categoria"].value_counts(dropna=False))

    # ============================================================
    # 2️⃣ Adiciona colunas para descritores Fukui/Fukushima
    # ============================================================
    new_cols = ["Avg_Fukui", "Avg_Fukushima", "Fukushima_NoNorm", "Fukushima_SumNorm"]
    for col in new_cols:
        if col not in df.columns:
            df[col] = None

    # ============================================================
    # 3️⃣ Loop sobre os PDBs e cálculo dos descritores
    # ============================================================
    for file in os.listdir(pdb_dir):
        if not file.endswith(".pdb"):
            continue
        pdb_path = os.path.join(pdb_dir, file)
        protein_name = os.path.splitext(file)[0]
        protein_output_dir = os.path.join(output_hocus, protein_name)

        fukui_file = os.path.join(
            protein_output_dir, "reactivity", "primordia", f"{protein_name}_PDB_RD", f"{protein_name}_fukui_pot_zero.pdb"
        )
        fukushima_file = os.path.join(
            protein_output_dir, "reactivity", "primordia", f"{protein_name}_PDB_RD", f"{protein_name}_fukushima.pdb"
        )
        pocket_file = os.path.join(pockets_dir, f"{protein_name}_pocket.pdb")

        if not all(os.path.exists(p) for p in [fukui_file, fukushima_file, pocket_file]):
            print(f"[AVISO] Arquivos ausentes para {protein_name}, pulando.")
            continue

        avg_fukui, avg_fukushima, fuku_no_norm, fuku_sum_norm = compute_average(
            pdb_path, pocket_file, fukui_file, fukushima_file
        )

        mask = df["ID"].astype(str).str.lower() == protein_name.lower()
        if mask.any():
            df.loc[mask, "Avg_Fukui"] = avg_fukui
            df.loc[mask, "Avg_Fukushima"] = avg_fukushima
            df.loc[mask, "Fukushima_NoNorm"] = fuku_no_norm
            df.loc[mask, "Fukushima_SumNorm"] = fuku_sum_norm
            print(f"[OK] {protein_name} atualizado com descritores.")
        else:
            print(f"[INFO] ID {protein_name} não encontrado no CSV, pulando.")

    # ============================================================
    # 4️⃣ Exporta CSV final com classificação e descritores
    # ============================================================
    df.to_csv(pdbbind_csv, index=False)
    print("\n✅ [CONCLUÍDO] CSV atualizado salvo com classificação e descritores em:")
    print(pdbbind_csv)


# ============================================================
# ENTRY POINT
# ============================================================
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python3 drug_pdbbind_csv.py <pdb_dir> <pockets_dir> <output_hocus> <pdbbind_csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
