#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc

sns.set_theme(
    style="white",
    palette="colorblind",
    font="Arial",
    font_scale=1.5,
    context="paper",
    rc={
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'lines.linewidth': 2,
        'axes.linewidth': 1.2,
    },
)


# ============================================================
# DRUGGABILITY CLASSIFICATION BASED ON KD/KI
# ============================================================
def classify_druggability(kd_value):
    """Classify target as Druggable (d), Borderline (b), or Non-druggable (n) based on KD/KI (µM)."""
    if pd.isna(kd_value):
        return np.nan
    try:
        kd_value = float(kd_value)
    except ValueError:
        return np.nan
    if kd_value > 10000:
        kd_value /= 1000.0
    if kd_value < 0.3:
        return "d"
    elif kd_value <= 1.0:
        return "b"
    else:
        return "n"


# ============================================================
# MULTICLASS ROC PLOT
# ============================================================
def plot_multiclass_roc(df, descriptor, label, output_dir):
    sub = df[["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[WARNING] No valid data for {descriptor}.")
        return

    colors = {"d": "#0072B2", "b": "#F0E442", "n": "#D55E00"}
    names = {"d": "Druggable", "b": "Borderline", "n": "Non-druggable"}

    plt.figure(figsize=(8, 6))
    aucs = {}

    for cat in ["d", "b", "n"]:
        y_true = (sub["Categoria"] == cat).astype(int)
        y_score = sub[descriptor]
        if len(np.unique(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        aucs[cat] = roc_auc
        n_pos = int(sum(y_true))
        n_neg = int(len(y_true) - n_pos)
        plt.plot(fpr, tpr, lw=2, color=colors[cat],
                 label=f"{names[cat]} (AUC={roc_auc:.2f}, N+={n_pos}, N-={n_neg})")

    plt.plot([0, 1], [0, 1], color="gray", lw=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend(loc="lower right")
    plt.tight_layout()
    out_path = os.path.join(output_dir, f"roc_multiclass_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Multiclass ROC saved: {out_path}")


# ============================================================
# BINARY ROC PLOT (D vs N)
# ============================================================
def plot_binary_roc(df, descriptor, label, output_dir):
    sub = df[df["Categoria"].isin(["d", "n"])][["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[WARNING] No valid data for D vs N in {descriptor}.")
        return

    y_true = (sub["Categoria"] == "d").astype(int)
    y_score = sub[descriptor]
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    n_pos = int(sum(y_true))
    n_neg = int(len(y_true) - n_pos)
    total_n = n_pos + n_neg

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, lw=2, color="#0072B2")
    plt.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.legend([f"AUC = {roc_auc:.2f}  |  N = {total_n}"], loc="lower right", frameon=False)
    plt.tight_layout()
    out_path = os.path.join(output_dir, f"roc_druggable_vs_non_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Binary ROC saved: {out_path}")


# ============================================================
# DENSITY PLOT
# ============================================================
def plot_density(df, descriptor, label, output_dir, binary=False):
    sub = df[["Categoria", descriptor]].dropna()
    if binary:
        sub = sub[sub["Categoria"].isin(["d", "n"])]

    colors = {"d": "#0072B2", "b": "#F0E442", "n": "#D55E00"}
    names = {"d": "Druggable", "b": "Borderline", "n": "Non-druggable"}

    plt.figure(figsize=(8, 6))

    for cat in ["d", "b", "n"]:
        subset = sub[sub["Categoria"] == cat][descriptor].dropna()
        if len(subset) <= 1:
            continue
        sns.kdeplot(subset, color=colors[cat], linewidth=2.5,
                    label=f"{names[cat]} (N={len(subset)})", cut=0,
                    common_norm=False)

    plt.xlabel(label)
    plt.ylabel("Probability density")
    plt.legend()
    plt.tight_layout()

    suffix = "_binary" if binary else "_multiclass"
    out_path = os.path.join(output_dir, f"density{suffix}_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Density plot saved: {out_path}")


# ============================================================
# MAIN
# ============================================================
def main(spreadsheet):
    print(f"[INFO] Reading file: {spreadsheet}")
    if not os.path.exists(spreadsheet):
        print(f"[ERROR] File not found: {spreadsheet}")
        sys.exit(1)

    if spreadsheet.endswith(".xlsx"):
        df = pd.read_excel(spreadsheet)
    else:
        df = pd.read_csv(spreadsheet)

    possible_columns = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    affinity_column = next((c for c in df.columns if any(alt.lower() == c.lower() for alt in possible_columns)), None)
    if affinity_column is None:
        print("[ERROR] Affinity column not found.")
        print("Available columns:", df.columns.tolist())
        sys.exit(1)

    df["Categoria"] = df[affinity_column].apply(classify_druggability)
    print("[INFO] Category distribution:")
    print(df["Categoria"].value_counts(dropna=False))

    descriptors = {
        "Avg_Fukui": "Fukui average (a.u.)",
        "Avg_Fukushima": "Fukushima average (a.u.)",
        "Fukushima_NoNorm": "Fukushima total (a.u.)",
        "Fukushima_SumNorm": "Fukushima normalized (a.u.)",
    }

    output_dir = "output_druggability_final"
    os.makedirs(output_dir, exist_ok=True)

    for desc, label in descriptors.items():
        if desc not in df.columns:
            print(f"[WARNING] Column {desc} not found.")
            continue
        print(f"\n[INFO] Processing {desc} ({label})")

        plot_multiclass_roc(df, desc, label, output_dir)
        plot_density(df, desc, label, output_dir, binary=False)
        plot_binary_roc(df, desc, label, output_dir)
        plot_density(df, desc, label, output_dir, binary=True)

    print("\n[DONE] All plots generated in:", output_dir)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_pdbbind_drug.py <pdbbind_spreadsheet.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
