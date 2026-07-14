import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from sklearn.metrics import roc_curve, auc

# -----------------------------------------
# Funções auxiliares
# -----------------------------------------
def select_best_cavities(df, descriptor):
    """
    Selects the cavity with the highest value of a given descriptor per protein.
    """
    df = df.copy()
    df = df.dropna(subset=[descriptor])
    best_df = df.loc[df.groupby("Protein_Code")[descriptor].idxmax()].reset_index(drop=True)
    return best_df


def plot_descriptor_analysis(df, descriptor, label, output_dir_violin, output_dir_roc):
    """
    Generates violin, histogram, and ROC plots for a given descriptor,
    using the best cavity per protein according to that descriptor.
    """
    os.makedirs(output_dir_violin, exist_ok=True)
    os.makedirs(output_dir_roc, exist_ok=True)

    palette = {"d": "#4caf50", "n": "#ef5350"}  # d = drogável, n = não drogável

    # Seleciona a melhor cavidade por proteína para o descritor atual
    best_df = select_best_cavities(df, descriptor)

    # Prepara dados
    sub = best_df[["Category", descriptor]].dropna()
    if sub["Category"].nunique() < 2:
        print(f"[AVISO] Descriptor {descriptor} não tem variação suficiente entre categorias.")
        return

    drug_vals = sub[sub["Category"] == "d"][descriptor].values
    non_vals = sub[sub["Category"] == "n"][descriptor].values

    # Welch t-test
    tstat, pval = ttest_ind(drug_vals, non_vals, equal_var=False)
    if pval < 0.001: stars = "***"
    elif pval < 0.01: stars = "**"
    elif pval < 0.05: stars = "*"
    else: stars = "ns"

    # ------------------------------
    # Gráfico de violino
    # ------------------------------
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        x="Category", y=descriptor, data=sub,
        order=["d", "n"], hue="Category",
        palette=palette, inner="box", linewidth=1.2, legend=False
    )

    max_y = sub[descriptor].max()
    y_pos = max_y * 1.1
    h = max_y * 0.03
    plt.plot([0, 0, 1, 1], [y_pos, y_pos + h, y_pos + h, y_pos], lw=1.2, c='black')
    plt.text(0.5, y_pos + h * 1.2, stars, ha='center', va='bottom', fontsize=12)

    plt.title(label, fontsize=14)
    plt.ylabel(label)
    plt.xlabel("Categoria (d = drogável, n = não drogável)")
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.text(0.98, 0.98, f"p={pval:.2e}", transform=plt.gca().transAxes,
             ha="right", va="top", fontsize=9, color="gray")

    plt.savefig(os.path.join(output_dir_violin, f"violin_{descriptor}.png"), bbox_inches="tight")
    plt.close()

    # ------------------------------
    # Histograma + ROC
    # ------------------------------
    plt.figure(figsize=(8, 6))
    for cat, color in [("d", "green"), ("n", "red")]:
        subset = sub.loc[sub["Category"] == cat, descriptor]
        plt.hist(subset, bins=25, alpha=0.6, label=f"{cat} (n={len(subset)})", color=color)
    plt.xlabel(label)
    plt.ylabel("Número de proteínas")
    plt.title(f"Distribuição - {label}")
    plt.legend()
    plt.savefig(os.path.join(output_dir_roc, f"hist_{descriptor}.png"))
    plt.close()

    # ROC curve
    y_true = sub["Category"].map({"d": 1, "n": 0}).values
    y_score = sub[descriptor].values
    if len(np.unique(y_score)) > 1:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
    else:
        fpr, tpr, roc_auc = [0, 1], [0, 1], np.nan

    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, color="navy", lw=2, label=f"AUC = {roc_auc:.2f}")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Falsos Positivos")
    plt.ylabel("Verdadeiros Positivos")
    plt.title(f"Curva ROC - {label}")
    plt.legend(loc="lower right")
    plt.savefig(os.path.join(output_dir_roc, f"roc_{descriptor}.png"))
    plt.close()

    print(f"[INFO] {descriptor}: p={pval:.4g}, AUC={roc_auc:.3f}, estrelas={stars}")


def main(summary_csv):
    df = pd.read_csv(summary_csv)
    if "Protein_Code" not in df.columns or "Category" not in df.columns:
        raise ValueError("O CSV deve conter as colunas 'Protein_Code' e 'Category'.")

    DESCRIPTORS = {
        "Volume": "Volume (Å³)",
        "Radius": "Raio (Å)",
        "Surface": "Área de superfície (Å²)",
        "Diameter": "Diâmetro (Å)",
        "Smallest_Size": "Menor dimensão (Å)",
        "Cylinder_Radius": "Raio do cilindro envolvente (Å)",
        "Hydrophobicity_Score": "Hidrofobicidade média",
        "Aromatic_Freq": "Frequência de resíduos aromáticos",
        "Otyr_Freq": "Frequência de resíduos de tirosina",
        "Carbon_Freq": "Frequência de átomos de carbono",
        "Fukui_Score": "Fukui médio",
        "Fukushima_Score": "Fukushima médio",
        "Fukushima_Sum": "Fukushima total (soma bruta)",
        "Fukushima_Sum_Norm": "Fukushima normalizado (átomos ativos)",
        "Score": "Pontuação clássica"
    }

    out_violin = "output_druggability_violin"
    out_roc = "output_druggability"

    for desc, label in DESCRIPTORS.items():
        if desc in df.columns:
            plot_descriptor_analysis(df, desc, label, out_violin, out_roc)
        else:
            print(f"[AVISO] Coluna {desc} não encontrada no CSV.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_by_best_descriptor.py <merged_cavities.csv>")
        sys.exit(1)
    summary_csv = sys.argv[1]
    main(summary_csv)
