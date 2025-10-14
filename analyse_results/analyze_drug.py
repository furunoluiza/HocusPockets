import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, gaussian_kde
from sklearn.metrics import roc_curve, auc

def select_best_cavities(df, score_col="Fukushima_Score"):
    """
    Selects the cavity with the highest Fukushima_Score per protein.
    """
    df = df.copy()
    df = df.dropna(subset=[score_col])
    best_df = (
        df.loc[df.groupby("Protein_Code")[score_col].idxmax()]
        .reset_index(drop=True)
    )
    print(f"[INFO] Selecionadas {len(best_df)} cavidades representativas (maior {score_col}).")
    return best_df


def plot_violin_by_category(df, output_dir="output_druggability_violin"):
    os.makedirs(output_dir, exist_ok=True)
    palette = {"d": "#4caf50", "n": "#ef5350"}  # d = drogável, n = não drogável

    DESCRIPTOR_LABELS = {
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

    descriptors_to_plot = list(DESCRIPTOR_LABELS.keys())

    for desc in descriptors_to_plot:
        if desc not in df.columns:
            continue

        sub = df[["Category", desc]].dropna()
        if sub["Category"].nunique() < 2:
            continue

        desc_label = DESCRIPTOR_LABELS.get(desc, desc)
        drug_vals = sub[sub["Category"] == "d"][desc].values
        non_vals = sub[sub["Category"] == "n"][desc].values

        # Teste t de Welch
        tstat, pval = ttest_ind(drug_vals, non_vals, equal_var=False)
        if pval < 0.001: stars = "***"
        elif pval < 0.01: stars = "**"
        elif pval < 0.05: stars = "*"
        else: stars = "ns"

        plt.figure(figsize=(8, 6))
        sns.violinplot(
            x="Category", y=desc, data=sub,
            order=["d", "n"],
            hue="Category",
            palette=palette,
            inner="box",
            linewidth=1.2,
            legend=False
        )

        max_y = sub[desc].max()
        y_pos = max_y * 1.1
        h = max_y * 0.03
        plt.plot([0, 0, 1, 1], [y_pos, y_pos + h, y_pos + h, y_pos], lw=1.2, c='black')
        plt.text(0.5, y_pos + h * 1.2, stars, ha='center', va='bottom', fontsize=12)

        plt.title(desc_label, fontsize=14)
        plt.ylabel(desc_label)
        plt.xlabel("Categoria da proteína (d = drogável, n = não drogável)")
        plt.grid(axis="y", linestyle="--", alpha=0.5)
        plt.text(0.98, 0.98, f"p={pval:.2e}", transform=plt.gca().transAxes,
                 ha="right", va="top", fontsize=9, color="gray")

        plt.savefig(os.path.join(output_dir, f"violin_{desc}.png"), bbox_inches="tight")
        plt.close()
        print(f"[INFO] Violin plot gerado para {desc_label} (p={pval:.4g}, {stars})")


def plot_distributions_and_roc(df, output_dir="output_druggability"):
    os.makedirs(output_dir, exist_ok=True)
    df["Grupo"] = df["Category"].map({"d": "Drogável", "n": "Não drogável"})

    descriptors = [
        "Volume", "Radius", "Surface", "Diameter", "Smallest_Size",
        "Cylinder_Radius", "Hydrophobicity_Score", "Aromatic_Freq", "Otyr_Freq",
        "Carbon_Freq", "Fukui_Score", "Fukushima_Score",
        "Fukushima_Sum", "Fukushima_Sum_Norm", "Score"
    ]

    for score_col in descriptors:
        if score_col not in df.columns:
            continue

        label_col = score_col.replace("_", " ")
        label_col = label_col[0].upper() + label_col[1:]

        # Histograma
        plt.figure(figsize=(8, 6))
        for grupo, cor in [("Drogável", "green"), ("Não drogável", "red")]:
            subset = df.loc[df["Grupo"] == grupo, score_col].dropna()
            plt.hist(subset, bins=25, alpha=0.6, label=f"{grupo} (n={len(subset)})", color=cor)
        plt.xlabel(label_col)
        plt.ylabel("Número de proteínas")
        plt.title(f"Distribuição - {label_col}")
        plt.legend()
        plt.savefig(os.path.join(output_dir, f"hist_{score_col}.png"))
        plt.close()

        # Curva ROC (agora para todos)
        try:
            y_true = df["Category"].map({"d": 1, "n": 0}).values
            y_score = df[score_col].values
            if len(np.unique(y_score)) < 2:
                continue
            fpr, tpr, _ = roc_curve(y_true, y_score)
            roc_auc = auc(fpr, tpr)
        except Exception:
            roc_auc = np.nan
            fpr, tpr = [0, 1], [0, 1]

        plt.figure(figsize=(6, 6))
        plt.plot(fpr, tpr, color="navy", lw=2, label=f"AUC = {roc_auc:.2f}")
        plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
        plt.xlabel("Falsos Positivos")
        plt.ylabel("Verdadeiros Positivos")
        plt.title(f"Curva ROC - {label_col}")
        plt.legend(loc="lower right")
        plt.savefig(os.path.join(output_dir, f"roc_{score_col}.png"))
        plt.close()
        print(f"[INFO] Curva ROC gerada para {score_col} (AUC={roc_auc:.3f})")


def main(summary_csv):
    df = pd.read_csv(summary_csv)
    if "Protein_Code" not in df.columns or "Category" not in df.columns:
        raise ValueError("A tabela deve conter as colunas 'Protein_Code' e 'Category'.")

    best_df = select_best_cavities(df, score_col="Fukushima_Score")

    plot_violin_by_category(best_df)
    plot_distributions_and_roc(best_df)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_by_best_cavity.py <merged_cavities.csv>")
        sys.exit(1)
    summary_csv = sys.argv[1]
    main(summary_csv)
