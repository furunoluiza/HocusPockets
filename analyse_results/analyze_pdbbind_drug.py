import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def assign_druggability_category(kd_value):
    """
    Classify target as druggable (d), borderline (b), or non-druggable (n)
    based on KD/KI affinity thresholds in µM.
    """
    if pd.isna(kd_value):
        return np.nan

    try:
        kd_value = float(kd_value)
    except ValueError:
        return np.nan

    # Convert from nM → µM if necessary (e.g., >10,000 nM = 10 µM)
    if kd_value > 10000:
        kd_value = kd_value / 1000.0

    if kd_value < 0.3:
        return "d"  # druggable
    elif kd_value <= 1.0:
        return "b"  # borderline
    else:
        return "n"  # non-druggable


def plot_roc_one_vs_rest(df, descriptor, label, output_dir):
    """
    Generates One-vs-Rest ROC curves for d, b, and n categories.
    """
    os.makedirs(output_dir, exist_ok=True)
    sub = df[["Category", descriptor]].dropna()

    if sub["Category"].nunique() < 2:
        print(f"[AVISO] Descriptor {descriptor} não tem variação suficiente entre categorias.")
        return

    y_score = sub[descriptor].values
    categories = ["d", "b", "n"]
    colors = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}
    aucs = {}

    plt.figure(figsize=(7, 6))

    for cat in categories:
        y_true = (sub["Category"] == cat).astype(int).values
        if len(np.unique(y_true)) < 2:
            continue

        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        aucs[cat] = roc_auc

        plt.plot(fpr, tpr, lw=2, color=colors[cat], label=f"{cat.upper()} (AUC={roc_auc:.2f})")

    plt.plot([0, 1], [0, 1], color="gray", linestyle="--", lw=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC One-vs-Rest — {label}")
    plt.legend(loc="lower right")
    plt.tight_layout()

    out_path = os.path.join(output_dir, f"roc_ovr_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()

    if aucs:
        macro_auc = np.mean(list(aucs.values()))
        print(f"[INFO] {descriptor}: Macro-AUC={macro_auc:.3f}  "
              f"(d={aucs.get('d', np.nan):.3f}, b={aucs.get('b', np.nan):.3f}, n={aucs.get('n', np.nan):.3f})")


def main(summary_file):

    df = pd.read_excel(summary_file)

    # Detect KD/KI column, considering µ or u
    possible_cols = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    found_col = None
    for col in df.columns:
        if any(alt.lower() == col.lower() for alt in possible_cols):
            found_col = col
            break

    if found_col is None:
        print("[ERRO] Coluna de KD/KI não encontrada (esperado 'kd/ki (µM)' ou 'kd/ki (uM)').")
        print("Colunas detectadas:", list(df.columns))
        sys.exit(1)
    else:
        print(f"[INFO] Coluna de afinidade detectada: '{found_col}'")

    # Classify based on KD/KI
    df["Category"] = df[found_col].apply(assign_druggability_category)
    print("\n[INFO] Contagem de classes:")
    print(df["Category"].value_counts(dropna=False))

    # Descriptor columns present in your data
    DESCRIPTORS = {
        "Avg_Fukui": "Fukui médio",
        "Avg_Fukushima": "Fukushima médio",
        "Fukushima_NoNorm": "Fukushima total (soma bruta)",
        "Fukushima_SumNorm": "Fukushima normalizado (átomos ativos)",
    }

    out_dir = "output_druggability_ovr"
    os.makedirs(out_dir, exist_ok=True)

    for desc, label in DESCRIPTORS.items():
        if desc in df.columns:
            plot_roc_one_vs_rest(df, desc, label, out_dir)
        else:
            print(f"[AVISO] Coluna {desc} não encontrada no arquivo.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_pdbbind_drug_multiclass.py <pdbbind_descriptors.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
