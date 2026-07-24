#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ttest_ind
from sklearn.metrics import roc_curve, auc, matthews_corrcoef

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

def plot_boxplots_descriptors(summary_csv, output_dir="output_resultados"):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(summary_csv)
    if "IsReference" not in df.columns:
        print("[ERRO] Coluna 'IsReference' não encontrada.")
        return

    # Forçar IsReference → booleano
    df["IsReference"] = df["IsReference"].astype(str).str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    ).fillna(False)
    df["Classe"] = df["IsReference"].map({True: "Positivo", False: "Negativo"})

    descriptors = [
        "Volume", "Radius", "Surface", "Diameter", "Smallest_Size",
        "Cylinder_radius", "Hydrophobicity_Score", "Aromatic_freq",
        "Otyr_freq", "Fukushima_Score", "Fukui_Score"
    ]

    palette = {"Positivo": "#4caf50", "Negativo": "#ef5350"}

    for desc in descriptors:
        if desc not in df.columns:
            continue

        sub = df[["Classe", desc]].dropna()
        if sub["Classe"].nunique() < 2:
            continue

        # Separar classes
        drug_vals = sub[sub["Classe"] == "Positivo"][desc].values
        non_vals = sub[sub["Classe"] == "Negativo"][desc].values

        # Teste t
        tstat, pval = ttest_ind(drug_vals, non_vals, equal_var=False)
        if pval < 0.001: stars = "***"
        elif pval < 0.01: stars = "**"
        elif pval < 0.05: stars = "*"
        else: stars = "ns"

        plt.figure(figsize=(6,5))
        ax = sns.boxplot(
            x="Classe", y=desc, data=sub,
            order=["Positivo", "Negativo"],
            notch=False, medianprops={"color": "red", "linewidth": 2},
            boxprops={"facecolor": "none"}
        )
        for patch, label in zip(ax.patches, ["Positivo", "Negativo"]):
            patch.set_facecolor(palette[label])

        # Barra de significância
        max_y = sub[desc].max()
        y_pos = max_y * 1.05
        h = max_y * 0.03
        ax.plot([0, 0, 1, 1], [y_pos, y_pos+h, y_pos+h, y_pos], lw=1.5, c='black')
        ax.text(0.5, y_pos+h*1.2, stars, ha='center', va='bottom', fontsize=12)

        n_drug = len(drug_vals)
        n_non = len(non_vals)
        plt.ylabel(desc)
        plt.xlabel("Classe")

        out_file = os.path.join(output_dir, f"boxplot_{desc}.png")
        plt.savefig(out_file, bbox_inches="tight")
        plt.close()
        print(f"[INFO] Boxplot {desc} gerado (p={pval:.4g}, {stars})")


DESCRIPTOR_LABELS = {
    "Volume": "Volume",
    "Radius": "Radius",
    "Surface": "Surface area",
    "Diameter": "Diameter",
    "Smallest_Size": "Smallest dimension",
    "Cylinder_Radius": "Minimum enclosing cylinder radius",
    "Hydrophobicity_Score": "Average hydrophobicity",
    "Aromatic_Freq": "Aromatic residue frequency",
    "Otyr_Freq": "Tyrosine (OH) frequency",
    "Carbon_Freq": "Carbon atom frequency",
    "Fukushima_Score": "Fukushima average",
    "Fukui_Score": "Fukui average",
    "Fukushima_Sum": "Fukushima total (sum)",
    "Fukushima_Sum_Norm": "Fukushima normalized (active atoms)",
}


def plot_violin_descriptors(summary_csv, output_dir="violin_plots"):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(summary_csv)
    if "IsReference" not in df.columns:
        print("[ERROR] Column 'IsReference' not found.")
        return

    df["IsReference"] = df["IsReference"].astype(str).str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    ).fillna(False)
    df["Class"] = df["IsReference"].map({True: "Reference (EC)", False: "Non-reference (NEC)"})

    descriptors_to_plot = [
        "Volume", "Radius", "Surface", "Diameter", "Smallest_Size",
        "Cylinder_Radius", "Hydrophobicity_Score", "Aromatic_Freq", "Carbon_Freq",
        "Otyr_Freq", "Fukui_Score", "Fukushima_Score", "Fukushima_Sum", "Fukushima_Sum_Norm"]

    for desc in descriptors_to_plot:
        if desc not in df.columns:
            continue

        desc_label = DESCRIPTOR_LABELS.get(desc, desc)

        sub = df[["Class", desc]].dropna()
        if sub["Class"].nunique() < 2:
            continue

        drug_vals = sub[sub["Class"] == "Reference (EC)"][desc].values
        non_vals = sub[sub["Class"] == "Non-reference (NEC)"][desc].values

        tstat, pval = ttest_ind(drug_vals, non_vals, equal_var=False)
        if pval < 0.001: stars = "***"
        elif pval < 0.01: stars = "**"
        elif pval < 0.05: stars = "*"
        else: stars = "ns"

        plt.figure(figsize=(8, 6))
        ax = sns.violinplot(
            x="Class", y=desc, data=sub,
            order=["Reference (EC)", "Non-reference (NEC)"],
            inner="box",
            linewidth=1.5,
            legend=False
        )

        max_y = sub[desc].max()
        y_pos = max_y * 1.15
        h = max_y * 0.03

        if y_pos > ax.get_ylim()[1]:
            ax.set_ylim(top=y_pos * 1.15)

        ax.plot([0, 0, 1, 1], [y_pos, y_pos + h, y_pos + h, y_pos], lw=1.5, c='black')
        ax.text(0.5, y_pos + h * 1.2, stars, ha='center', va='bottom', fontsize=14)

        n_ref = len(drug_vals)
        n_non = len(non_vals)
        p_value_str = f"p={pval:.2e}"

        plt.text(0.98, 0.98, f"Reference n={n_ref}\nNon-reference n={n_non}\n{p_value_str}",
                 transform=plt.gca().transAxes,
                 ha="right", va="top", fontsize=11, color="gray")

        plt.ylabel(desc_label)
        plt.xlabel("Class")

        out_file = os.path.join(output_dir, f"violinplot_{desc}.png")
        plt.savefig(out_file, bbox_inches="tight")
        plt.close()
        print(f"[INFO] Violin plot {desc_label} saved (p={pval:.4g}, {stars})")

def plot_distributions_and_roc(summary_csv, output_dir="output_distributions"):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(summary_csv)

    df["Group"] = df["IsReference"].map({True: "Reference (EC)", False: "Non-reference (NEC)"})

    descriptors = [
        ("Fukui_Score", "Fukui average (a.u.)"),
        ("Fukushima_Score", "Fukushima average (a.u.)"),
        ("Fukushima_Sum", "Fukushima total (a.u.)"),
        ("Fukushima_Sum_Norm", "Fukushima normalized (a.u.)"),
        ("Score", "Classic score (a.u.)"),
        ("Carbon_Freq", "Carbon atom frequency"),
        ("Volume", "Volume (Å³)"),
        ("Radius", "Radius (Å)"),
        ("Surface", "Surface area (Å²)"),
        ("Diameter", "Diameter (Å)"),
        ("Smallest_Size", "Smallest dimension (Å)"),
        ("Cylinder_Radius", "Cylinder radius (Å)"),
        ("Hydrophobicity_Score", "Average hydrophobicity (a.u.)"),
        ("Aromatic_Freq", "Aromatic residue frequency"),
        ("Otyr_Freq", "Tyrosine (OH) frequency")
    ]

    groups = [
        ("Reference (EC)", "EC", "#0072B2"),
        ("Non-reference (NEC)", "NEC", "#D55E00")
    ]

    for score_col, label_col in descriptors:
        if score_col not in df.columns:
            print(f"[WARNING] Column '{score_col}' not found. Skipping.")
            continue

        # ----------------------------------------------------------------------
        # HISTOGRAM
        # ----------------------------------------------------------------------
        plt.figure(figsize=(8, 6))

        for grupo_real, nome_legenda, cor in groups:
            subset = df.loc[df["Group"] == grupo_real, score_col].dropna()
            if len(subset) == 0:
                continue
            sns.histplot(subset, bins=30, alpha=0.6, color=cor,
                         label=f"{nome_legenda} (N = {len(subset)})", element="step")

        plt.xlabel(label_col)
        plt.ylabel("Number of cavities")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"hist_{score_col}.png"), dpi=300)
        plt.close()
        print(f"[OK] Histogram saved: {score_col}")

        # ----------------------------------------------------------------------
        # KDE (Density curves)
        # ----------------------------------------------------------------------
        plt.figure(figsize=(8, 6))

        for grupo_real, nome_legenda, cor in groups:
            subset = df.loc[df["Group"] == grupo_real, score_col].dropna()
            if len(subset) > 1:
                sns.kdeplot(subset, color=cor, linewidth=2.5,
                            label=f"{nome_legenda} (N = {len(subset)})", cut=0,
                            common_norm=False)

        plt.xlabel(label_col)
        plt.ylabel("Probability density")
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"density_{score_col}.png"), dpi=300)
        plt.close()
        print(f"[OK] KDE curve saved: {score_col}")

        # ROC
        y_true = df["IsReference"].astype(int).values
        y_score = df[score_col].values
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        plt.figure(figsize=(6, 6))
        plt.plot(fpr, tpr, color="#0072B2", lw=2.5)
        plt.plot([0, 1], [0, 1], color="gray", linestyle="--", lw=1.5)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend([f"AUC = {roc_auc:.3f}  |  N = {len(y_true)}"], loc="lower right", frameon=False)
        plt.savefig(os.path.join(output_dir, f"roc_{score_col}.png"))
        plt.close()


def plot_cumulative_from_scores(summary_csv, output_dir="output_cumulative"):
    """
    Generates grouped bar charts (% of proteins) for Top-1/2/3 cumulative
    hit rate using scores from all_cavities_summary.csv.
    """
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(summary_csv)

    df["IsReference"] = (
        df["IsReference"]
        .astype(str).str.strip().str.lower()
        .map({"true": True, "false": False, "1": True, "0": False})
        .fillna(False)
    )

    methods = [
        ("Score", "Classic score (a.u.)", "classic_score_cumulative_bar.png"),
        ("Fukui_Score", "Fukui average",  "fukui_cumulative_bar.png"),
        ("Fukushima_Score", "Fukushima average", "fukushima_cumulative_bar.png"),
        ("Fukushima_Sum", "Fukushima total (a.u.)", "fukushima_sum_cumulative_bar.png"),
        ("Fukushima_Sum_Norm", "Fukushima normalized (a.u.)", "fukushima_sum_norm_cumulative_bar.png")
    ]

    for score_col, label, outname in methods:
        if score_col not in df.columns:
            print(f"[WARNING] Column {score_col} not found in {summary_csv}; skipping {label}.")
            continue

        df_score = df.copy()
        df_score[score_col] = pd.to_numeric(df_score[score_col], errors="coerce")
        df_score = df_score.dropna(subset=[score_col, "Protein", "Cavity_ID"])

        ks = [1, 2, 3]
        hits = {k: 0 for k in ks}
        considered_proteins = 0

        for protein, sub in df_score.groupby("Protein", dropna=True):
            idx = sub.groupby("Cavity_ID")[score_col].idxmax()
            reps = sub.loc[idx].copy()

            if reps.empty:
                continue

            considered_proteins += 1

            if "RowID" in reps.columns:
                reps = reps.sort_values([score_col, "RowID"], ascending=[False, True])
            else:
                reps = reps.sort_values(score_col, ascending=False)

            for k in ks:
                topk = reps.head(k)
                if topk["IsReference"].any():
                    hits[k] += 1

        if considered_proteins == 0:
            print(f"[WARNING] No valid proteins for {label}.")
            continue

        hits_list   = [100.0 * hits[k] / considered_proteins for k in ks]
        misses_list = [100.0 * (considered_proteins - hits[k]) / considered_proteins for k in ks]

        bar_data = pd.DataFrame({
            "Top": [f"Top {k}" for k in ks] * 2,
            "Group": ["EC"] * 3 + ["NEC"] * 3,
            "Proteins (%)": hits_list + misses_list,
        })

        plt.figure(figsize=(8, 6))
        ax = sns.barplot(data=bar_data, x="Top", y="Proteins (%)", hue="Group")

        for container in ax.containers:
            for bar in container:
                yval = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2, yval + 1,
                    f"{yval:.1f}%",
                    ha="center", va="bottom", fontsize=12, color="black"
                )

        ax.set_ylim(0, 105)
        ax.legend(loc="upper left", title="Group")

        ax.text(
            0.98, 0.95, f"N = {considered_proteins}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=11, color="gray"
        )

        out_path = os.path.join(output_dir, outname)
        plt.savefig(out_path)
        plt.close()
        print(f"[INFO] Cumulative bar chart ({label}) saved: {out_path}")


def plot_binary_classification_table(summary_csv, output_dir="output_classification"):
    """
    Generates a CSV with binary classification metrics (AUC, accuracy, sensitivity,
    specificity, F1, MCC, enrichment) for each quantum descriptor,
    using the optimal threshold by Youden index.
    """
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(summary_csv)

    df["IsReference"] = (
        df["IsReference"]
        .astype(str).str.strip().str.lower()
        .map({"true": True, "false": False, "1": True, "0": False})
        .fillna(False)
    )

    quantum_descriptors = ["Fukui_Score", "Fukushima_Score", "Fukushima_Sum"]
    y_true = df["IsReference"].astype(int).values

    results = []
    for col in quantum_descriptors:
        if col not in df.columns:
            print(f"[WARNING] Column {col} not found. Skipping.")
            continue

        y_score = pd.to_numeric(df[col], errors="coerce").values
        mask = ~np.isnan(y_score)
        if mask.sum() == 0:
            continue

        y_t = y_true[mask]
        y_s = y_score[mask]

        if len(np.unique(y_t)) < 2:
            print(f"[WARNING] {col} does not have both classes for binary classification.")
            continue

        # ROC + Youden threshold
        fpr, tpr, thresholds = roc_curve(y_t, y_s)
        roc_auc = auc(fpr, tpr)
        youden = tpr - fpr
        best_idx = np.argmax(youden)
        best_threshold = thresholds[best_idx]

        # Binary classification
        y_pred = (y_s >= best_threshold).astype(int)

        # Confusion matrix
        tp = np.sum((y_pred == 1) & (y_t == 1))
        fp = np.sum((y_pred == 1) & (y_t == 0))
        tn = np.sum((y_pred == 0) & (y_t == 0))
        fn = np.sum((y_pred == 0) & (y_t == 1))

        n_total = tp + fp + tn + fn
        accuracy  = (tp + tn) / n_total if n_total > 0 else 0
        sensitivity    = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity   = tn / (tn + fp) if (tn + fp) > 0 else 0
        precision  = tp / (tp + fp) if (tp + fp) > 0 else 0
        f1        = 2 * precision * sensitivity / (precision + sensitivity) if (precision + sensitivity) > 0 else 0
        mcc       = matthews_corrcoef(y_t, y_pred)
        prevalence = (tp + fn) / n_total if n_total > 0 else 1
        enrichment  = precision / prevalence if prevalence > 0 and precision > 0 else 1.0

        results.append({
            "Descriptor": col,
            "AUC": round(roc_auc, 4),
            "Threshold": round(best_threshold, 4),
            "Accuracy": round(accuracy, 4),
            "Sensitivity": round(sensitivity, 4),
            "Specificity": round(specificity, 4),
            "F1": round(f1, 4),
            "MCC": round(mcc, 4),
            "Enrichment": round(enrichment, 4),
        })

    if not results:
        print("[WARNING] No quantum descriptors found for binary classification table.")
        return

    out_path = os.path.join(output_dir, "binary_classification_table.csv")
    pd.DataFrame(results).to_csv(out_path, index=False)
    print(f"[INFO] Binary classification table saved: {out_path}")

def main(summary_csv):
    plot_violin_descriptors(summary_csv)
    plot_distributions_and_roc(summary_csv)
    plot_cumulative_from_scores(summary_csv)
    plot_binary_classification_table(summary_csv)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python analyze_summary.py <all_cavities_summary.csv>")
        sys.exit(1)
    summary_csv = sys.argv[1]
    main(summary_csv)
