#!/usr/bin/env python3
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import gaussian_kde, ttest_ind
from sklearn.metrics import roc_curve, auc

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
        plt.title(f"{desc}\n(Positivo n={n_drug}, Negativo n={n_non})")
        plt.ylabel(desc)
        plt.xlabel("Classe")

        out_file = os.path.join(output_dir, f"boxplot_{desc}.png")
        plt.savefig(out_file, bbox_inches="tight")
        plt.close()
        print(f"[INFO] Boxplot {desc} gerado (p={pval:.4g}, {stars})")

def plot_violin_descriptors(summary_csv, output_dir="output_resultados_violino"):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(summary_csv)
    if "IsReference" not in df.columns:
        print("[ERRO] Coluna 'IsReference' não encontrada.")
        return

    # Mapeamento para traduzir o nome das colunas no gráfico
    # Chave: Nome da coluna no CSV (variável 'desc')
    # Valor: Rótulo a ser usado no gráfico
    DESCRIPTOR_LABELS = {
        "Volume": "Volume",
        "Radius": "Raio",
        "Surface": "Área de superfície",
        "Diameter": "Diâmetro",
        "Smallest_Size": "Menor dimensão",
        "Cylinder_Radius": "Raio do menor cilindro envolvente",
        "Hydrophobicity_Score": "Hidrofobicidade média",
        "Aromatic_Freq": "Frequência de resíduos aromáticos",
        "Otyr_Freq": "Frequência de resíduos de tirosina",
        "Carbon_Freq": "Frequência de átomos de carbono",
        "Fukushima_Score": "Fukushima Média",
        "Fukui_Score": "Fukui Média",
        "Fukushima_Sum": "Fukushima Total (soma bruta)",
        "Fukushima_Sum_Norm": "Fukushima Normalizado (átomos ativos)"
    }
    
    # As colunas que não têm tradução específica no pedido (Fukushima_Score, Fukui_Score) foram mantidas como estão.

    # 1. Preparação de Dados
    
    # Forçar IsReference → booleano e criar coluna 'Classe'
    df["IsReference"] = df["IsReference"].astype(str).str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    ).fillna(False)
    df["Classe"] = df["IsReference"].map({True: "Positivo", False: "Negativo"})

    descriptors_to_plot = [
        "Volume", "Radius", "Surface", "Diameter", "Smallest_Size",
        "Cylinder_Radius", "Hydrophobicity_Score", "Aromatic_Freq", "Carbon_Freq",
        "Otyr_Freq", "Fukui_Score", "Fukushima_Score", "Fukushima_Sum", "Fukushima_Sum_Norm"]

    palette = {"Positivo": "#4caf50", "Negativo": "#ef5350"}

    # 2. Geração dos Gráficos
    for desc in descriptors_to_plot:
        if desc not in df.columns:
            continue

        # Obtém o rótulo traduzido para o gráfico
        desc_label = DESCRIPTOR_LABELS.get(desc, desc)

        sub = df[["Classe", desc]].dropna()
        if sub["Classe"].nunique() < 2:
            continue

        # Separar classes para o teste t
        drug_vals = sub[sub["Classe"] == "Positivo"][desc].values
        non_vals = sub[sub["Classe"] == "Negativo"][desc].values

        # Teste t de Student (Welch's t-test, assuming unequal variances)
        tstat, pval = ttest_ind(drug_vals, non_vals, equal_var=False)
        
        # Atribuição de estrelas de significância
        if pval < 0.001: stars = "***"
        elif pval < 0.01: stars = "**"
        elif pval < 0.05: stars = "*"
        else: stars = "ns"

        # 3. Plotagem do Gráfico de Violino
        plt.figure(figsize=(8, 6))
        
        # Gráfico de Violino (substituindo boxplot)
        ax = sns.violinplot(
            x="Classe", y=desc, data=sub,
            order=["Positivo", "Negativo"],
            hue="Classe",
            palette=palette,
            inner="box", # Desenha um boxplot dentro do violino
            linewidth=1.5,
            legend=False
        )

        # Barra de significância (mantida)
        max_y = sub[desc].max()
        y_pos = max_y * 1.15
        h = max_y * 0.03
        
        # Ajuste de altura para evitar sobreposição em casos extremos
        if y_pos > ax.get_ylim()[1]:
            ax.set_ylim(top=y_pos * 1.15)
            
        ax.plot([0, 0, 1, 1], [y_pos, y_pos + h, y_pos + h, y_pos], lw=1.5, c='black')
        ax.text(0.5, y_pos + h * 1.2, stars, ha='center', va='bottom', fontsize=12)

        # 4. Rótulos, Títulos e Salvamento
        n_drug = len(drug_vals)
        n_non = len(non_vals)
        
        # Título com o nome do descritor e a contagem de amostras (discreta)
        plt.title(f"{desc_label}", fontsize=14)

        p_value_str = f"p={pval:.2e}"
        
        # Anotação discreta da contagem (Positivos e Negativos)
        plt.text(0.98, 0.98, f"Positivo n={n_drug}\nNegativo n={n_non}\n{p_value_str}", 
                 transform=plt.gca().transAxes,
                 ha="right", va="top", fontsize=9, color="gray")
        
        plt.ylabel(desc_label)
        plt.xlabel("Classe")

        plt.grid(axis="y", linestyle="--", alpha=0.5)

        out_file = os.path.join(output_dir, f"violinplot_{desc}.png")
        plt.savefig(out_file, bbox_inches="tight")
        plt.close()
        print(f"[INFO] Gráfico de violino {desc_label} gerado (p={pval:.4g}, {stars})")

def plot_distributions_and_roc(summary_csv, output_dir="output_resultados"):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(summary_csv)
    df["Grupo"] = df["IsReference"].map({True: "Positivo", False: "Negativo"})

    for score_col, label_col in [
        ("Fukui_Score", "Fukui médio (a.u.)"),                         # adimensional
        ("Fukushima_Score", "Fukushima médio (a.u.)"),                 # adimensional
        ("Fukushima_Sum", "Fukushima total (a.u.)"),                   # adimensional
        ("Fukushima_Sum_Norm", "Fukushima normalizado (a.u.)"),        # adimensional
        ("Score", "Pontuação clássica (a.u.)"),                        # adimensional
        ("Carbon_Freq", "Frequência de átomos de carbono (%)"),        # porcentagem
        ("Volume", "Volume (Å³)"),                                     # angstrom cúbico
        ("Radius", "Raio (Å)"),                                        # angstrom
        ("Surface", "Área de superfície (Å²)"),                        # angstrom quadrado
        ("Diameter", "Diâmetro (Å)"),                                  # angstrom
        ("Smallest_Size", "Menor dimensão (Å)"),                       # angstrom
        ("Cylinder_Radius", "Menor raio do cilindro envolvente (Å)"),  # angstrom
        ("Hydrophobicity_Score", "Hidrofobicidade média (a.u.)"),      # adimensional (escala arbitrária)
        ("Aromatic_Freq", "Frequência de resíduos aromáticos (%)"),     # porcentagem
        ("Otyr_Freq", "Frequência de resíduos de tirosina (%)")         # porcentagem
    ]:

        if score_col not in df.columns:
            continue

        # Histograma
        plt.figure(figsize=(8,6))
        total_n = 0
        for grupo, cor in [("Positivo", "green"), ("Negativo", "red")]:
            subset = df.loc[df["Grupo"] == grupo, score_col].dropna()
            total_n += len(subset)
            plt.hist(subset, bins=30, alpha=0.6, label=f"{grupo} (n={len(subset)})", color=cor)
        plt.xlabel(label_col)
        plt.ylabel("Número de cavidades")
        plt.title(f"Histograma - {label_col}")
        plt.legend()
        plt.text(0.98, 0.95, f"Total n={total_n}", transform=plt.gca().transAxes,
                 ha="right", va="top", fontsize=9, color="gray")
        plt.savefig(os.path.join(output_dir, f"hist_{score_col}.png"))
        plt.close()

        # KDE
        plt.figure(figsize=(8,6))
        for grupo, cor in [("CCE", "green"), ("CNCE", "red")]:
            subset = df.loc[df["Grupo"] == grupo, score_col].dropna()
            if len(subset) > 1:
                kde = gaussian_kde(subset)
                x_vals = np.linspace(subset.min(), subset.max(), 300)
                densidade = kde(x_vals)
                plt.plot(x_vals, densidade, color=cor, lw=2,
                         label=f"{grupo} (N = {len(subset)})")
        plt.xlabel(label_col)
        plt.ylabel("Densidade de probabilidade")
        plt.title(label_col)
        plt.legend(frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"density_{score_col}.png"))
        plt.close()

        # ROC
        y_true = df["IsReference"].astype(int).values
        y_score = df[score_col].values
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        plt.figure(figsize=(6,6))
        plt.plot(fpr, tpr, color="navy", lw=2)
        plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
        plt.xlabel("Falsos Positivos")
        plt.ylabel("Verdadeiros Positivos")
        plt.title(label_col)
        plt.legend([f"AUC = {roc_auc:.2f}  |  N = {total_n}"], loc="lower right", frameon=False)
        plt.savefig(os.path.join(output_dir, f"roc_{score_col}.png"))
        plt.close()


def plot_cumulative_from_scores(summary_csv, output_dir="output_resultados"):
    """
    Gera gráficos de barras (em %) de acerto cumulativo até Top1/Top2/Top3
    usando os scores do all_cavities_summary.csv.
    Regras:
      - Ranking por Protein, separadamente para cada score (Fukui/Fukushima).
      - Antes de ranquear, mantém Cavity_IDs únicos, pegando a linha de MAIOR score
        daquele método por Cavity_ID (para evitar repetir a mesma cavidade).
      - 'Acerto' para Top-k se QUALQUER cavidade entre as Top-k tem IsReference=True.
    """
    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(summary_csv)

    # Robustez: força IsReference para bool
    df["IsReference"] = (
        df["IsReference"]
        .astype(str).str.strip().str.lower()
        .map({"true": True, "false": False, "1": True, "0": False})
        .fillna(False)
    )

    methods = [
        ("Score", "Pontuação clássica (a.u.)", "classico_cumulativo_bar.png"),
        ("Fukui_Score", "Fukui médio",  "fukui_cumulativo_bar.png"),
        ("Fukushima_Score", "Fukushima médio", "fukushima_cumulativo_bar.png"),
        ("Fukushima_Sum", "Fukushima total (a.u.)", "fukushima_sum_cumulativo_bar.png"),
        ("Fukushima_Sum_Norm", "Fukushima normalizado (a.u.)", "fukushima_sum_norm_cumulativo_bar.png")
    ]

    for score_col, label, outname in methods:
        if score_col not in df.columns:
            print(f"[AVISO] Coluna {score_col} não encontrada em {summary_csv}; pulando {label}.")
            continue

        # converte score para numérico e ignora NaN
        df_score = df.copy()
        df_score[score_col] = pd.to_numeric(df_score[score_col], errors="coerce")
        df_score = df_score.dropna(subset=[score_col, "Protein", "Cavity_ID"])

        # counts por Top-k
        ks = [1, 2, 3]
        hits = {k: 0 for k in ks}
        considered_proteins = 0

        for protein, sub in df_score.groupby("Protein", dropna=True):
            # dedup por Cavity_ID mantendo MAIOR score desse método para cada cavidade
            # idxmax entrega o índice da linha com maior score dentro de cada Cavity_ID
            # (tie-break natural fica pelo primeiro máximo encontrado)
            idx = sub.groupby("Cavity_ID")[score_col].idxmax()
            reps = sub.loc[idx].copy()

            # Se não sobrou nada útil para ranquear, pula esta proteína
            if reps.empty:
                continue

            considered_proteins += 1

            # Ordena por score desc; tie-break opcional por RowID para estabilidade
            if "RowID" in reps.columns:
                reps = reps.sort_values([score_col, "RowID"], ascending=[False, True])
            else:
                reps = reps.sort_values(score_col, ascending=False)

            # Para cada k, verifica se alguma Top-k é referência
            for k in ks:
                topk = reps.head(k)
                if topk["IsReference"].any():
                    hits[k] += 1

        if considered_proteins == 0:
            print(f"[AVISO] Nenhuma proteína válida para {label}.")
            continue

        # Calcula % (e 'Errado' como o complemento)
        hits_list   = [100.0 * hits[k] / considered_proteins for k in ks]
        misses_list = [100.0 * (considered_proteins - hits[k]) / considered_proteins for k in ks]

        # ========= Plot no mesmo estilo do seu script anterior =========
        x = np.arange(len(ks))   # posições no eixo X
        width = 0.25             # barras mais "delicadas"

        plt.figure(figsize=(8, 6))
        bars1 = plt.bar(x - width/2, hits_list,   width, label="CCE", color="#4caf50")
        bars2 = plt.bar(x + width/2, misses_list, width, label="CNCE",  color="#ef5350")

        # Rótulos de % acima das barras
        for bars in (bars1, bars2):
            for bar in bars:
                yval = bar.get_height()
                plt.text(
                    bar.get_x() + bar.get_width()/2, yval + 1,
                    f"{yval:.1f}%",
                    ha="center", va="bottom", fontsize=9
                )

        plt.xticks(x, [f"Top {k}" for k in ks])
        plt.ylabel("Proteínas (%)")
        plt.ylim(0, 105)
        plt.title(label)
        plt.legend(loc="upper left")

        # N total discreto (mesmo lugar do gráfico anterior)
        plt.text(
            0.02, 0.88, f"N = {considered_proteins}",
            transform=plt.gca().transAxes, ha="left", va="top",
            fontsize=9, color="gray"
        )

        out_path = os.path.join(output_dir, outname)
        plt.savefig(out_path)
        plt.close()
        print(f"[INFO] Gráfico cumulativo ({label}) gerado em: {out_path}")


def main(summary_csv):
    plot_violin_descriptors(summary_csv)
    plot_distributions_and_roc(summary_csv)
    plot_cumulative_from_scores(summary_csv)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_summary.py <all_cavities_summary.csv>")
        sys.exit(1)
    summary_csv = sys.argv[1]
    main(summary_csv)
