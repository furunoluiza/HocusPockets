import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc


def classificar_drogabilidade(kd_value):
    """Classifica o alvo como Drogável (d), Borderline (b) ou Não Drogável (n) com base em KD/KI (µM)."""
    if pd.isna(kd_value):
        return np.nan
    try:
        kd_value = float(kd_value)
    except ValueError:
        return np.nan
    if kd_value > 10000:  # conversão nM → µM
        kd_value /= 1000.0
    if kd_value < 0.3:
        return "d"
    elif kd_value <= 1.0:
        return "b"
    else:
        return "n"


def plotar_roc_multiclasse(df, descriptor, label, pasta_saida):
    """Curvas ROC (One-vs-Rest) para d, b, n."""
    sub = df[["Categoria", descriptor]].dropna()
    if sub["Categoria"].nunique() < 2:
        return

    y_score = sub[descriptor].values
    cores = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}
    nomes = {"d": "Drogável", "b": "Borderline", "n": "Não Drogável"}
    aucs = {}

    plt.figure(figsize=(7, 6))
    for cat in ["d", "b", "n"]:
        y_true = (sub["Categoria"] == cat).astype(int).values
        if len(np.unique(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        aucs[cat] = roc_auc
        N = int(sum(y_true))
        plt.plot(fpr, tpr, lw=2, color=cores[cat],
                 label=f"{nomes[cat]} (AUC={roc_auc:.2f}, N={N})")

    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curvas ROC (Um contra o Resto)\nDescritor: {label}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_saida, f"roc_multiclasse_{descriptor}.png"), dpi=300)
    plt.close()


def plotar_roc_drogavel_borderline_vs_nao(df, descriptor, label, pasta_saida):
    """Curva ROC: (Drogável + Borderline) vs Não Drogável."""
    sub = df[df["Categoria"].isin(["d", "b", "n"])][["Categoria", descriptor]].dropna()
    y_true = sub["Categoria"].isin(["d", "b"]).astype(int).values
    y_score = sub[descriptor].values
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    N = int(sum(y_true))

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, lw=2, color="#2196f3", label=f"Drogável + Borderline (AUC={roc_auc:.2f}, N={N})")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curva ROC — Drogável + Borderline vs Não Drogável\nDescritor: {label}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_saida, f"roc_drogavel_borderline_vs_nao_{descriptor}.png"), dpi=300)
    plt.close()
    print(f"[INFO] ROC (D+B vs N) — {descriptor}: AUC={roc_auc:.3f}, N={N}")


def plotar_roc_drogavel_vs_nao(df, descriptor, label, pasta_saida):
    """Curva ROC binária pura: Drogável vs Não Drogável (eliminação dos borderline)."""
    sub = df[df["Categoria"].isin(["d", "n"])][["Categoria", descriptor]].dropna()
    if sub.empty:
        return
    y_true = (sub["Categoria"] == "d").astype(int).values
    y_score = sub[descriptor].values
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    N = int(sum(y_true))

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, lw=2, color="#4caf50", label=f"Drogável vs Não Drogável (AUC={roc_auc:.2f}, N={N})")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curva ROC — Drogável vs Não Drogável\nDescritor: {label}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_saida, f"roc_drogavel_vs_nao_{descriptor}.png"), dpi=300)
    plt.close()
    print(f"[INFO] ROC (D vs N) — {descriptor}: AUC={roc_auc:.3f}, N={N}")


def plotar_densidade(df, descriptor, label, pasta_saida):
    """Gráfico de densidade (distribuição dos valores do descritor por categoria)."""
    sub = df[["Categoria", descriptor]].dropna()
    plt.figure(figsize=(7, 5))
    sns.kdeplot(data=sub, x=descriptor, hue="Categoria",
                fill=True, common_norm=False,
                palette={"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"},
                alpha=0.5)
    plt.xlabel(label)
    plt.ylabel("Densidade de probabilidade")
    plt.title(f"Distribuição dos valores — {label}")
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_saida, f"densidade_{descriptor}.png"), dpi=300)
    plt.close()


def main(planilha):
    # --- Leitura ---
    if planilha.endswith(".xlsx"):
        df = pd.read_excel(planilha)
    else:
        df = pd.read_csv(planilha)

    # --- Detecta coluna KD/KI (µM ou uM) ---
    colunas_possiveis = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    col_afinidade = next((c for c in df.columns if any(alt.lower() == c.lower() for alt in colunas_possiveis)), None)
    if col_afinidade is None:
        print("[ERRO] Nenhuma coluna de KD/KI encontrada. Colunas disponíveis:", df.columns.tolist())
        sys.exit(1)

    print(f"[INFO] Coluna de afinidade detectada: {col_afinidade}")
    df["Categoria"] = df[col_afinidade].apply(classificar_drogabilidade)
    print("\n[INFO] Distribuição das categorias:")
    print(df["Categoria"].value_counts(dropna=False))

    # --- Descritores ---
    DESCRITORES = {
        "Avg_Fukui": "Fukui médio",
        "Avg_Fukushima": "Fukushima médio",
        "Fukushima_NoNorm": "Fukushima total (soma bruta)",
        "Fukushima_SumNorm": "Fukushima normalizado (átomos ativos)",
    }

    pasta_saida = "output_druggability_final"
    os.makedirs(pasta_saida, exist_ok=True)

    for desc, label in DESCRITORES.items():
        if desc not in df.columns:
            print(f"[AVISO] Coluna {desc} não encontrada.")
            continue

        # 1️⃣ ROC One-vs-Rest (3 curvas)
        plotar_roc_multiclasse(df, desc, label, pasta_saida)
        # 2️⃣ ROC Drogável + Borderline vs Não Drogável
        plotar_roc_drogavel_borderline_vs_nao(df, desc, label, pasta_saida)
        # 3️⃣ ROC Drogável vs Não Drogável
        plotar_roc_drogavel_vs_nao(df, desc, label, pasta_saida)
        # 4️⃣ Gráfico de densidade
        plotar_densidade(df, desc, label, pasta_saida)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_final.py <planilha_pdbbind.xlsx>")
        sys.exit
