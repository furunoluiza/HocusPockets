#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc


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
# GRÁFICO ROC MULTICLASSE
# ============================================================
def plotar_roc_multiclasse(df, descriptor, label, pasta_saida):
    """Plota as 3 curvas ROC (d, b, n vs resto) no mesmo gráfico."""
    sub = df[["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados válidos para {descriptor}.")
        return

    cores = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}
    nomes = {"d": "Drogável", "b": "Borderline", "n": "Não Drogável"}

    plt.figure(figsize=(7, 6))
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
        plt.plot(fpr, tpr, lw=2, color=cores[cat],
                 label=f"{nomes[cat]} (AUC={roc_auc:.2f}, N+={n_pos}, N-={n_neg})")

    plt.plot([0, 1], [0, 1], color="gray", lw=1)
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curvas ROC — Multiclasse (Um contra o Resto)\nDescritor: {label}")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    out_path = os.path.join(pasta_saida, f"roc_multiclasse_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico multiclasse salvo: {out_path}")


# ============================================================
# GRÁFICO ROC BINÁRIO (D vs N)
# ============================================================
def plotar_roc_drogavel_vs_nao(df, descriptor, label, pasta_saida):
    """Plota ROC binária: Drogável vs Não Drogável (sem borderline)."""
    sub = df[df["Categoria"].isin(["d", "n"])][["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados válidos para D vs N em {descriptor}.")
        return

    y_true = (sub["Categoria"] == "d").astype(int)
    y_score = sub[descriptor]
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    n_pos = int(sum(y_true))
    n_neg = int(len(y_true) - n_pos)

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, lw=2, color="#2e7d32",
             label=f"Drogável vs Não Drogável (AUC={roc_auc:.2f}, N+={n_pos}, N-={n_neg})")
    plt.plot([0, 1], [0, 1], color="gray", lw=1)
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curva ROC — Drogável vs Não Drogável\nDescritor: {label}")
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    out_path = os.path.join(pasta_saida, f"roc_drogavel_vs_nao_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico D vs N salvo: {out_path}")


# ============================================================
# GRÁFICO DE DENSIDADE
# ============================================================
def plotar_densidade(df, descriptor, label, pasta_saida, binario=False):
    """Plota gráfico de densidade (sem preenchimento), com legenda e N."""
    sub = df[["Categoria", descriptor]].dropna()
    if binario:
        sub = sub[sub["Categoria"].isin(["d", "n"])]

    categorias = {
        "d": "Drogável",
        "b": "Borderline",
        "n": "Não Drogável"
    }
    cores = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}

    plt.figure(figsize=(7, 5))
    for cat, nome in categorias.items():
        subset = sub[sub["Categoria"] == cat][descriptor]
        if len(subset) == 0:
            continue
        sns.kdeplot(subset, color=cores[cat], lw=2, fill=False,
                    label=f"{nome} (N={len(subset)})")
    plt.xlabel(label)
    plt.ylabel("Densidade de probabilidade")
    titulo = f"Distribuição dos valores — {label}"
    if binario:
        titulo += "\n(Drogável vs Não Drogável)"
    plt.title(titulo)
    plt.legend()
    plt.tight_layout()

    suffix = "_binary" if binario else "_multiclasse"
    out_path = os.path.join(pasta_saida, f"densidade{suffix}_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico de densidade salvo: {out_path}")


# ============================================================
# MAIN
# ============================================================
def main(planilha):
    print(f"[INICIANDO] Lendo arquivo: {planilha}")
    if not os.path.exists(planilha):
        print(f"[ERRO] Arquivo não encontrado: {planilha}")
        sys.exit(1)

    # Lê planilha
    if planilha.endswith(".xlsx"):
        df = pd.read_excel(planilha)
    else:
        df = pd.read_csv(planilha)

    # Detecta coluna KD/KI
    colunas_possiveis = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    col_afinidade = next((c for c in df.columns if any(alt.lower() == c.lower() for alt in colunas_possiveis)), None)
    if col_afinidade is None:
        print("[ERRO] Coluna de afinidade não encontrada.")
        print("Colunas disponíveis:", df.columns.tolist())
        sys.exit(1)

    # Classifica categorias
    df["Categoria"] = df[col_afinidade].apply(classificar_drogabilidade)
    print("[INFO] Distribuição de categorias:")
    print(df["Categoria"].value_counts(dropna=False))

    descritores = {
        "Avg_Fukui": "Fukui médio",
        "Avg_Fukushima": "Fukushima médio",
        "Fukushima_NoNorm": "Fukushima total (soma bruta)",
        "Fukushima_SumNorm": "Fukushima normalizado (átomos ativos)"
    }

    pasta_saida = "output_druggability_final"
    os.makedirs(pasta_saida, exist_ok=True)

    for desc, label in descritores.items():
        if desc not in df.columns:
            print(f"[AVISO] Coluna {desc} não encontrada.")
            continue
        print(f"\n[PROCESSANDO] {desc} ({label})")

        # ROC multiclasse
        plotar_roc_multiclasse(df, desc, label, pasta_saida)
        plotar_densidade(df, desc, label, pasta_saida, binario=False)

        # ROC binária (D vs N)
        plotar_roc_drogavel_vs_nao(df, desc, label, pasta_saida)
        plotar_densidade(df, desc, label, pasta_saida, binario=True)

    print("\n✅ [CONCLUÍDO] Todos os gráficos foram gerados em:", pasta_saida)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_final.py <planilha_pdbbind.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
