#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
    """Curvas ROC One-vs-Rest (d, b, n)."""
    sub = df[["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados válidos para {descriptor}.")
        return

    y_score = sub[descriptor].values
    cores = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}
    nomes = {"d": "Drogável", "b": "Borderline", "n": "Não Drogável"}

    plt.figure(figsize=(7, 6))
    for cat in ["d", "b", "n"]:
        y_true = (sub["Categoria"] == cat).astype(int).values
        if len(np.unique(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        N = int(sum(y_true))
        plt.plot(fpr, tpr, lw=2, color=cores[cat],
                 label=f"{nomes[cat]} (AUC={roc_auc:.2f}, N={N})")

    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curvas ROC (Um contra o Resto)\nDescritor: {label}")
    plt.legend(loc="lower right")
    plt.tight_layout()

    out_path = os.path.join(pasta_saida, f"roc_multiclasse_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico salvo: {out_path}")


def plotar_roc_drogavel_borderline_vs_nao(df, descriptor, label, pasta_saida):
    """Curva ROC (D+B) vs N."""
    sub = df[df["Categoria"].isin(["d", "b", "n"])][["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados para (D+B) vs N em {descriptor}.")
        return
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
    out_path = os.path.join(pasta_saida, f"roc_drogavel_borderline_vs_nao_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico salvo: {out_path}")


def plotar_roc_drogavel_vs_nao(df, descriptor, label, pasta_saida):
    """Curva ROC binária pura: Drogável vs Não Drogável."""
    sub = df[df["Categoria"].isin(["d", "n"])][["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados para (D vs N) em {descriptor}.")
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
    out_path = os.path.join(pasta_saida, f"roc_drogavel_vs_nao_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico salvo: {out_path}")


def plotar_densidade(df, descriptor, label, pasta_saida):
    """Gráfico de densidade das distribuições."""
    sub = df[["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Sem dados para densidade de {descriptor}.")
        return
    plt.figure(figsize=(7, 5))
    sns.kdeplot(
        data=sub,
        x=descriptor,
        hue="Categoria",
        fill=True,
        common_norm=False,
        palette={"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"},
        alpha=0.5
    )
    plt.xlabel(label)
    plt.ylabel("Densidade de probabilidade")
    plt.title(f"Distribuição dos valores — {label}")
    plt.tight_layout()
    out_path = os.path.join(pasta_saida, f"densidade_{descriptor}.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Gráfico salvo: {out_path}")


def main(planilha):
    print(f"[INICIANDO] Arquivo: {planilha}")
    if not os.path.exists(planilha):
        print(f"[ERRO] Arquivo não encontrado: {planilha}")
        sys.exit(1)

    # Leitura
    if planilha.endswith(".xlsx"):
        df = pd.read_excel(planilha)
    else:
        df = pd.read_csv(planilha)

    # Detecta coluna kd/ki
    colunas_possiveis = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    col_afinidade = next((c for c in df.columns if any(alt.lower() == c.lower() for alt in colunas_possiveis)), None)
    if col_afinidade is None:
        print("[ERRO] Nenhuma coluna de KD/KI encontrada.")
        print("Colunas disponíveis:", df.columns.tolist())
        sys.exit(1)

    print(f"[INFO] Coluna de afinidade detectada: {col_afinidade}")
    df["Categoria"] = df[col_afinidade].apply(classificar_drogabilidade)
    print("\n[INFO] Distribuição das categorias:")
    print(df["Categoria"].value_counts(dropna=False))

    descritores = {
        "Avg_Fukui": "Fukui médio",
        "Avg_Fukushima": "Fukushima médio",
        "Fukushima_NoNorm": "Fukushima total (soma bruta)",
        "Fukushima_SumNorm": "Fukushima normalizado (átomos ativos)",
    }

    pasta_saida = "output_druggability_final"
    os.makedirs(pasta_saida, exist_ok=True)

    for desc, label in descritores.items():
        if desc not in df.columns:
            print(f"[AVISO] Coluna {desc} não encontrada no arquivo.")
            continue
        print(f"\n[PROCESSANDO] {desc} ({label})")
        plotar_roc_multiclasse(df, desc, label, pasta_saida)
        plotar_roc_drogavel_borderline_vs_nao(df, desc, label, pasta_saida)
        plotar_roc_drogavel_vs_nao(df, desc, label, pasta_saida)
        plotar_densidade(df, desc, label, pasta_saida)

    print("\n[CONCLUÍDO] Todos os gráficos foram gerados em:", pasta_saida)


if __name__ == "__main__":
    print("[DEBUG] Entrando no bloco principal.")
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_final.py <planilha_pdbbind.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
