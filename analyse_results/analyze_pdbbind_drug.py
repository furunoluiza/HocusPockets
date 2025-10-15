#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc


# =============================
# 1️⃣ CLASSIFICAÇÃO DE ALVOS
# =============================
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


# =============================
# 2️⃣ GRÁFICOS ROC
# =============================
def plotar_roc(df, y_true, y_score, titulo, nome_saida, pasta, cor, label_legenda, n_pos, n_neg):
    """Função genérica para plotar uma curva ROC."""
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(6, 6))
    plt.plot(fpr, tpr, color=cor, lw=2,
             label=f"{label_legenda} (AUC={roc_auc:.2f}, N+={n_pos}, N-={n_neg})")
    plt.plot([0, 1], [0, 1], color="gray", lw=1)
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(titulo)
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(pasta, nome_saida), dpi=300)
    plt.close()
    return roc_auc


# =============================
# 3️⃣ GRÁFICO DE DENSIDADE
# =============================
def plotar_densidade(df, descriptor, label, pasta_saida):
    """Gráfico de densidade por categoria, sem preenchimento."""
    sub = df[["Categoria", descriptor]].dropna()

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
        sns.kdeplot(subset, color=cores[cat], lw=2, label=f"{nome} (N={len(subset)})", fill=False)
    plt.xlabel(label)
    plt.ylabel("Densidade de probabilidade")
    plt.title(f"Distribuição dos valores — {label}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(pasta_saida, f"densidade_{descriptor}.png"), dpi=300)
    plt.close()


# =============================
# 4️⃣ EXECUÇÃO PRINCIPAL
# =============================
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

    # Detecta coluna de KD/KI
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

        sub = df[["Categoria", desc]].dropna()

        # 1️⃣ ROC Multiclasse (d, b, n)
        for cat, nome, cor in [
            ("d", "Drogável", "#4caf50"),
            ("b", "Borderline", "#ffa726"),
            ("n", "Não Drogável", "#ef5350")
        ]:
            y_true = (sub["Categoria"] == cat).astype(int)
            y_score = sub[desc]
            n_pos = int(sum(y_true))
            n_neg = int(len(y_true) - n_pos)
            if n_pos > 1 and n_neg > 1:
                auc_val = plotar_roc(
                    df, y_true, y_score,
                    f"Curva ROC — {nome} vs Resto\nDescritor: {label}",
                    f"roc_{cat}_vs_resto_{desc}.png",
                    pasta_saida, cor, nome, n_pos, n_neg
                )
                print(f"   {nome}: AUC={auc_val:.3f}, N+={n_pos}, N-={n_neg}")

        # 2️⃣ Drogável + Borderline vs Não Drogável
        sub_db = sub[sub["Categoria"].isin(["d", "b", "n"])]
        y_true_db = sub_db["Categoria"].isin(["d", "b"]).astype(int)
        y_score_db = sub_db[desc]
        n_pos_db = int(sum(y_true_db))
        n_neg_db = int(len(y_true_db) - n_pos_db)
        if n_pos_db > 1 and n_neg_db > 1:
            auc_val = plotar_roc(
                df, y_true_db, y_score_db,
                f"Curva ROC — Drogável + Borderline vs Não Drogável\nDescritor: {label}",
                f"roc_drogavel_borderline_vs_nao_{desc}.png",
                pasta_saida, "#2196f3", "Drogável + Borderline", n_pos_db, n_neg_db
            )
            print(f"   Drogável+Borderline: AUC={auc_val:.3f}, N+={n_pos_db}, N-={n_neg_db}")

        # 3️⃣ Drogável vs Não Drogável (eliminação dos borderline)
        sub_dn = sub[sub["Categoria"].isin(["d", "n"])]
        y_true_dn = (sub_dn["Categoria"] == "d").astype(int)
        y_score_dn = sub_dn[desc]
        n_pos_dn = int(sum(y_true_dn))
        n_neg_dn = int(len(y_true_dn) - n_pos_dn)
        if n_pos_dn > 1 and n_neg_dn > 1:
            auc_val = plotar_roc(
                df, y_true_dn, y_score_dn,
                f"Curva ROC — Drogável vs Não Drogável\nDescritor: {label}",
                f"roc_drogavel_vs_nao_{desc}.png",
                pasta_saida, "#2e7d32", "Drogável", n_pos_dn, n_neg_dn
            )
            print(f"   Drogável vs Não Drogável: AUC={auc_val:.3f}, N+={n_pos_dn}, N-={n_neg_dn}")

        # 4️⃣ Gráfico de densidade
        plotar_densidade(df, desc, label, pasta_saida)

    print("\n✅ [CONCLUÍDO] Todos os gráficos foram gerados em:", pasta_saida)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_final.py <planilha_pdbbind.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
