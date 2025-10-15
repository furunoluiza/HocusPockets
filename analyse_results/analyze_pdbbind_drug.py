#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc


# ==========================================================
# CLASSIFICAÇÃO DOS ALVOS PELA AFINIDADE KD/KI
# ==========================================================
def classificar_drogabilidade(kd_value):
    """Classifica o alvo como Drogável (d), Borderline (b) ou Não Drogável (n)."""
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


# ==========================================================
# FUNÇÃO DE PLOTAGEM DA CURVA ROC BINÁRIA
# ==========================================================
def plotar_roc_binaria(sub, classe_pos, classe_neg, descriptor, label, pasta_saida, nome_fig):
    """Gera uma curva ROC comparando duas classes específicas."""
    if sub.empty:
        print(f"[AVISO] Nenhum dado disponível para {nome_fig}.")
        return

    # Define as classes de interesse
    sub = sub[sub["Categoria"].isin([classe_pos, classe_neg])].dropna()
    if sub.empty:
        print(f"[AVISO] Nenhum dado válido para {classe_pos} vs {classe_neg} em {descriptor}.")
        return

    y_true = (sub["Categoria"] == classe_pos).astype(int).values
    y_score = sub[descriptor].values

    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    N_pos = int(sum(y_true))
    N_neg = int(len(y_true) - N_pos)

    # Mapeia as cores e nomes
    cores = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}
    nomes = {"d": "Drogável", "b": "Borderline", "n": "Não Drogável"}

    plt.figure(figsize=(7, 6))
    plt.plot(fpr, tpr, lw=2, color=cores[classe_pos],
             label=f"{nomes[classe_pos]} (AUC={roc_auc:.2f}, N={N_pos})")
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--")
    plt.xlabel("Taxa de Falsos Positivos")
    plt.ylabel("Taxa de Verdadeiros Positivos")
    plt.title(f"Curva ROC — {nomes[classe_pos]} vs {nomes[classe_neg]}\nDescritor: {label}")
    plt.legend(loc="lower right", fontsize=10)
    plt.tight_layout()

    caminho = os.path.join(pasta_saida, f"roc_{nome_fig}_{descriptor}.png")
    plt.savefig(caminho, dpi=300)
    plt.close()

    print(f"[OK] Curva ROC salva: {caminho}")
    print(f"     AUC={roc_auc:.3f} | N({nomes[classe_pos]})={N_pos} | N({nomes[classe_neg]})={N_neg}")


# ==========================================================
# FUNÇÃO DE PLOTAGEM DA DENSIDADE (KDE)
# ==========================================================
def plotar_densidade(df, descriptor, label, pasta_saida):
    """Gera o gráfico de densidade (sobreposição das distribuições por categoria)."""
    sub = df[["Categoria", descriptor]].dropna()
    if sub.empty:
        print(f"[AVISO] Nenhum dado para densidade de {descriptor}.")
        return

    plt.figure(figsize=(8, 5))
    sns.kdeplot(
        data=sub,
        x=descriptor,
        hue="Categoria",
        fill=True,
        common_norm=False,
        palette={"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"},
        alpha=0.4,
    )

    # Adiciona linhas verticais nas médias
    medias = sub.groupby("Categoria")[descriptor].mean()
    for cat, media in medias.items():
        cor = {"d": "#4caf50", "b": "#ffa726", "n": "#ef5350"}[cat]
        plt.axvline(media, color=cor, linestyle="--", lw=1.5)

    nomes_legenda = {
        "d": "Drogável",
        "b": "Borderline",
        "n": "Não Drogável"
    }
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(handles, [nomes_legenda.get(lbl, lbl) for lbl in labels], fontsize=10)

    plt.xlabel(label)
    plt.ylabel("Densidade de probabilidade")
    plt.title(f"Distribuição dos valores do descritor — {label}")
    plt.tight_layout()
    caminho = os.path.join(pasta_saida, f"densidade_{descriptor}.png")
    plt.savefig(caminho, dpi=300)
    plt.close()
    print(f"[OK] Gráfico de densidade salvo: {caminho}")


# ==========================================================
# FUNÇÃO PRINCIPAL
# ==========================================================
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

    pasta_saida = "output_druggability_comparacoes"
    os.makedirs(pasta_saida, exist_ok=True)

    for desc, label in descritores.items():
        if desc not in df.columns:
            print(f"[AVISO] Coluna {desc} não encontrada.")
            continue

        print(f"\n[PROCESSANDO] {desc} ({label})")

        # 1️⃣ Densidade
        plotar_densidade(df, desc, label, pasta_saida)

        # 2️⃣ Drogável vs Borderline
        plotar_roc_binaria(df, "d", "b", desc, label, pasta_saida, "drogavel_vs_borderline")

        # 3️⃣ Drogável + Borderline vs Não Drogável
        df_db = df.copy()
        df_db["Categoria"] = df_db["Categoria"].replace({"b": "d"})
        plotar_roc_binaria(df_db, "d", "n", desc, label, pasta_saida, "drogavel_borderline_vs_nao")

        # 4️⃣ Drogável vs Não Drogável (eliminação de borderline)
        df_puro = df[df["Categoria"].isin(["d", "n"])]
        plotar_roc_binaria(df_puro, "d", "n", desc, label, pasta_saida, "drogavel_vs_nao")

    print("\n[CONCLUÍDO] Todos os gráficos foram gerados em:", pasta_saida)


# ==========================================================
# EXECUÇÃO PRINCIPAL
# ==========================================================
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_druggability_comparacoes.py <planilha_pdbbind.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
