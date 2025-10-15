import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc


def classificar_drogabilidade(kd_value):
    """
    Classifica o alvo como Drogável (d), Borderline (b) ou Não Drogável (n)
    com base no valor de KD/KI em µM.
    """
    if pd.isna(kd_value):
        return np.nan

    try:
        kd_value = float(kd_value)
    except ValueError:
        return np.nan

    # Converte nM → µM se necessário
    if kd_value > 10000:
        kd_value = kd_value / 1000.0

    if kd_value < 0.3:
        return "d"
    elif kd_value <= 1.0:
        return "b"
    else:
        return "n"


def plotar_roc_multiclasse(df, descriptor, label, pasta_saida):
    """
    Gera curvas ROC One-vs-Rest para Drogável, Borderline e Não Drogável,
    além de uma curva combinada (Drogável + Borderline) vs Não Drogável.
    """
    os.makedirs(pasta_saida, exist_ok=True)

    sub = df[["Categoria", descriptor]].dropna()
    if sub["Categoria"].nunique() < 2:
        print(f"[AVISO] O descritor {descriptor} não tem variação suficiente entre categorias.")
        return

    y_score = sub[descriptor].values
    classes = ["d", "b", "n"]
    cores = {
        "d": "#4caf50",  # verde
        "b": "#ffa726",  # laranja
        "n": "#ef5350",  # vermelho
        "db": "#2196f3",  # azul (Drogável + Borderline)
    }
    nomes = {
        "d": "Drogável",
        "b": "Borderline",
        "n": "Não Drogável",
        "db": "Drogável + Borderline",
    }

    aucs = {}
    ns = {}

    plt.figure(figsize=(7, 6))

    # --- Curvas individuais One-vs-Rest ---
    for cat in classes:
        y_true = (sub["Categoria"] == cat).astype(int).values
        if len(np.unique(y_true)) < 2:
            continue
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        aucs[cat] = roc_auc
        ns[cat] = int(sum(y_true))
        plt.plot(fpr, tpr, lw=2, color=cores[cat],
                 label=f"{nomes[cat]} (AUC={roc_auc:.2f}, N={ns[cat]})")

    # --- Curva Drogável + Borderline vs Não Drogável ---
    y_true_db = sub["Categoria"].isin(["d", "b"]).astype(int).values
    if len(np.unique(y_true_db)) > 1:
        fpr_db, tpr_db, _ = roc_curve(y_true_db, y_score)
        auc_db = auc(fpr_db, tpr_db)
        aucs["db"] = auc_db
        ns["db"] = int(sum(y_true_db))
        plt.plot(fpr_db, tpr_db, lw=2, color=cores["db"],
                 label=f"{nomes['db']} (AUC={auc_db:.2f}, N={ns['db']})")

    # --- Gráfico base ---
    plt.plot([0, 1], [0, 1], color="gray", linestyle="--", lw=1)
    plt.xlabel("Taxa de Falsos Positivos", fontsize=12)
    plt.ylabel("Taxa de Verdadeiros Positivos", fontsize=12)
    plt.title(f"Curvas ROC (Um contra o Resto)\nDescritor: {label}", fontsize=13)
    plt.legend(loc="lower right", fontsize=9)
    plt.tight_layout()

    caminho_saida = os.path.join(pasta_saida, f"roc_multiclasse_{descriptor}.png")
    plt.savefig(caminho_saida, dpi=300)
    plt.close()

    # --- Impressão dos resultados ---
    if aucs:
        macro_auc = np.mean([v for k, v in aucs.items() if k != "db"])
        print(f"\n[RESULTADO] {descriptor} ({label})")
        for k in ["d", "b", "n", "db"]:
            if k in aucs:
                print(f" - {nomes[k]:<22} → AUC={aucs[k]:.3f}, N={ns[k]}")
        print(f"   Média (macro-AUC): {macro_auc:.3f}")


def main(arquivo_planilha):
    # Detecta extensão e carrega
    if arquivo_planilha.endswith(".xlsx"):
        df = pd.read_excel(arquivo_planilha)
    else:
        df = pd.read_csv(arquivo_planilha)

    # Detecta coluna de afinidade (µM ou uM)
    colunas_possiveis = ["kd/ki (µM)", "kd/ki (uM)", "KD/Ki (µM)", "KD/Ki (uM)"]
    coluna_afinidade = None
    for col in df.columns:
        if any(alt.lower() == col.lower() for alt in colunas_possiveis):
            coluna_afinidade = col
            break

    if coluna_afinidade is None:
        print("[ERRO] Nenhuma coluna de KD/KI encontrada ('kd/ki (µM)' ou 'kd/ki (uM)').")
        print("Colunas disponíveis:", list(df.columns))
        sys.exit(1)
    else:
        print(f"[INFO] Coluna de afinidade detectada: '{coluna_afinidade}'")

    # Classificação
    df["Categoria"] = df[coluna_afinidade].apply(classificar_drogabilidade)
    print("\n[INFO] Distribuição das categorias:")
    print(df["Categoria"].value_counts(dropna=False))

    # Descritores
    DESCRITORES = {
        "Avg_Fukui": "Fukui médio",
        "Avg_Fukushima": "Fukushima médio",
        "Fukushima_NoNorm": "Fukushima total (soma bruta)",
        "Fukushima_SumNorm": "Fukushima normalizado (átomos ativos)",
    }

    pasta_saida = "output_druggability_multiclasse"
    os.makedirs(pasta_saida, exist_ok=True)

    for desc, label in DESCRITORES.items():
        if desc in df.columns:
            plotar_roc_multiclasse(df, desc, label, pasta_saida)
        else:
            print(f"[AVISO] Coluna {desc} não encontrada na planilha.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Uso: python analyze_pdbbind_drug_multiclasse.py <planilha_pdbbind.xlsx>")
        sys.exit(1)
    main(sys.argv[1])
