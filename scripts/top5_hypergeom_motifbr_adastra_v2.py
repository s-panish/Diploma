#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================

INPUT_FILE = "motifbreakr_adastra_with_cell_lines.tsv"
OUTPUT_FILE = "tf_top5_frequency_motifbr_adastra.tsv"
TOP5_PAIRS_FILE = "top5_tf_cell_lines_motifbr.tsv"
TF_CELL_AGG_FILE = "tf_cell_line_aggregated_motifbr.tsv"

# Минимальные пороги для включения пары TF × cell_line в анализ.
# Для новой таблицы на ~8931 строк K >= 100 может быть слишком строгим.
MIN_K = 3
MIN_k = 1

# Колонки, по которым определяется изменение мотива
SCORE_COLS = [
    "comparison_motif_score_0.001",
    "comparison_motif_score_0.0005",
    "comparison_motif_score_0.0001"
]

def fdr_bh(pvals):
    """
    Коррекция p-value методом Benjamini-Hochberg.
    """
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)

    if m == 0:
        return pvals

    order = np.argsort(pvals)
    sorted_pvals = pvals[order]

    qvals_sorted = np.empty(m, dtype=float)
    qvals_sorted[-1] = sorted_pvals[-1]

    for i in range(m - 2, -1, -1):
        qvals_sorted[i] = min(
            sorted_pvals[i] * m / (i + 1),
            qvals_sorted[i + 1]
        )

    qvals = np.empty(m, dtype=float)
    qvals[order] = qvals_sorted

    return np.minimum(qvals, 1.0)


print("1. Загрузка данных...")
df = pd.read_csv(INPUT_FILE, sep="\t", low_memory=False)

print(f"   Всего строк во входной таблице: {len(df):,}")

required_cols = ["SNP_id", "providerId", "cell_line"] + SCORE_COLS
missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    raise ValueError(
        "Во входной таблице отсутствуют необходимые колонки: "
        + ", ".join(missing_cols)
    )

print(f"   Уникальных SNP_id: {df['SNP_id'].nunique():,}")
print(f"   Уникальных providerId: {df['providerId'].nunique():,}")
print(f"   Уникальных cell_line: {df['cell_line'].nunique():,}")


# 2. Подготовка score-колонок
print("2. Определение вариантов, изменяющих мотив...")

for col in SCORE_COLS:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# -1 = ослабление/потеря мотива
#  1 = усиление/появление мотива
#  0 = нет значимого изменения по данному порогу
df["is_broken"] = df[SCORE_COLS].abs().eq(1).any(axis=1)

print(f"   Строк с признаком изменения мотива: {df['is_broken'].sum():,}")

# 3. Извлечение базового имени TF из providerId
print("3. Извлечение имени TF из providerId...")

df["TF"] = (
    df["providerId"]
    .astype(str)
    .str.split(".")
    .str[0]
)

print(f"   Уникальных TF после извлечения имени: {df['TF'].nunique():,}")


# 4. Удаление повторов, возникающих из-за структуры объединённой таблицы
print("4. Удаление технических повторов...")

before_dedup = len(df)

# Основная единица анализа:
# один вариант × один мотив × одна клеточная линия
dedup_cols = ["SNP_id", "providerId", "TF", "cell_line"]

df_unique = (
    df[dedup_cols + ["is_broken"]]
    .drop_duplicates()
    .copy()
)

after_dedup = len(df_unique)

print(f"   Было строк: {before_dedup:,}")
print(f"   Осталось уникальных SNP_id × providerId × cell_line: {after_dedup:,}")


# 5. Агрегация по каждой паре TF × cell_line
print("5. Расчёт k, K и ratio для каждой пары TF × cell_line...")

tf_cell_agg = (
    df_unique
    .groupby(["TF", "cell_line"], as_index=False)
    .agg(
        k=("is_broken", "sum"),
        K=("is_broken", "size"),
        unique_snps=("SNP_id", "nunique"),
        unique_motifs=("providerId", "nunique")
    )
)

tf_cell_agg["k"] = tf_cell_agg["k"].astype(int)
tf_cell_agg["K"] = tf_cell_agg["K"].astype(int)

print(f"   Пар TF × cell_line до фильтрации: {len(tf_cell_agg):,}")
print(f"   Уникальных cell_line до фильтрации: {tf_cell_agg['cell_line'].nunique():,}")
print(f"   Уникальных TF до фильтрации: {tf_cell_agg['TF'].nunique():,}")


# 6. Фильтрация малых групп
print(f"6. Применение фильтра: K >= {MIN_K} и k >= {MIN_k}...")

tf_cell_agg = tf_cell_agg[
    (tf_cell_agg["K"] >= MIN_K) &
    (tf_cell_agg["k"] >= MIN_k)
].copy()

print(f"   Пар TF × cell_line после фильтрации: {len(tf_cell_agg):,}")
print(f"   Уникальных cell_line после фильтрации: {tf_cell_agg['cell_line'].nunique():,}")
print(f"   Уникальных TF после фильтрации: {tf_cell_agg['TF'].nunique():,}")

if len(tf_cell_agg) == 0:
    raise ValueError(
        "После фильтрации не осталось данных. "
        "Уменьшите MIN_K и/или MIN_k."
    )

tf_cell_agg["ratio"] = tf_cell_agg["k"] / tf_cell_agg["K"]

tf_cell_agg = tf_cell_agg.sort_values(
    ["ratio", "k", "K"],
    ascending=[False, False, False]
).reset_index(drop=True)

tf_cell_agg.to_csv(TF_CELL_AGG_FILE, sep="\t", index=False)

print(f"   Агрегированная таблица TF × cell_line сохранена в: {TF_CELL_AGG_FILE}")


# 7. Определение глобального топ-5% по ratio
print("7. Расчёт глобального топ-5% по ratio = k / K...")

threshold_95 = np.percentile(tf_cell_agg["ratio"], 95)

tf_cell_agg["in_top5"] = tf_cell_agg["ratio"] >= threshold_95

n = int(tf_cell_agg["in_top5"].sum())  # число пар TF × cell_line в топ-5%
N = int(len(tf_cell_agg))              # общее число пар TF × cell_line после фильтрации

print(f"   Порог 95-го перцентиля: {threshold_95:.6f}")
print(f"   n, число пар в топ-5%: {n:,}")
print(f"   N, всего пар после фильтрации: {N:,}")


# 8. Сохранение всех пар TF × cell_line, попавших в топ-5%
print("8. Сохранение пар TF × cell_line, попавших в топ-5%...")

top5_percent_pairs = (
    tf_cell_agg[tf_cell_agg["in_top5"]]
    .sort_values(["ratio", "k", "K"], ascending=[False, False, False])
    .reset_index(drop=True)
)

top5_cols = [
    "TF",
    "cell_line",
    "k",
    "K",
    "ratio",
    "unique_snps",
    "unique_motifs",
    "in_top5"
]

top5_percent_pairs[top5_cols].to_csv(
    TOP5_PAIRS_FILE,
    sep="\t",
    index=False
)

print(f"   Сохранено {len(top5_percent_pairs):,} пар в: {TOP5_PAIRS_FILE}")
print("\n   Первые 10 пар из топ-5%:")
print(
    top5_percent_pairs[
        ["TF", "cell_line", "k", "K", "ratio"]
    ]
    .head(10)
    .to_string(index=False)
)


# 9. Агрегация по TF
print("9. Подсчёт частоты попадания каждого TF в топ-5%...")

tf_agg = (
    tf_cell_agg
    .groupby("TF", as_index=False)
    .agg(
        k=("in_top5", "sum"),
        K=("in_top5", "size")
    )
)

tf_agg["k"] = tf_agg["k"].astype(int)
tf_agg["K"] = tf_agg["K"].astype(int)

# Здесь:
# k = число cell_line, где TF попал в топ-5%
# K = общее число cell_line, где TF прошёл фильтр и участвовал в анализе
# n = общее число пар TF × cell_line в топ-5%
# N = общее число пар TF × cell_line после фильтрации

tf_agg["n"] = n
tf_agg["N"] = N


# 10. Гипергеометрический тест
print("10. Расчёт p-value гипергеометрическим тестом...")

tf_agg["p_value"] = tf_agg.apply(
    lambda row: hypergeom.sf(
        row["k"] - 1,
        row["N"],
        row["n"],
        row["K"]
    ),
    axis=1
)


# 11. Fold enrichment
print("11. Расчёт fold enrichment...")

background_frequency = n / N

tf_agg["observed_frequency"] = tf_agg["k"] / tf_agg["K"]
tf_agg["background_frequency"] = background_frequency

tf_agg["fold_enrichment"] = np.where(
    tf_agg["background_frequency"] > 0,
    tf_agg["observed_frequency"] / tf_agg["background_frequency"],
    np.nan
)


# 12. FDR-коррекция
print("12. FDR-коррекция Benjamini-Hochberg...")

tf_agg["p_value_fdr"] = fdr_bh(tf_agg["p_value"].values)


# 13. Сортировка и сохранение основной таблицы
print("13. Сортировка и сохранение итоговой таблицы...")

tf_agg = tf_agg.sort_values(
    ["p_value_fdr", "p_value", "fold_enrichment"],
    ascending=[True, True, False]
).reset_index(drop=True)

out_cols = [
    "TF",
    "k",
    "K",
    "n",
    "N",
    "observed_frequency",
    "background_frequency",
    "p_value",
    "p_value_fdr",
    "fold_enrichment"
]

tf_agg[out_cols].to_csv(OUTPUT_FILE, sep="\t", index=False)


print("\nГотово.")
print(f"Итоговая таблица сохранена в: {OUTPUT_FILE}")
print(f"Таблица пар TF × cell_line сохранена в: {TF_CELL_AGG_FILE}")
print(f"Таблица топ-5% пар сохранена в: {TOP5_PAIRS_FILE}")
print(f"TF в итоговой таблице: {len(tf_agg):,}")
print(f"Значимых TF при FDR < 0.05: {(tf_agg['p_value_fdr'] < 0.05).sum():,}")
print(f"Медианный fold_enrichment: {tf_agg['fold_enrichment'].median():.3f}")