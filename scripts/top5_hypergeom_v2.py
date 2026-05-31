#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================
INPUT_FILE = 'tf_snv_analysis_results/tf_snv_tissue_pioneer_class_family_filtered.tsv'
OUTPUT_FILE = 'table2_tf_top5_frequency_enrichment_v2.tsv'
TOP5_PAIRS_FILE = 'top5_tf_tissues_v2.tsv'   # файл для топ-5% пар TF×tissue
TISSUE_COL = 'tissue'
# =============================================

print("1. Загрузка данных...")
df = pd.read_csv(INPUT_FILE, sep='\t')
df.columns = [col.strip() for col in df.columns]

required_cols = [
    'TF',
    'tissue',
    'N(tf+t)',
    'K(tf+t)',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family'
]
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    raise ValueError(
        f"В файле отсутствуют обязательные столбцы: {missing_cols}\n"
        f"Найденные столбцы: {list(df.columns)}"
    )

# приведение типов
df['N(tf+t)'] = pd.to_numeric(df['N(tf+t)'], errors='coerce')
df['K(tf+t)'] = pd.to_numeric(df['K(tf+t)'], errors='coerce')

# удаляем строки с пропусками в критичных столбцах
df = df.dropna(subset=['TF', 'tissue', 'N(tf+t)', 'K(tf+t)']).copy()

# пересчитываем долю для надежности
df['prop_K_N'] = np.where(
    df['N(tf+t)'] > 0,
    df['K(tf+t)'] / df['N(tf+t)'],
    0.0
)

print(f"   Всего комбинаций TF×ткань: {len(df)}")
print(f"   Уникальных ТФ: {df['TF'].nunique()}")

# 2. Определение глобального порога топ-5% по prop_K_N
print("2. Расчёт глобального топ-5%...")
threshold_95 = np.percentile(df['prop_K_N'], 95)
df['in_top5'] = df['prop_K_N'] >= threshold_95

n_top = int(df['in_top5'].sum())   # n: всего TF×tissue в топ-5%
N_total = int(len(df))             # N: всего TF×tissue в датасете

print(f"   Порог 95-го перцентиля: {threshold_95:.6f}")
print(f"   n (в топ-5%): {n_top}")
print(f"   N (всего): {N_total}")

# 2.1. Сохранение всех пар TF×tissue, попавших в топ-5%
print("2.1. Сохранение всех пар TF×tissue, попавших в топ-5%...")
top5_percent_pairs = df[df['in_top5']].copy()
top5_percent_pairs = top5_percent_pairs.sort_values('prop_K_N', ascending=False).reset_index(drop=True)

save_cols = [
    'TF',
    'tissue',
    'N(tf+t)',
    'K(tf+t)',
    'prop_K_N',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'in_top5'
]
top5_percent_pairs[save_cols].to_csv(TOP5_PAIRS_FILE, sep='\t', index=False)

print(f"   Сохранено {len(top5_percent_pairs)} пар в {TOP5_PAIRS_FILE}")
print(top5_percent_pairs[['TF', 'tissue', 'prop_K_N']].head(10).to_string(index=False))

# 3. Агрегация по каждому ТФ
print("3. Подсчёт k и K для каждого ТФ...")
tf_agg = df.groupby('TF')['in_top5'].agg(
    k='sum',          # сколько раз TF попал в топ-5%
    total_obs='count' # общее число наблюдений для TF
).reset_index()

tf_agg['K'] = tf_agg['total_obs'] - tf_agg['k']  # сколько раз TF НЕ попал в топ-5%
tf_agg.drop(columns=['total_obs'], inplace=True)

# подтягиваем метаданные TF
tf_meta = df[[
    'TF',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family'
]].drop_duplicates(subset=['TF']).copy()

# проверка на случай, если у одного TF метаданные почему-то различаются
meta_check = df.groupby('TF').agg(
    n_is_pioneer=('is_pioneer', 'nunique'),
    n_class=('tfclass:class', 'nunique'),
    n_family=('tfclass:family', 'nunique')
).reset_index()

bad_meta = meta_check[
    (meta_check['n_is_pioneer'] > 1) |
    (meta_check['n_class'] > 1) |
    (meta_check['n_family'] > 1)
]

if len(bad_meta) > 0:
    print("ВНИМАНИЕ: у некоторых TF различаются метаданные:")
    print(bad_meta.to_string(index=False))

tf_agg = tf_agg.merge(tf_meta, on='TF', how='left')

# добавляем глобальные константы
tf_agg['n'] = n_top
tf_agg['N'] = N_total

# 4. Гипергеометрический тест: P(X >= k)
print("4. Расчёт p-value...")
tf_agg['p_value'] = tf_agg.apply(
    lambda row: hypergeom.sf(
        row['k'] - 1,
        row['N'],
        row['n'],
        row['k'] + row['K']
    ),
    axis=1
)

# 5. Fold enrichment
obs_freq = tf_agg['k'] / (tf_agg['k'] + tf_agg['K'])
bg_freq = n_top / N_total if N_total > 0 else np.nan
tf_agg['fold_enrichment'] = np.where(
    (tf_agg['k'] + tf_agg['K']) > 0,
    obs_freq / bg_freq,
    0.0
)

# 6. FDR-коррекция Benjamini-Hochberg
def fdr_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    if m == 0:
        return pvals

    order = np.argsort(pvals)
    sorted_pvals = pvals[order]

    qvals = np.empty(m, dtype=float)
    qvals[-1] = min(sorted_pvals[-1], 1.0)

    for i in range(m - 2, -1, -1):
        qvals[i] = min(sorted_pvals[i] * m / (i + 1), qvals[i + 1])

    qvals = np.minimum(qvals, 1.0)

    result = np.empty(m, dtype=float)
    result[order] = qvals
    return result

tf_agg['p_value_fdr'] = fdr_bh(tf_agg['p_value'].values)

# 7. Сортировка и сохранение основной таблицы
print("5. Сортировка и сохранение основной таблицы...")
tf_agg = tf_agg.sort_values('p_value_fdr', ascending=True).reset_index(drop=True)

out_cols = [
    'TF',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'k',
    'K',
    'n',
    'N',
    'p_value',
    'p_value_fdr',
    'fold_enrichment'
]

tf_agg[out_cols].to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f"\nГотово! Результат сохранён в {OUTPUT_FILE}")
print(f"ТФ в таблице: {len(tf_agg)}")
print(f"Значимых при FDR < 0.05: {(tf_agg['p_value_fdr'] < 0.05).sum()}")
print(f"Медианный fold_enrichment: {tf_agg['fold_enrichment'].median():.3f}")