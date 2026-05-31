#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================
INPUT_FILE = 'tf_snv_analysis_results/tf_snv_tissue_pioneer_class_family.tsv'
OUTPUT_TISSUES = 'tissues_summary_filtered_2.tsv'
OUTPUT_FILTERED_LIST = 'filtered_tf_tissue_list_v2.tsv'
OUTPUT_ENRICHMENT = 'table2_hypergeom_enrichment_v2.tsv'
MIN_THRESHOLD = 1000  # порог для суммарных K(tf+t) и N(tf+t) по ткани
# =============================================

print("1. Загрузка исходных данных...")
df = pd.read_csv(INPUT_FILE, sep='\t')

# чистим названия столбцов
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

# переводим числовые столбцы
df['N(tf+t)'] = pd.to_numeric(df['N(tf+t)'], errors='coerce')
df['K(tf+t)'] = pd.to_numeric(df['K(tf+t)'], errors='coerce')

# убираем строки с пропусками в критичных столбцах
df = df.dropna(subset=['TF', 'tissue', 'N(tf+t)', 'K(tf+t)']).copy()

print(f"   Строк: {len(df)} | Уникальных тканей: {df['tissue'].nunique()} | Уникальных TF: {df['TF'].nunique()}")

# 2. Агрегация по тканям
print("2. Суммирование K(tf+t) и N(tf+t) по тканям...")
tissue_agg = df.groupby('tissue', as_index=False).agg(
    K_sum=('K(tf+t)', 'sum'),
    N_sum=('N(tf+t)', 'sum')
)

# 3. Фильтрация тканей по порогу
print(f"3. Фильтрация тканей (K_sum >= {MIN_THRESHOLD}, N_sum >= {MIN_THRESHOLD})...")
valid_tissues = tissue_agg[
    (tissue_agg['K_sum'] >= MIN_THRESHOLD) &
    (tissue_agg['N_sum'] >= MIN_THRESHOLD)
].copy()

if len(valid_tissues) == 0:
    raise ValueError("Ни одна ткань не прошла фильтрацию. Уменьшите MIN_THRESHOLD.")

print(f"   Отобрано тканей: {len(valid_tissues)} из {len(tissue_agg)}")

valid_tissues.to_csv(OUTPUT_TISSUES, sep='\t', index=False)
print(f"   Таблица тканей сохранена: {OUTPUT_TISSUES}")

# 4. Фильтрация исходной таблицы по отобранным тканям
print("4. Фильтрация исходной таблицы по отобранным тканям...")
valid_tissue_names = valid_tissues['tissue'].tolist()
df_filtered = df[df['tissue'].isin(valid_tissue_names)].copy()

# пересчитываем K/N
df_filtered['prop_K_N'] = np.where(
    df_filtered['N(tf+t)'] > 0,
    df_filtered['K(tf+t)'] / df_filtered['N(tf+t)'],
    0.0
)

df_filtered.to_csv(OUTPUT_FILTERED_LIST, sep='\t', index=False)
print(f"   Отфильтрованный список TF×tissue сохранён: {OUTPUT_FILTERED_LIST}")
print(f"   Строк после фильтрации: {len(df_filtered)}")

# 5. Отбор топ-5% комбинаций по K/N
print("5. Отбор топ-5% комбинаций по prop_K_N...")
threshold_95 = np.percentile(df_filtered['prop_K_N'], 95)
df_top5 = df_filtered[df_filtered['prop_K_N'] >= threshold_95].copy()

print(f"   Порог 95-го перцентиля: {threshold_95:.6f}")
print(f"   Строк в топ-5%: {len(df_top5)}")

if len(df_top5) == 0:
    raise ValueError("В топ-5% не осталось данных.")

# 6. Агрегация по TF
print("6. Расчёт обогащения по TF...")
tf_agg = df_top5.groupby('TF', as_index=False).agg(
    a=('K(tf+t)', 'sum'),
    total_sites=('N(tf+t)', 'sum')
)

tf_agg['b'] = tf_agg['total_sites'] - tf_agg['a']
tf_agg = tf_agg.drop(columns=['total_sites'])

# справочная таблица с метаданными TF
tf_meta = df_top5[[
    'TF',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family'
]].drop_duplicates(subset=['TF']).copy()

# если вдруг для одного TF метаданные различаются, можно проверить это отдельно
meta_check = df_top5.groupby('TF').agg(
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
    print("ВНИМАНИЕ: у некоторых TF метаданные различаются между строками:")
    print(bad_meta.to_string(index=False))

# объединяем статистику с метаданными
tf_agg = tf_agg.merge(tf_meta, on='TF', how='left')

# глобальные суммы по топ-5%
total_k_bg = tf_agg['a'].sum()
total_b_bg = tf_agg['b'].sum()
total_sites_bg = total_k_bg + total_b_bg

tf_agg['c'] = total_k_bg - tf_agg['a']
tf_agg['d'] = total_b_bg - tf_agg['b']

# гипергеометрический тест
tf_agg['p_value'] = tf_agg.apply(
    lambda row: hypergeom.sf(
        row['a'] - 1,
        total_sites_bg,
        total_k_bg,
        row['a'] + row['b']
    ),
    axis=1
)

# fold enrichment
obs_freq = tf_agg['a'] / (tf_agg['a'] + tf_agg['b'])
bg_freq = total_k_bg / total_sites_bg if total_sites_bg > 0 else np.nan
tf_agg['fold_enrichment'] = np.where(
    (tf_agg['a'] + tf_agg['b'] > 0) & (bg_freq > 0),
    obs_freq / bg_freq,
    0.0
)

# FDR Benjamini-Hochberg
def fdr_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    if m == 0:
        return pvals

    order = np.argsort(pvals)
    ranked_pvals = pvals[order]

    qvals = np.empty(m, dtype=float)
    qvals[-1] = min(ranked_pvals[-1], 1.0)

    for i in range(m - 2, -1, -1):
        qvals[i] = min(ranked_pvals[i] * m / (i + 1), qvals[i + 1])

    qvals = np.minimum(qvals, 1.0)

    result = np.empty(m, dtype=float)
    result[order] = qvals
    return result

tf_agg['p_value_fdr'] = fdr_bh(tf_agg['p_value'].values)

# 7. Сортировка и сохранение
print("7. Сортировка и сохранение результата...")
tf_agg = tf_agg.sort_values('p_value_fdr', ascending=True).reset_index(drop=True)

out_cols = [
    'TF',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'a',
    'b',
    'c',
    'd',
    'p_value',
    'p_value_fdr',
    'fold_enrichment'
]

tf_agg[out_cols].to_csv(OUTPUT_ENRICHMENT, sep='\t', index=False)

print(f"\nГотово! Результат сохранён в {OUTPUT_ENRICHMENT}")
print(f"ТФ в итоговой таблице: {len(tf_agg)}")
print(f"Значимых при FDR < 0.05: {(tf_agg['p_value_fdr'] < 0.05).sum()}")
print(f"Медианный fold_enrichment: {tf_agg['fold_enrichment'].median():.3f}")