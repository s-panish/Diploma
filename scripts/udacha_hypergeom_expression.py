#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================
INPUT_FILE = 'tf_snv_tissue_pioneer_class_family_filtered_with_expression.tsv'

OUTPUT_FILE = 'table_tf_top5_frequency_enrichment_expression_ge1.tsv'
TOP5_PAIRS_FILE = 'top5_tf_tissues_expression_ge1.tsv'

EXPRESSION_THRESHOLD = 1
# =============================================


print('1. Загрузка данных...')
df = pd.read_csv(INPUT_FILE, sep='\t')
df.columns = [col.strip() for col in df.columns]

required_cols = [
    'TF',
    'tissue',
    'N(tf+t)',
    'K(tf+t)',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'expression_value'
]

missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    raise ValueError(
        f'В файле отсутствуют обязательные столбцы: {missing_cols}\n'
        f'Найденные столбцы: {list(df.columns)}'
    )

# Приведение числовых столбцов к числовому формату
df['N(tf+t)'] = pd.to_numeric(df['N(tf+t)'], errors='coerce')
df['K(tf+t)'] = pd.to_numeric(df['K(tf+t)'], errors='coerce')
df['expression_value'] = pd.to_numeric(df['expression_value'], errors='coerce')

# Удаляем строки с пропусками в базовых критичных столбцах
df = df.dropna(subset=['TF', 'tissue', 'N(tf+t)', 'K(tf+t)']).copy()

print(f'   Всего строк после удаления пропусков в базовых столбцах: {len(df)}')
print(f'   Уникальных TF до фильтрации по экспрессии: {df["TF"].nunique()}')

# Фильтрация по экспрессии
df = df[df['expression_value'] >= EXPRESSION_THRESHOLD].copy()

if len(df) == 0:
    raise ValueError(
        f'После фильтрации expression_value >= {EXPRESSION_THRESHOLD} не осталось строк.'
    )

print(f'   Строк после фильтрации expression_value >= {EXPRESSION_THRESHOLD}: {len(df)}')
print(f'   Уникальных TF после фильтрации по экспрессии: {df["TF"].nunique()}')
print(f'   Уникальных тканей после фильтрации по экспрессии: {df["tissue"].nunique()}')

# Пересчёт доли K/N
df['prop_K_N'] = np.where(
    df['N(tf+t)'] > 0,
    df['K(tf+t)'] / df['N(tf+t)'],
    np.nan
)

df = df.dropna(subset=['prop_K_N']).copy()

if len(df) == 0:
    raise ValueError('После расчёта prop_K_N не осталось строк с корректными значениями.')

# 2. Определение глобального порога топ-5% по prop_K_N
print('2. Расчёт глобального топ-5% среди строк с expression_value >= 1...')

threshold_95 = np.percentile(df['prop_K_N'], 95)

# В топ-5% попадают все строки с prop_K_N >= 95-го перцентиля
df['in_top5'] = df['prop_K_N'] >= threshold_95

n_top = int(df['in_top5'].sum())   # n: всего TF×tissue в топ-5%
N_total = int(len(df))             # N: всего TF×tissue в анализируемом наборе

print(f'   Порог 95-го перцентиля: {threshold_95:.6f}')
print(f'   n, строк в топ-5%: {n_top}')
print(f'   N, всего строк после фильтрации по экспрессии: {N_total}')
print(f'   Реальная доля строк в топе: {n_top / N_total:.4f}')

# 2.1. Сохранение всех пар TF×tissue, попавших в топ-5%
print('2.1. Сохранение всех пар TF×tissue, попавших в топ-5%...')

top5_percent_pairs = df[df['in_top5']].copy()
top5_percent_pairs = top5_percent_pairs.sort_values(
    'prop_K_N',
    ascending=False
).reset_index(drop=True)

save_cols = [
    'TF',
    'tissue',
    'N(tf+t)',
    'K(tf+t)',
    'prop_K_N',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'expression_value',
    'in_top5'
]

if 'expression_unit' in df.columns:
    save_cols.append('expression_unit')

if 'expression_source' in df.columns:
    save_cols.append('expression_source')

top5_percent_pairs[save_cols].to_csv(
    TOP5_PAIRS_FILE,
    sep='\t',
    index=False
)

print(f'   Сохранено {len(top5_percent_pairs)} пар в файл: {TOP5_PAIRS_FILE}')
print(top5_percent_pairs[['TF', 'tissue', 'prop_K_N', 'expression_value']].head(10).to_string(index=False))

# 3. Агрегация по каждому TF
print('3. Подсчёт k и K для каждого TF...')

tf_agg = df.groupby('TF')['in_top5'].agg(
    k='sum',
    total_obs='count'
).reset_index()

# k — сколько раз TF попал в топ-5%
# K — сколько раз TF не попал в топ-5%
tf_agg['K'] = tf_agg['total_obs'] - tf_agg['k']
tf_agg.drop(columns=['total_obs'], inplace=True)

tf_agg['k'] = tf_agg['k'].astype(int)
tf_agg['K'] = tf_agg['K'].astype(int)

# Подтягиваем метаданные TF
tf_meta = df[
    [
        'TF',
        'is_pioneer',
        'tfclass:class',
        'tfclass:family'
    ]
].drop_duplicates(subset=['TF']).copy()

# Проверка, не различаются ли метаданные у одного TF
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
    print('ВНИМАНИЕ: у некоторых TF различаются метаданные:')
    print(bad_meta.to_string(index=False))

tf_agg = tf_agg.merge(tf_meta, on='TF', how='left')

# Добавляем глобальные параметры гипергеометрического теста
tf_agg['n'] = n_top
tf_agg['N'] = N_total

# 4. Гипергеометрический тест
print('4. Расчёт p-value гипергеометрического теста...')

tf_agg['p_value'] = tf_agg.apply(
    lambda row: hypergeom.sf(
        int(row['k']) - 1,
        int(row['N']),
        int(row['n']),
        int(row['k'] + row['K'])
    ),
    axis=1
)

# 5. Частоты и fold enrichment
print('5. Расчёт observed_frequency, background_frequency и fold_enrichment...')

tf_agg['observed_frequency'] = tf_agg['k'] / (tf_agg['k'] + tf_agg['K'])

background_frequency = n_top / N_total if N_total > 0 else np.nan
tf_agg['background_frequency'] = background_frequency

tf_agg['fold_enrichment'] = np.where(
    tf_agg['background_frequency'] > 0,
    tf_agg['observed_frequency'] / tf_agg['background_frequency'],
    np.nan
)

# 6. FDR-коррекция Benjamini-Hochberg
print('6. FDR-коррекция Benjamini-Hochberg...')

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
print('7. Сортировка и сохранение результата...')

tf_agg = tf_agg.sort_values(
    ['p_value_fdr', 'p_value', 'fold_enrichment'],
    ascending=[True, True, False]
).reset_index(drop=True)

out_cols = [
    'TF',
    'is_pioneer',
    'tfclass:class',
    'tfclass:family',
    'k',
    'K',
    'n',
    'N',
    'observed_frequency',
    'background_frequency',
    'p_value',
    'p_value_fdr',
    'fold_enrichment'
]

tf_agg[out_cols].to_csv(
    OUTPUT_FILE,
    sep='\t',
    index=False
)

print('\nГотово!')
print(f'Основная таблица сохранена в: {OUTPUT_FILE}')
print(f'Пары TF×tissue из топ-5% сохранены в: {TOP5_PAIRS_FILE}')
print(f'TF в итоговой таблице: {len(tf_agg)}')
print(f'Значимых TF при FDR < 0.05: {(tf_agg["p_value_fdr"] < 0.05).sum()}')
print(f'Медианный fold_enrichment: {tf_agg["fold_enrichment"].median():.3f}')

print('\nТоп-10 TF по FDR:')
print(
    tf_agg[
        [
            'TF',
            'k',
            'K',
            'observed_frequency',
            'background_frequency',
            'fold_enrichment',
            'p_value',
            'p_value_fdr'
        ]
    ].head(10).to_string(index=False)
)