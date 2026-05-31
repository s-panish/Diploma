#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================

INPUT_FILE = 'motifbreakr_adastra_with_cell_lines_changed_expression_ge1.tsv'

OUTPUT_FILE = 'table_tf_broken_site_enrichment_adastra_expression_ge1.tsv'
BROKEN_SITES_FILE = 'broken_sites_adastra_expression_ge1.tsv'
TF_CELL_LINE_COUNTS_FILE = 'tf_cell_line_broken_counts_adastra_expression_ge1.tsv'

EXPRESSION_THRESHOLD = 1

COMPARISON_COLUMNS = [
    'comparison_motif_score_0.001',
    'comparison_motif_score_0.0005',
    'comparison_motif_score_0.0001'
]

BROKEN_VALUES = {-1, 1}

# =============================================


print('1. Загрузка данных...')

df = pd.read_csv(INPUT_FILE, sep='\t')
df.columns = [col.strip() for col in df.columns]

required_cols = [
    'geneSymbol',
    'cell_line',
    'expression_value'
] + COMPARISON_COLUMNS

missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    raise ValueError(
        f'В файле отсутствуют обязательные столбцы: {missing_cols}\n'
        f'Найденные столбцы: {list(df.columns)}'
    )

print(f'   Всего строк в исходном файле: {len(df)}')
print(f'   Уникальных TF в исходном файле: {df["geneSymbol"].nunique()}')
print(f'   Уникальных клеточных/тканевых контекстов: {df["cell_line"].nunique()}')

# Приведение expression_value к числовому формату
df['expression_value'] = pd.to_numeric(df['expression_value'], errors='coerce')

# Приведение comparison-столбцов к числовому формату
for col in COMPARISON_COLUMNS:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Удаляем строки с пропусками в критичных столбцах
df = df.dropna(subset=['geneSymbol', 'cell_line', 'expression_value']).copy()

print(f'   Строк после удаления пропусков в geneSymbol/cell_line/expression_value: {len(df)}')

# Фильтрация по экспрессии
df = df[df['expression_value'] >= EXPRESSION_THRESHOLD].copy()

if len(df) == 0:
    raise ValueError(
        f'После фильтрации expression_value >= {EXPRESSION_THRESHOLD} не осталось строк.'
    )

print(f'   Строк после фильтрации expression_value >= {EXPRESSION_THRESHOLD}: {len(df)}')
print(f'   Уникальных TF после фильтрации по экспрессии: {df["geneSymbol"].nunique()}')
print(f'   Уникальных клеточных/тканевых контекстов после фильтрации: {df["cell_line"].nunique()}')

# 2. Определение сломанных сайтов
print('2. Определение сломанных сайтов...')

broken_mask = pd.Series(False, index=df.index)

for col in COMPARISON_COLUMNS:
    broken_mask = broken_mask | df[col].isin(BROKEN_VALUES)

df['is_broken_site'] = broken_mask

n_broken = int(df['is_broken_site'].sum())   # n: всего сломанных сайтов
N_total = int(len(df))                       # N: всего сайтов в анализируемом наборе

if n_broken == 0:
    raise ValueError(
        'В файле не найдено ни одного сломанного сайта со значениями -1 или 1 '
        'в comparison_motif_score_* столбцах.'
    )

print(f'   Всего сайтов после фильтрации по экспрессии: {N_total}')
print(f'   Сломанных сайтов: {n_broken}')
print(f'   Доля сломанных сайтов: {n_broken / N_total:.6f}')

if n_broken == N_total:
    print(
        '   ВНИМАНИЕ: все строки в файле являются сломанными сайтами. '
        'Гипергеометрический тест в таком случае статистически неинформативен, '
        'поскольку отсутствует фон из несломанных сайтов.'
    )

# 2.1. Сохранение всех сломанных сайтов
print('2.1. Сохранение таблицы сломанных сайтов...')

broken_sites = df[df['is_broken_site']].copy()

sort_cols = ['geneSymbol', 'cell_line']
broken_sites = broken_sites.sort_values(sort_cols).reset_index(drop=True)

broken_sites.to_csv(
    BROKEN_SITES_FILE,
    sep='\t',
    index=False
)

print(f'   Сохранено сломанных сайтов: {len(broken_sites)}')
print(f'   Файл со сломанными сайтами: {BROKEN_SITES_FILE}')

# 3. Агрегация по каждому TF
print('3. Подсчёт k и K для каждого TF...')

tf_agg = df.groupby('geneSymbol')['is_broken_site'].agg(
    k='sum',
    total_sites='count'
).reset_index()

tf_agg = tf_agg.rename(columns={'geneSymbol': 'TF'})

# k — сколько сломанных сайтов у данного TF
# K — сколько несломанных сайтов у данного TF
tf_agg['K'] = tf_agg['total_sites'] - tf_agg['k']

tf_agg['k'] = tf_agg['k'].astype(int)
tf_agg['K'] = tf_agg['K'].astype(int)
tf_agg['total_sites'] = tf_agg['total_sites'].astype(int)

# Добавляем глобальные параметры гипергеометрического теста
tf_agg['n'] = n_broken
tf_agg['N'] = N_total

# 4. Гипергеометрический тест
print('4. Расчёт p-value гипергеометрического теста...')

# Для каждого TF:
# N — общее число сайтов в анализе
# n — общее число сломанных сайтов
# k + K — число всех сайтов данного TF
# k — число сломанных сайтов данного TF
#
# p-value показывает вероятность получить не меньше k сломанных сайтов
# у данного TF при случайном распределении сломанных сайтов по всем строкам.

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

background_frequency = n_broken / N_total if N_total > 0 else np.nan
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
    qvals[-1] = min(sorted_pvals[-1] * m / m, 1.0)

    for i in range(m - 2, -1, -1):
        qvals[i] = min(sorted_pvals[i] * m / (i + 1), qvals[i + 1])

    qvals = np.minimum(qvals, 1.0)

    result = np.empty(m, dtype=float)
    result[order] = qvals

    return result

tf_agg['p_value_fdr'] = fdr_bh(tf_agg['p_value'].values)

# 7. Дополнительная таблица по TF × cell_line
print('7. Подсчёт сломанных сайтов по TF × cell_line...')

tf_cell_line_counts = df.groupby(['geneSymbol', 'cell_line']).agg(
    total_sites=('is_broken_site', 'count'),
    broken_sites=('is_broken_site', 'sum'),
    mean_expression=('expression_value', 'mean')
).reset_index()

tf_cell_line_counts = tf_cell_line_counts.rename(
    columns={
        'geneSymbol': 'TF'
    }
)

tf_cell_line_counts['intact_sites'] = (
    tf_cell_line_counts['total_sites'] - tf_cell_line_counts['broken_sites']
)

tf_cell_line_counts['broken_fraction'] = np.where(
    tf_cell_line_counts['total_sites'] > 0,
    tf_cell_line_counts['broken_sites'] / tf_cell_line_counts['total_sites'],
    np.nan
)

tf_cell_line_counts = tf_cell_line_counts.sort_values(
    ['broken_fraction', 'broken_sites', 'total_sites'],
    ascending=[False, False, False]
).reset_index(drop=True)

tf_cell_line_counts.to_csv(
    TF_CELL_LINE_COUNTS_FILE,
    sep='\t',
    index=False
)

print(f'   Таблица TF × cell_line сохранена в: {TF_CELL_LINE_COUNTS_FILE}')

# 8. Сортировка и сохранение основной таблицы
print('8. Сортировка и сохранение результата...')

tf_agg = tf_agg.sort_values(
    ['p_value_fdr', 'p_value', 'fold_enrichment', 'k'],
    ascending=[True, True, False, False]
).reset_index(drop=True)

out_cols = [
    'TF',
    'k',
    'K',
    'total_sites',
    'n',
    'N',
    'observed_frequency',
    'background_frequency',
    'fold_enrichment',
    'p_value',
    'p_value_fdr'
]

tf_agg[out_cols].to_csv(
    OUTPUT_FILE,
    sep='\t',
    index=False
)

print('\nГотово!')
print(f'Основная таблица сохранена в: {OUTPUT_FILE}')
print(f'Сломанные сайты сохранены в: {BROKEN_SITES_FILE}')
print(f'Таблица TF × cell_line сохранена в: {TF_CELL_LINE_COUNTS_FILE}')
print(f'TF в итоговой таблице: {len(tf_agg)}')
print(f'Значимых TF при FDR < 0.05: {(tf_agg["p_value_fdr"] < 0.05).sum()}')
print(f'Фоновая доля сломанных сайтов: {background_frequency:.6f}')
print(f'Медианный fold_enrichment: {tf_agg["fold_enrichment"].median():.3f}')

print('\nТоп-10 TF по FDR:')
print(
    tf_agg[
        [
            'TF',
            'k',
            'K',
            'total_sites',
            'observed_frequency',
            'background_frequency',
            'fold_enrichment',
            'p_value',
            'p_value_fdr'
        ]
    ].head(10).to_string(index=False)
)