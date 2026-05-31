#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

INPUT_FILE = 'total_tf_snv_pioneer_class_family.tsv'
OUTPUT_FILE = 'mannwhitney_by_tf_family.tsv'
SKIPPED_FILE = 'mannwhitney_by_tf_family_skipped.tsv'

MIN_GROUP_SIZE = 3


def parse_bool(value):
    if pd.isna(value):
        return np.nan

    value = str(value).strip().lower()

    if value in ['true', '1', 'yes', 'y', 't']:
        return True

    if value in ['false', '0', 'no', 'n', 'f']:
        return False

    return np.nan


def benjamini_hochberg(pvalues):
    pvalues = np.array(pvalues, dtype=float)
    qvalues = np.full(len(pvalues), np.nan)

    valid = ~np.isnan(pvalues)
    valid_pvalues = pvalues[valid]

    m = len(valid_pvalues)

    if m == 0:
        return qvalues

    order = np.argsort(valid_pvalues)
    sorted_pvalues = valid_pvalues[order]

    sorted_qvalues = sorted_pvalues * m / np.arange(1, m + 1)
    sorted_qvalues = np.minimum.accumulate(sorted_qvalues[::-1])[::-1]
    sorted_qvalues = np.clip(sorted_qvalues, 0, 1)

    corrected = np.empty(m)
    corrected[order] = sorted_qvalues

    qvalues[valid] = corrected

    return qvalues


print('1. Загрузка данных...')

df = pd.read_csv(INPUT_FILE, sep='\t', dtype=str)
df.columns = [str(col).strip() for col in df.columns]

required_cols = ['TF', 'is_pioneer', 'K/N', 'TF_family']
missing_cols = [col for col in required_cols if col not in df.columns]

if len(missing_cols) > 0:
    raise ValueError(f'В таблице отсутствуют обязательные столбцы: {missing_cols}')

print(f'   Всего строк в исходной таблице: {len(df)}')

print('2. Подготовка данных...')

df['is_pioneer'] = df['is_pioneer'].apply(parse_bool)

df['K/N'] = (
    df['K/N']
    .astype(str)
    .str.strip()
    .str.replace(',', '.', regex=False)
    .str.replace('%', '', regex=False)
)

df['K/N'] = pd.to_numeric(df['K/N'], errors='coerce')

df['TF_family'] = df['TF_family'].astype(str).str.strip()
df.loc[df['TF_family'].isin(['', 'nan', 'None']), 'TF_family'] = np.nan

before_filter = len(df)

df = df.dropna(subset=['is_pioneer', 'K/N', 'TF_family'])

after_filter = len(df)

print(f'   Строк до удаления NA/некорректных значений: {before_filter}')
print(f'   Строк после удаления NA/некорректных значений: {after_filter}')
print(f'   Уникальных семейств TF: {df["TF_family"].nunique()}')

results = []
skipped = []

print('3. Тест Манна–Уитни по семействам TF...')

for family, family_df in df.groupby('TF_family'):
    pioneer_values = family_df.loc[family_df['is_pioneer'] == True, 'K/N'].dropna()
    non_pioneer_values = family_df.loc[family_df['is_pioneer'] == False, 'K/N'].dropna()

    n_pioneer = len(pioneer_values)
    n_non_pioneer = len(non_pioneer_values)

    if n_pioneer < MIN_GROUP_SIZE or n_non_pioneer < MIN_GROUP_SIZE:
        skipped.append({
            'TF_family': family,
            'n_pioneer': n_pioneer,
            'n_non_pioneer': n_non_pioneer,
            'reason': f'размер одной из групп меньше {MIN_GROUP_SIZE}'
        })
        continue

    stat, pvalue = mannwhitneyu(
        pioneer_values,
        non_pioneer_values,
        alternative='two-sided',
        method='auto'
    )

    median_pioneer = pioneer_values.median()
    median_non_pioneer = non_pioneer_values.median()

    mean_pioneer = pioneer_values.mean()
    mean_non_pioneer = non_pioneer_values.mean()

    diff_median = median_pioneer - median_non_pioneer
    diff_mean = mean_pioneer - mean_non_pioneer

    rank_biserial = (2 * stat) / (n_pioneer * n_non_pioneer) - 1

    if diff_median > 0:
        direction = 'выше у пионерных TF'
    elif diff_median < 0:
        direction = 'выше у не-пионерных TF'
    else:
        direction = 'медианы равны'

    results.append({
        'TF_family': family,
        'n_pioneer': n_pioneer,
        'n_non_pioneer': n_non_pioneer,
        'median_pioneer_KN': median_pioneer,
        'median_non_pioneer_KN': median_non_pioneer,
        'diff_median_pioneer_minus_non_pioneer': diff_median,
        'mean_pioneer_KN': mean_pioneer,
        'mean_non_pioneer_KN': mean_non_pioneer,
        'diff_mean_pioneer_minus_non_pioneer': diff_mean,
        'mannwhitney_U': stat,
        'p_value': pvalue,
        'rank_biserial_correlation': rank_biserial,
        'direction_by_median': direction
    })

results_df = pd.DataFrame(results)

if len(results_df) > 0:
    results_df['p_adj_BH'] = benjamini_hochberg(results_df['p_value'])
    results_df = results_df.sort_values(
        by=['p_adj_BH', 'p_value'],
        ascending=[True, True]
    )

skipped_df = pd.DataFrame(skipped)

results_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
skipped_df.to_csv(SKIPPED_FILE, sep='\t', index=False)

print('4. Готово.')
print(f'   Минимальный размер каждой группы: {MIN_GROUP_SIZE}')
print(f'   Семейств протестировано: {len(results_df)}')
print(f'   Семейств пропущено: {len(skipped_df)}')
print(f'   Основной файл результатов: {OUTPUT_FILE}')
print(f'   Файл с пропущенными семействами: {SKIPPED_FILE}')

if len(results_df) > 0:
    print('')
    print('Топ-10 семейств по скорректированному p-value:')
    print(
        results_df[
            [
                'TF_family',
                'n_pioneer',
                'n_non_pioneer',
                'median_pioneer_KN',
                'median_non_pioneer_KN',
                'diff_median_pioneer_minus_non_pioneer',
                'p_value',
                'p_adj_BH',
                'direction_by_median'
            ]
        ].head(10).to_string(index=False)
    )