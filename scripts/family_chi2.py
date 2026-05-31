#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

INPUT_FILE = 'total_tf_snv_pioneer_class_family.tsv'

OUTPUT_FILE = 'chi_square_by_tf_family.tsv'
SKIPPED_FILE = 'chi_square_by_tf_family_skipped.tsv'

MIN_TF_PER_GROUP = 3


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

required_cols = ['TF', 'N(tf)', 'K(tf)', 'is_pioneer', 'TF_family']
missing_cols = [col for col in required_cols if col not in df.columns]

if len(missing_cols) > 0:
    raise ValueError(f'В таблице отсутствуют обязательные столбцы: {missing_cols}')

print(f'   Всего строк: {len(df)}')

print('2. Подготовка данных...')

df['is_pioneer'] = df['is_pioneer'].apply(parse_bool)

df['N(tf)'] = pd.to_numeric(
    df['N(tf)'].astype(str).str.strip().str.replace(',', '.', regex=False),
    errors='coerce'
)

df['K(tf)'] = pd.to_numeric(
    df['K(tf)'].astype(str).str.strip().str.replace(',', '.', regex=False),
    errors='coerce'
)

df['TF_family'] = df['TF_family'].astype(str).str.strip()
df.loc[df['TF_family'].isin(['', 'nan', 'None']), 'TF_family'] = np.nan

before_filter = len(df)

df = df.dropna(subset=['TF', 'N(tf)', 'K(tf)', 'is_pioneer', 'TF_family'])

df['N(tf)'] = df['N(tf)'].astype(int)
df['K(tf)'] = df['K(tf)'].astype(int)

df = df[df['N(tf)'] >= 0]
df = df[df['K(tf)'] >= 0]
df = df[df['K(tf)'] <= df['N(tf)']]

after_filter = len(df)

print(f'   Строк до фильтрации: {before_filter}')
print(f'   Строк после фильтрации: {after_filter}')
print(f'   Уникальных семейств TF: {df["TF_family"].nunique()}')

results = []
skipped = []

print('3. χ²-тест по семействам TF...')

for family, family_df in df.groupby('TF_family'):
    pioneer_df = family_df[family_df['is_pioneer'] == True]
    non_pioneer_df = family_df[family_df['is_pioneer'] == False]

    n_tf_pioneer = pioneer_df['TF'].nunique()
    n_tf_non_pioneer = non_pioneer_df['TF'].nunique()

    if n_tf_pioneer < MIN_TF_PER_GROUP or n_tf_non_pioneer < MIN_TF_PER_GROUP:
        skipped.append({
            'TF_family': family,
            'n_tf_pioneer': n_tf_pioneer,
            'n_tf_non_pioneer': n_tf_non_pioneer,
            'reason': f'меньше {MIN_TF_PER_GROUP} TF в одной из групп'
        })
        continue

    K_pioneer = int(pioneer_df['K(tf)'].sum())
    N_pioneer = int(pioneer_df['N(tf)'].sum())
    no_snv_pioneer = N_pioneer - K_pioneer

    K_non_pioneer = int(non_pioneer_df['K(tf)'].sum())
    N_non_pioneer = int(non_pioneer_df['N(tf)'].sum())
    no_snv_non_pioneer = N_non_pioneer - K_non_pioneer

    if N_pioneer == 0 or N_non_pioneer == 0:
        skipped.append({
            'TF_family': family,
            'n_tf_pioneer': n_tf_pioneer,
            'n_tf_non_pioneer': n_tf_non_pioneer,
            'reason': 'нулевое суммарное N(tf) в одной из групп'
        })
        continue

    table = np.array([
        [K_pioneer, no_snv_pioneer],
        [K_non_pioneer, no_snv_non_pioneer]
    ])

    if np.any(table < 0):
        skipped.append({
            'TF_family': family,
            'n_tf_pioneer': n_tf_pioneer,
            'n_tf_non_pioneer': n_tf_non_pioneer,
            'reason': 'получилось отрицательное значение N-K'
        })
        continue

    if np.any(table.sum(axis=0) == 0) or np.any(table.sum(axis=1) == 0):
        skipped.append({
            'TF_family': family,
            'n_tf_pioneer': n_tf_pioneer,
            'n_tf_non_pioneer': n_tf_non_pioneer,
            'reason': 'нулевая строка или колонка в таблице сопряженности'
        })
        continue

    chi2, p_value, dof, expected = chi2_contingency(
        table,
        correction=False
    )

    fisher_or, fisher_p_value = fisher_exact(
        table,
        alternative='two-sided'
    )

    prop_pioneer = K_pioneer / N_pioneer
    prop_non_pioneer = K_non_pioneer / N_non_pioneer

    diff_prop = prop_pioneer - prop_non_pioneer

    odds_ratio = (
        ((K_pioneer + 0.5) / (no_snv_pioneer + 0.5)) /
        ((K_non_pioneer + 0.5) / (no_snv_non_pioneer + 0.5))
    )

    total_count = table.sum()
    phi = np.sqrt(chi2 / total_count)

    min_expected = expected.min()

    if diff_prop > 0:
        direction = 'выше у пионерных TF'
    elif diff_prop < 0:
        direction = 'выше у не-пионерных TF'
    else:
        direction = 'доли равны'

    results.append({
        'TF_family': family,

        'n_tf_pioneer': n_tf_pioneer,
        'n_tf_non_pioneer': n_tf_non_pioneer,

        'K_pioneer': K_pioneer,
        'N_pioneer': N_pioneer,
        'no_snv_pioneer': no_snv_pioneer,
        'proportion_pioneer': prop_pioneer,

        'K_non_pioneer': K_non_pioneer,
        'N_non_pioneer': N_non_pioneer,
        'no_snv_non_pioneer': no_snv_non_pioneer,
        'proportion_non_pioneer': prop_non_pioneer,

        'diff_proportion_pioneer_minus_non_pioneer': diff_prop,

        'chi2': chi2,
        'dof': dof,
        'p_value_chi2': p_value,

        'odds_ratio': odds_ratio,
        'fisher_exact_p_value': fisher_p_value,

        'phi_effect_size': phi,
        'min_expected_count': min_expected,
        'expected_K_pioneer': expected[0, 0],
        'expected_no_snv_pioneer': expected[0, 1],
        'expected_K_non_pioneer': expected[1, 0],
        'expected_no_snv_non_pioneer': expected[1, 1],

        'direction': direction
    })

results_df = pd.DataFrame(results)

if len(results_df) > 0:
    results_df['p_adj_BH_chi2'] = benjamini_hochberg(results_df['p_value_chi2'])

    results_df = results_df.sort_values(
        by=['p_adj_BH_chi2', 'p_value_chi2'],
        ascending=[True, True]
    )

skipped_df = pd.DataFrame(skipped)

results_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
skipped_df.to_csv(SKIPPED_FILE, sep='\t', index=False)

print('4. Готово.')
print(f'   Минимальное число TF в каждой группе: {MIN_TF_PER_GROUP}')
print(f'   Семейств протестировано: {len(results_df)}')
print(f'   Семейств пропущено: {len(skipped_df)}')
print(f'   Основной файл результатов: {OUTPUT_FILE}')
print(f'   Файл с пропущенными семействами: {SKIPPED_FILE}')

if len(results_df) > 0:
    print('')
    print('Топ семейств по p_adj_BH_chi2:')
    print(
        results_df[
            [
                'TF_family',
                'n_tf_pioneer',
                'n_tf_non_pioneer',
                'proportion_pioneer',
                'proportion_non_pioneer',
                'diff_proportion_pioneer_minus_non_pioneer',
                'odds_ratio',
                'p_value_chi2',
                'p_adj_BH_chi2',
                'fisher_exact_p_value',
                'min_expected_count',
                'direction'
            ]
        ].head(20).to_string(index=False)
    )