#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt


# =========================================================
# НАСТРОЙКИ
# =========================================================

INPUT_FILE = 'motifbreakr_analyzed_001_log.tsv'

# Таблица с колонками:
# TF    is_pioneer
PIONEER_ANNOTATION_FILE = 'chi2_motifbr_per_TF_with_pioneer.tsv'

OUTPUT_TF_FRACTIONS = 'motifbreakr_broken_site_fraction_by_tf.tsv'
OUTPUT_SUMMARY = 'permutation_test_broken_sites_pioneer_vs_nonpioneer_summary.tsv'
OUTPUT_PERMUTATIONS = 'permutation_test_broken_sites_pioneer_vs_nonpioneer_null_distribution.tsv'
OUTPUT_PLOT = 'pioneer_vs_nonpioneer_mean_median_with_errorbars.png'

N_PERMUTATIONS = 100000
RANDOM_SEED = 42

# Цвета
PIONEER_COLOR = "#9370DB"
NON_PIONEER_COLOR = "skyblue"

# Варианты:
# 'two-sided'
# 'greater'
# 'less'
ALTERNATIVE = 'two-sided'

TF_COLUMN_MOTIFBREAKR = 'geneSymbol'
TF_COLUMN_ANNOTATION = 'TF'

COMPARISON_COLUMNS = [
    'comparison_motif_score_0.001',
    'comparison_motif_score_0.0005',
    'comparison_motif_score_0.0001'
]

# motifbreakR:
#  2  = strong
# -2  = strong decrease
BROKEN_VALUES = {-2, 2}

# =========================================================


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Удаляет технический индексный столбец
    и очищает названия колонок.
    """

    first_col = str(df.columns[0]).strip()

    if first_col.startswith('Unnamed') or first_col == '':
        df = df.drop(columns=df.columns[0])

    df.columns = [str(c).strip() for c in df.columns]

    return df


def parse_bool_column(series: pd.Series) -> pd.Series:
    """
    Преобразует колонку is_pioneer в bool.
    """

    if series.dtype == bool:
        return series

    mapping = {
        'true': True,
        'false': False,
        '1': True,
        '0': False,
        'yes': True,
        'no': False,
        'y': True,
        'n': False,
    }

    parsed = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
        .map(mapping)
    )

    if parsed.isna().any():

        bad_values = sorted(
            series.loc[parsed.isna()]
            .astype(str)
            .unique()
        )

        raise ValueError(
            'Не удалось распознать значения '
            'в колонке is_pioneer:\n'
            + ', '.join(bad_values)
        )

    return parsed.astype(bool)


def calculate_p_value(
    observed_stat,
    permuted_stats,
    alternative
):
    """
    Эмпирический p-value.
    """

    if alternative == 'two-sided':

        extreme_count = np.sum(
            np.abs(permuted_stats)
            >= abs(observed_stat)
        )

    elif alternative == 'greater':

        extreme_count = np.sum(
            permuted_stats >= observed_stat
        )

    elif alternative == 'less':

        extreme_count = np.sum(
            permuted_stats <= observed_stat
        )

    else:

        raise ValueError(
            "ALTERNATIVE должен быть "
            "'two-sided', 'greater' или 'less'"
        )

    p_value = (
        (extreme_count + 1)
        / (len(permuted_stats) + 1)
    )

    return p_value


def plot_mean_median_with_errorbars(
    pioneer_values,
    nonpioneer_values,
    output_file
):
    """
    Boxplot с отображением:
    - распределения
    - среднего
    - медианы
    """

    plt.rcParams['font.family'] = 'Arial'

    fig, ax = plt.subplots(
        figsize=(7, 6),
        dpi=300
    )

    box = ax.boxplot(
        [pioneer_values, nonpioneer_values],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
        medianprops=dict(
            color='black',
            linewidth=2
        ),
        whiskerprops=dict(
            linewidth=1.5
        ),
        capprops=dict(
            linewidth=1.5
        ),
        boxprops=dict(
            linewidth=1.5
        )
    )

    colors = [
        PIONEER_COLOR,
        NON_PIONEER_COLOR
    ]

    for patch, color in zip(box['boxes'], colors):

        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    pioneer_mean = np.mean(pioneer_values)
    nonpioneer_mean = np.mean(nonpioneer_values)

    pioneer_median = np.median(pioneer_values)
    nonpioneer_median = np.median(nonpioneer_values)

    # Средние
    ax.scatter(
        [1, 2],
        [pioneer_mean, nonpioneer_mean],
        color='red',
        s=80,
        zorder=5,
        label='Среднее'
    )

    # Медианы
    ax.scatter(
        [1, 2],
        [pioneer_median, nonpioneer_median],
        color='black',
        s=80,
        zorder=5,
        label='Медиана'
    )

    # Подписи
    ax.text(
        1.08,
        pioneer_mean,
        f'среднее = {pioneer_mean:.3f}',
        fontsize=10,
        va='center'
    )

    ax.text(
        2.08,
        nonpioneer_mean,
        f'среднее = {nonpioneer_mean:.3f}',
        fontsize=10,
        va='center'
    )

    ax.text(
        1.08,
        pioneer_median,
        f'медиана = {pioneer_median:.3f}',
        fontsize=10,
        va='center'
    )

    ax.text(
        2.08,
        nonpioneer_median,
        f'медиана = {nonpioneer_median:.3f}',
        fontsize=10,
        va='center'
    )

    ax.set_xticks([1, 2])

    ax.set_xticklabels([
        'Пионерные ТФ',
        'Не-пионерные ТФ'
    ])

    ax.set_ylabel(
        'Доля сломанных сайтов'
    )

    ax.set_title(
        'Распределение доли сломанных сайтов'
    )

    ax.grid(
        axis='y',
        linestyle='--',
        alpha=0.4
    )

    ax.legend(
        frameon=False,
        fontsize=10,
        loc='upper right'
    )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    plt.savefig(
        output_file,
        bbox_inches='tight'
    )

    plt.close()

    print(f'   График сохранён: {output_file}')


def main():

    # =====================================================
    # ЗАГРУЗКА motifbreakR
    # =====================================================

    print('1. Загрузка motifbreakR-таблицы...')

    input_path = Path(INPUT_FILE)

    if not input_path.exists():

        raise FileNotFoundError(
            f'Файл не найден: {INPUT_FILE}'
        )

    df = pd.read_csv(
        input_path,
        sep='\t',
        engine='python'
    )

    df = clean_columns(df)

    required_motif_cols = (
        [TF_COLUMN_MOTIFBREAKR]
        + COMPARISON_COLUMNS
    )

    missing_motif_cols = [
        c
        for c in required_motif_cols
        if c not in df.columns
    ]

    if missing_motif_cols:

        raise ValueError(
            'Отсутствуют обязательные колонки:\n'
            + ', '.join(missing_motif_cols)
        )

    print(f'   Всего строк: {len(df)}')

    # =====================================================
    # ОЧИСТКА TF
    # =====================================================

    print('2. Очистка TF...')

    df[TF_COLUMN_MOTIFBREAKR] = (
        df[TF_COLUMN_MOTIFBREAKR]
        .astype(str)
        .str.strip()
    )

    df = df[
        df[TF_COLUMN_MOTIFBREAKR].notna()
    ].copy()

    df = df[
        df[TF_COLUMN_MOTIFBREAKR] != ''
    ].copy()

    df = df[
        df[TF_COLUMN_MOTIFBREAKR]
        .str.lower() != 'nan'
    ].copy()

    if len(df) == 0:

        raise ValueError(
            'После фильтрации не осталось строк.'
        )

    print(f'   Строк после очистки TF: {len(df)}')

    # =====================================================
    # ОПРЕДЕЛЕНИЕ СЛОМАННЫХ САЙТОВ
    # =====================================================

    print('3. Определение сломанных сайтов...')

    for col in COMPARISON_COLUMNS:

        df[col] = pd.to_numeric(
            df[col],
            errors='coerce'
        )

    comparison_matrix = df[COMPARISON_COLUMNS]

    df['is_broken_site'] = (
        comparison_matrix
        .isin(BROKEN_VALUES)
        .any(axis=1)
    )

    n_total_rows = len(df)

    n_broken_rows = int(
        df['is_broken_site'].sum()
    )

    broken_fraction = (
        n_broken_rows / n_total_rows
    )

    print(f'   Всего строк: {n_total_rows}')
    print(f'   Сломанных сайтов: {n_broken_rows}')
    print(f'   Доля: {broken_fraction:.6f}')

    # =====================================================
    # АГРЕГАЦИЯ ПО TF
    # =====================================================

    print('4. Расчёт доли сломанных сайтов по TF...')

    tf_fraction = (
        df
        .groupby(
            TF_COLUMN_MOTIFBREAKR,
            as_index=False
        )
        .agg(
            N_sites=(
                TF_COLUMN_MOTIFBREAKR,
                'size'
            ),
            K_broken_sites=(
                'is_broken_site',
                'sum'
            )
        )
        .rename(
            columns={
                TF_COLUMN_MOTIFBREAKR: 'TF'
            }
        )
    )

    tf_fraction['K_broken_sites'] = (
        tf_fraction['K_broken_sites']
        .astype(int)
    )

    tf_fraction['broken_site_fraction'] = (
        tf_fraction['K_broken_sites']
        / tf_fraction['N_sites']
    )

    print(f'   TF после агрегации: {len(tf_fraction)}')

    # =====================================================
    # ЗАГРУЗКА АННОТАЦИИ pioneer
    # =====================================================

    print('5. Загрузка pioneer-аннотации...')

    annotation_path = Path(
        PIONEER_ANNOTATION_FILE
    )

    if not annotation_path.exists():

        raise FileNotFoundError(
            f'Файл не найден:\n'
            f'{PIONEER_ANNOTATION_FILE}'
        )

    annotation = pd.read_csv(
        annotation_path,
        sep='\t',
        engine='python'
    )

    annotation = clean_columns(annotation)

    required_annotation_cols = [
        TF_COLUMN_ANNOTATION,
        'is_pioneer'
    ]

    missing_annotation_cols = [
        c
        for c in required_annotation_cols
        if c not in annotation.columns
    ]

    if missing_annotation_cols:

        raise ValueError(
            'В annotation отсутствуют колонки:\n'
            + ', '.join(missing_annotation_cols)
        )

    annotation[TF_COLUMN_ANNOTATION] = (
        annotation[TF_COLUMN_ANNOTATION]
        .astype(str)
        .str.strip()
    )

    annotation['is_pioneer'] = (
        parse_bool_column(
            annotation['is_pioneer']
        )
    )

    annotation_unique = (
        annotation[
            [TF_COLUMN_ANNOTATION, 'is_pioneer']
        ]
        .drop_duplicates(
            subset=[TF_COLUMN_ANNOTATION]
        )
        .rename(
            columns={
                TF_COLUMN_ANNOTATION: 'TF'
            }
        )
    )

    tf_fraction = pd.merge(
        tf_fraction,
        annotation_unique,
        on='TF',
        how='left'
    )

    tf_fraction_for_test = (
        tf_fraction
        .dropna(subset=['is_pioneer'])
        .copy()
    )

    tf_fraction_for_test['is_pioneer'] = (
        tf_fraction_for_test['is_pioneer']
        .astype(bool)
    )

    analysis_col = 'broken_site_fraction'

    pioneer_values = (
        tf_fraction_for_test
        .loc[
            tf_fraction_for_test['is_pioneer'],
            analysis_col
        ]
        .to_numpy()
    )

    nonpioneer_values = (
        tf_fraction_for_test
        .loc[
            ~tf_fraction_for_test['is_pioneer'],
            analysis_col
        ]
        .to_numpy()
    )

    if len(pioneer_values) == 0:

        raise ValueError(
            'Нет pioneer TF.'
        )

    if len(nonpioneer_values) == 0:

        raise ValueError(
            'Нет non-pioneer TF.'
        )

    print(f'   Pioneer TF: {len(pioneer_values)}')
    print(f'   Non-pioneer TF: {len(nonpioneer_values)}')

    # =====================================================
    # ПЕРМУТАЦИОННЫЙ ТЕСТ
    # =====================================================

    print('6. Пермутационный тест...')

    observed_stat = (
        np.median(pioneer_values)
        - np.median(nonpioneer_values)
    )

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    all_values = np.concatenate([
        pioneer_values,
        nonpioneer_values
    ])

    n_pioneer = len(pioneer_values)

    permuted_stats = np.empty(
        N_PERMUTATIONS
    )

    for i in range(N_PERMUTATIONS):

        shuffled = rng.permutation(
            all_values
        )

        perm_pioneer = shuffled[:n_pioneer]

        perm_nonpioneer = shuffled[n_pioneer:]

        permuted_stats[i] = (
            np.median(perm_pioneer)
            - np.median(perm_nonpioneer)
        )

    p_value = calculate_p_value(
        observed_stat,
        permuted_stats,
        ALTERNATIVE
    )

    print(
        f'   Observed statistic: '
        f'{observed_stat:.6f}'
    )

    print(
        f'   p-value: '
        f'{p_value:.6g}'
    )

    # =====================================================
    # СОХРАНЕНИЕ TF ТАБЛИЦЫ
    # =====================================================

    print('7. Сохранение таблиц...')

    tf_fraction.to_csv(
        OUTPUT_TF_FRACTIONS,
        sep='\t',
        index=False
    )

    pd.DataFrame({
        'permuted_statistic': permuted_stats
    }).to_csv(
        OUTPUT_PERMUTATIONS,
        sep='\t',
        index=False
    )

    summary_df = pd.DataFrame([{
        'observed_statistic': observed_stat,
        'p_value': p_value,
        'n_pioneer': len(pioneer_values),
        'n_nonpioneer': len(nonpioneer_values),
        'median_pioneer': np.median(pioneer_values),
        'median_nonpioneer': np.median(nonpioneer_values),
        'mean_pioneer': np.mean(pioneer_values),
        'mean_nonpioneer': np.mean(nonpioneer_values)
    }])

    summary_df.to_csv(
        OUTPUT_SUMMARY,
        sep='\t',
        index=False
    )

    print(f'   Сохранён: {OUTPUT_TF_FRACTIONS}')
    print(f'   Сохранён: {OUTPUT_PERMUTATIONS}')
    print(f'   Сохранён: {OUTPUT_SUMMARY}')

    # =====================================================
    # ГРАФИК
    # =====================================================

    print('8. Построение графика...')

    plot_mean_median_with_errorbars(
        pioneer_values,
        nonpioneer_values,
        OUTPUT_PLOT
    )

    print('\nГотово.')


if __name__ == '__main__':
    main()