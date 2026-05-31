#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


# ===================== НАСТРОЙКИ =====================

INPUT_FILE = Path('motifbreakr_adastra_with_cell_lines.tsv')
TF_LIST_FILE = Path('tf_list.tsv')
CONTEXT_ALIASES_FILE = Path('context_aliases.tsv')

HPA_FILES = [
    {
        'source': 'cell_line',
        'path': Path('HPA/rna_celline.tsv'),
        'context_cols': ['Cell line'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'immune_cell',
        'path': Path('HPA/rna_immune_cell.tsv'),
        'context_cols': ['Immune cell'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'brain_pfc_subregion',
        'path': Path('HPA/rna_pfc_brain_hpa.tsv'),
        'context_cols': ['Subregion'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'brain_subregion',
        'path': Path('HPA/rna_brain_hpa.tsv'),
        'context_cols': ['Subregion'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'brain_region',
        'path': Path('HPA/rna_brain_region_hpa.tsv'),
        'context_cols': ['Brain region'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'tissue_consensus',
        'path': Path('HPA/rna_tissue_consensus.tsv'),
        'context_cols': ['Tissue'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'cell_line_cancer',
        'path': Path('HPA/rna_cell_line_cancer.tsv'),
        'context_cols': ['Cancer'],
        'value_col': 'nTPM',
        'unit': 'nTPM'
    },
    {
        'source': 'single_cell_type',
        'path': Path('HPA/rna_single_cell_type.tsv'),
        'context_cols': ['Cell type'],
        'value_col': 'nCPM',
        'unit': 'nCPM'
    },
    {
        'source': 'single_cell_type_group',
        'path': Path('HPA/rna_single_cell_type_group.tsv'),
        'context_cols': ['Cell type group', 'Cell type'],
        'value_col': 'nCPM',
        'unit': 'nCPM'
    }
]

COMPARISON_COLUMNS = [
    'comparison_motif_score_0.001',
    'comparison_motif_score_0.0005',
    'comparison_motif_score_0.0001'
]

CHANGED_VALUES = {-1, 1}
EXPRESSION_THRESHOLD = 1.0

USE_ONLY_CHANGED_SITES = False

OUTPUT_FILE = Path('motifbreakR_adastra_expression/motifbreakr_adastra_with_cell_lines_changed_expression_ge1.tsv')
UNMATCHED_FILE = Path('motifbreakR_adastra_expression/motifbreakr_adastra_with_cell_lines_changed_unmatched_expression.tsv')
LOW_EXPRESSION_FILE = Path('motifbreakR_adastra_expression/motifbreakr_adastra_with_cell_lines_changed_expression_lt1.tsv')
DIAGNOSTICS_FILE = Path('motifbreakR_adastra_expression/motifbreakr_adastra_with_cell_lines_changed_expression_diagnostics.tsv')

# ======================================================


DEFAULT_CONTEXT_ALIASES = [
    ('A3_Jurkat_clone_A3__Childhood_T_acute_lymphoblastic_leukemia', 'cell_line', 'JURKAT'),
    ('HeLa_S3__cervical_adenocarcinoma', 'cell_line', 'HeLa'),
    ('TALL1__T-Acute_Lymphoblastic_Leukemia', 'cell_line', 'TALL-1 [Human adult T-ALL]'),
    ('heart', 'tissue_consensus', 'heart muscle'),
    ('periphery_of_retina', 'tissue_consensus', 'retina'),
    ('pancreatic_islets', 'single_cell_type', 'pancreatic islet cells'),
    ('pancreatic_islets', 'tissue_consensus', 'pancreas'),
    ('MDM_monocyte_derived_macrophages', 'single_cell_type', 'macrophages'),
    ('MDM_monocyte_derived_macrophages_', 'single_cell_type', 'macrophages'),
    ('WA01__H1__human_embryonic_stem_cells', 'cell_line', 'H1'),
    ('WA01__H1__human_embryonic_stem_cells_', 'cell_line', 'H1'),
    ('WA01__H1__human_embryonic_stem_cells', 'single_cell_type', 'embryonic stem cells'),
    ('WA01__H1__human_embryonic_stem_cells_', 'single_cell_type', 'embryonic stem cells'),
    ('embryonic_stem_cells-derived_progenitors', 'single_cell_type', 'embryonic stem cells'),
    ('CD14+_monocytes', 'single_cell_type', 'monocytes'),
    ('CD14+_monocytes', 'immune_cell', 'monocytes'),
    ('CD19+_B_cells', 'single_cell_type', 'B-cells'),
    ('CD19+_B_cells', 'immune_cell', 'B-cells'),
    ('K562__myelogenous_leukemia', 'cell_line', 'K-562')
]

DEFAULT_TF_ALIASES = {
    'AP2D': 'TFAP2D',
    'BMAL1': 'BMAL1',
    'ZF64A': 'ZFP64',
    'ZN740': 'ZNF740',
    'ZN219': 'ZNF219'
}


def clean_value(x):
    if pd.isna(x):
        return ''
    return str(x).strip()


def strip_tail_underscores(x):
    return re.sub(r'_+$', '', clean_value(x))


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f'Не найден файл: {path}')

    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = [clean_value(c) for c in df.columns]
    return df


def to_float(x):
    x = clean_value(x).replace(',', '.')

    if x == '':
        return pd.NA

    try:
        return float(x)
    except ValueError:
        return pd.NA


def find_first_column(df, possible_cols, table_name):
    for col in possible_cols:
        if col in df.columns:
            return col

    raise ValueError(
        f'В таблице {table_name} не найден ни один из столбцов {possible_cols}. '
        f'Доступные столбцы: {list(df.columns)}'
    )


def normalize_gene(x):
    return clean_value(x).upper()


def normalize_general_context(x):
    x = strip_tail_underscores(x)
    x = x.replace('_', ' ')
    x = x.replace('-', ' ')
    x = x.replace('/', ' ')
    x = x.replace('\\', ' ')
    x = re.sub(r'[^A-Za-zА-Яа-я0-9]+', ' ', x)
    x = x.lower()
    x = re.sub(r'\s+', ' ', x).strip()
    return x


def split_context(x):
    x = strip_tail_underscores(x)
    parts = [p for p in re.split(r'__+', x) if clean_value(p) != '']

    if len(parts) == 0:
        return [x]

    return parts


def normalize_cell_line_context(x):
    x = strip_tail_underscores(x)
    parts = split_context(x)

    if len(parts) >= 2:
        x = parts[0]

    x = x.lower()
    x = re.sub(r'[^a-z0-9]+', '', x)
    return x


def add_unique(values, value):
    value = clean_value(value)

    if value != '' and value not in values:
        values.append(value)


def context_keys(context, source):
    keys = []

    if source == 'cell_line':
        add_unique(keys, normalize_cell_line_context(context))

        for part in split_context(context):
            part_key = part.lower()
            part_key = re.sub(r'[^a-z0-9]+', '', part_key)
            add_unique(keys, part_key)

        return keys

    base = normalize_general_context(context)
    add_unique(keys, base)

    if ' cells' in base:
        add_unique(keys, base.replace(' cells', ' cell'))

    if ' cell' in base and ' cells' not in base:
        add_unique(keys, base.replace(' cell', ' cells'))

    if ' lymphocytes' in base:
        add_unique(keys, base.replace(' lymphocytes', ' lymphocyte'))

    if ' lymphocyte' in base and ' lymphocytes' not in base:
        add_unique(keys, base.replace(' lymphocyte', ' lymphocytes'))

    manual_general = {
        'pancreatic islets': ['pancreatic islet cells', 'pancreatic islet', 'pancreas'],
        'pancreatic islet': ['pancreatic islet cells', 'pancreatic islets', 'pancreas'],
        'mdm monocyte derived macrophages': ['macrophages', 'monocyte derived macrophages'],
        'monocyte derived macrophages': ['macrophages'],
        'human embryonic stem cells': ['embryonic stem cells'],
        'h1 human embryonic stem cells': ['embryonic stem cells'],
        'wa01 h1 human embryonic stem cells': ['embryonic stem cells'],
        'embryonic stem cells derived progenitors': ['embryonic stem cells'],
        'hesc derived neurons': ['neurons'],
        'b cells': ['b cell', 'b-cells', 'b-cell'],
        'cd19 b cells': ['b cells', 'b cell', 'b-cells', 'b-cell'],
        't cells': ['t cell', 't-cells', 't-cell'],
        'cd4 t cells': ['cd4 t cell', 'cd4 t-cells', 'cd4 t-cell'],
        'cd8 t cells': ['cd8 t cell', 'cd8 t-cells', 'cd8 t-cell'],
        'monocytes': ['monocyte', 'classical monocyte', 'classical monocytes'],
        'cd14 monocytes': ['monocytes', 'monocyte', 'classical monocytes'],
        'natural killer cell': ['nk cell', 'nk-cell', 'nk cells', 'nk-cells'],
        'mature natural killer': ['natural killer cell', 'nk cell', 'nk cells'],
        'regulatory t cell': ['regulatory t cells', 't reg', 't-reg'],
        'fibroblast': ['fibroblasts'],
        'skin fibroblast': ['fibroblasts'],
        'foreskin fibroblast': ['fibroblasts'],
        'lung fibroblasts': ['fibroblasts'],
        'hepatocyte': ['hepatocytes'],
        'myocyte': ['myocytes'],
        'myoblast': ['myoblasts'],
        'trophoblast cell': ['trophoblast cells']
    }

    if base in manual_general:
        for value in manual_general[base]:
            add_unique(keys, normalize_general_context(value))
            add_unique(keys, value)

    return keys


def make_changed_mask(df):
    missing = [c for c in COMPARISON_COLUMNS if c not in df.columns]

    if missing:
        raise ValueError(
            f'В основной таблице отсутствуют comparison-столбцы: {missing}. '
            f'Доступные столбцы: {list(df.columns)}'
        )

    mask = pd.Series(False, index=df.index)

    for col in COMPARISON_COLUMNS:
        values = pd.to_numeric(
            df[col].astype(str).str.strip().str.replace(',', '.', regex=False),
            errors='coerce'
        )
        mask = mask | values.isin(CHANGED_VALUES)

    return mask


def load_tf_map():
    tf_map = {}

    if TF_LIST_FILE.exists():
        tf_df = read_tsv(TF_LIST_FILE)

        required = ['curated:uniprot_id', 'auto:gene_symbol']
        missing = [c for c in required if c not in tf_df.columns]

        if missing:
            raise ValueError(
                f'В {TF_LIST_FILE} отсутствуют столбцы: {missing}. '
                f'Доступные столбцы: {list(tf_df.columns)}'
            )

        tf_df = tf_df[required].copy()
        tf_df['curated:uniprot_id'] = tf_df['curated:uniprot_id'].apply(clean_value)
        tf_df['auto:gene_symbol'] = tf_df['auto:gene_symbol'].apply(clean_value)
        tf_df = tf_df[
            (tf_df['curated:uniprot_id'] != '') &
            (tf_df['auto:gene_symbol'] != '')
        ].copy()
        tf_df = tf_df.drop_duplicates(subset=['curated:uniprot_id'], keep='first')

        tf_map.update(dict(zip(tf_df['curated:uniprot_id'], tf_df['auto:gene_symbol'])))
        print(f'   Загружено TF-соответствий из {TF_LIST_FILE}: {len(tf_map)}')
    else:
        print(f'   {TF_LIST_FILE} не найден: geneSymbol будет использоваться напрямую.')

    tf_map.update(DEFAULT_TF_ALIASES)
    print(f'   Всего TF-соответствий с ручными заменами: {len(tf_map)}')

    return tf_map


def map_tf(tf, tf_map):
    tf = clean_value(tf)

    if tf == '':
        return ''

    if tf in tf_map:
        return tf_map[tf]

    upper_tf = tf.upper()

    for key, value in tf_map.items():
        if clean_value(key).upper() == upper_tf:
            return value

    return tf


def load_context_aliases():
    aliases = {}

    def add_alias(input_context, source, hpa_context):
        input_context = strip_tail_underscores(input_context)
        source = clean_value(source)
        hpa_context = clean_value(hpa_context)

        if input_context == '' or source == '' or hpa_context == '':
            return

        for key in [input_context, normalize_general_context(input_context)]:
            if key not in aliases:
                aliases[key] = []
            aliases[key].append((source, hpa_context))

    for input_context, source, hpa_context in DEFAULT_CONTEXT_ALIASES:
        add_alias(input_context, source, hpa_context)

    if CONTEXT_ALIASES_FILE.exists():
        alias_df = read_tsv(CONTEXT_ALIASES_FILE)
        required = ['input_context', 'hpa_source', 'hpa_context']
        missing = [c for c in required if c not in alias_df.columns]

        if missing:
            raise ValueError(
                f'В {CONTEXT_ALIASES_FILE} отсутствуют столбцы: {missing}. '
                f'Доступные столбцы: {list(alias_df.columns)}'
            )

        alias_df = alias_df[required].copy()

        for _, row in alias_df.iterrows():
            add_alias(row['input_context'], row['hpa_source'], row['hpa_context'])

        print(f'   Дополнительные соответствия контекстов из файла: {len(alias_df)}')
    else:
        print(f'   {CONTEXT_ALIASES_FILE} не найден: используются только соответствия из скрипта.')

    print(f'   Ключей в словаре контекстов: {len(aliases)}')
    return aliases


def load_hpa_table(info):
    df = read_tsv(info['path'])
    context_col = find_first_column(df, info['context_cols'], info['path'].name)

    required = ['Gene name', context_col, info['value_col']]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(
            f'В {info["path"]} отсутствуют столбцы: {missing}. '
            f'Доступные столбцы: {list(df.columns)}'
        )

    df = df[required].copy()
    df = df.rename(columns={context_col: 'HPA_context', info['value_col']: 'expression_value'})

    df['Gene name'] = df['Gene name'].apply(clean_value)
    df['HPA_context'] = df['HPA_context'].apply(clean_value)
    df['expression_value'] = df['expression_value'].apply(clean_value)
    df['expression_value_numeric'] = df['expression_value'].apply(to_float)

    df = df[
        (df['Gene name'] != '') &
        (df['HPA_context'] != '') &
        (df['expression_value'] != '') &
        (df['expression_value_numeric'].notna())
    ].copy()

    df['gene_key'] = df['Gene name'].apply(normalize_gene)
    df['source'] = info['source']
    df['expression_unit'] = info['unit']

    return df


def build_expression_lookup():
    lookup = {}

    for info in HPA_FILES:
        expr_df = load_hpa_table(info)
        print(f'   {info["path"]}: {len(expr_df)} строк')

        for _, row in expr_df.iterrows():
            for context_key in context_keys(row['HPA_context'], row['source']):
                key = (row['source'], row['gene_key'], context_key)
                value = float(row['expression_value_numeric'])

                record = {
                    'expression_value': row['expression_value'],
                    'expression_value_numeric': value,
                    'expression_unit': row['expression_unit'],
                    'expression_source': row['source'],
                    'matched_hpa_context': row['HPA_context']
                }

                if key not in lookup:
                    lookup[key] = record
                else:
                    if value > lookup[key]['expression_value_numeric']:
                        lookup[key] = record

    return lookup


def add_attempt(attempts, seen, source, context):
    for context_key in context_keys(context, source):
        key = (source, context_key)

        if key not in seen:
            seen.add(key)
            attempts.append((source, context_key))


def build_attempts(raw_context, aliases):
    attempts = []
    seen = set()

    raw_context = strip_tail_underscores(raw_context)
    norm_context = normalize_general_context(raw_context)

    for key in [raw_context, norm_context]:
        if key in aliases:
            for source, hpa_context in aliases[key]:
                add_attempt(attempts, seen, source, hpa_context)

    parts = split_context(raw_context)

    if '__' in raw_context:
        left_part = parts[0]
        right_part = parts[-1]

        add_attempt(attempts, seen, 'cell_line', left_part)
        add_attempt(attempts, seen, 'cell_line', raw_context)

        for part in parts:
            add_attempt(attempts, seen, 'cell_line', part)

        for source in [
            'cell_line_cancer',
            'tissue_consensus',
            'brain_pfc_subregion',
            'brain_subregion',
            'brain_region'
        ]:
            add_attempt(attempts, seen, source, right_part)

        add_attempt(attempts, seen, 'immune_cell', left_part)

        for part in parts:
            add_attempt(attempts, seen, 'single_cell_type', part)
            add_attempt(attempts, seen, 'single_cell_type_group', part)

        add_attempt(attempts, seen, 'single_cell_type', right_part)
        add_attempt(attempts, seen, 'single_cell_type_group', right_part)

    else:
        for info in HPA_FILES:
            add_attempt(attempts, seen, info['source'], raw_context)

    return attempts


def main():
    print('1. Загрузка основной таблицы')
    df = read_tsv(INPUT_FILE)

    gene_col = find_first_column(df, ['geneSymbol', 'TF', 'Gene name'], INPUT_FILE.name)
    context_col = find_first_column(df, ['cell_line', 'tissue', 'Tissue'], INPUT_FILE.name)

    print(f'   Файл: {INPUT_FILE}')
    print(f'   Строк всего: {len(df)}')
    print(f'   Столбец TF: {gene_col}')
    print(f'   Столбец контекста: {context_col}')

    if USE_ONLY_CHANGED_SITES:
        print('\n2. Фильтрация изменённых сайтов по comparison-столбцам')
        changed_mask = make_changed_mask(df)
        work_df = df.loc[changed_mask].copy()
        print(f'   Значения для отбора: {sorted(CHANGED_VALUES)}')
        print(f'   Строк после отбора изменённых сайтов: {len(work_df)}')
    else:
        print('\n2. Использование всех строк без фильтрации по comparison-столбцам')
        work_df = df.copy()
        print(f'   Строк взято в анализ: {len(work_df)}')

    print('\n3. Загрузка TF-словаря')
    tf_map = load_tf_map()

    work_df['_hpa_gene_symbol'] = work_df[gene_col].apply(lambda x: map_tf(x, tf_map))
    work_df['_gene_key'] = work_df['_hpa_gene_symbol'].apply(normalize_gene)

    print('\n4. Загрузка ручных соответствий контекстов')
    aliases = load_context_aliases()

    print('\n5. Загрузка HPA и создание индекса экспрессии')
    expression_lookup = build_expression_lookup()
    print(f'   Ключей в индексе HPA: {len(expression_lookup)}')

    hpa_genes = set()
    hpa_contexts = {}

    for source, gene_key, context_key in expression_lookup:
        hpa_genes.add(gene_key)

        if source not in hpa_contexts:
            hpa_contexts[source] = set()

        hpa_contexts[source].add(context_key)

    print('\n6. Поиск экспрессии для выбранных строк')

    expression_values = []
    expression_values_numeric = []
    expression_units = []
    expression_sources = []
    matched_contexts = []
    statuses = []
    checked_genes = []

    for _, row in work_df.iterrows():
        gene_key = row['_gene_key']
        checked_gene = row['_hpa_gene_symbol']
        raw_context = row[context_col]

        checked_genes.append(checked_gene)

        if gene_key == '':
            expression_values.append(pd.NA)
            expression_values_numeric.append(pd.NA)
            expression_units.append('')
            expression_sources.append('')
            matched_contexts.append('')
            statuses.append('empty geneSymbol')
            continue

        attempts = build_attempts(raw_context, aliases)
        hit = None

        for source, context_key in attempts:
            lookup_key = (source, gene_key, context_key)

            if lookup_key in expression_lookup:
                hit = expression_lookup[lookup_key]
                break

        if hit is not None:
            expression_values.append(hit['expression_value'])
            expression_values_numeric.append(hit['expression_value_numeric'])
            expression_units.append(hit['expression_unit'])
            expression_sources.append(hit['expression_source'])
            matched_contexts.append(hit['matched_hpa_context'])

            if hit['expression_value_numeric'] >= EXPRESSION_THRESHOLD:
                statuses.append(f'matched expression >= {EXPRESSION_THRESHOLD:g}')
            else:
                statuses.append(f'matched expression < {EXPRESSION_THRESHOLD:g}')
        else:
            expression_values.append(pd.NA)
            expression_values_numeric.append(pd.NA)
            expression_units.append('')
            expression_sources.append('')
            matched_contexts.append('')

            gene_exists = gene_key in hpa_genes
            context_exists = False

            for source, context_key in attempts:
                if source in hpa_contexts and context_key in hpa_contexts[source]:
                    context_exists = True
                    break

            if not gene_exists:
                statuses.append('Gene name not found in HPA')
            elif not context_exists:
                statuses.append('context not found in checked HPA sources')
            else:
                statuses.append('gene and context exist, but exact pair not found')

    work_df['expression_value'] = expression_values
    work_df['expression_unit'] = expression_units
    work_df['expression_source'] = expression_sources
    work_df['matched_hpa_context'] = matched_contexts
    work_df['_expression_value_numeric'] = expression_values_numeric
    work_df['_checked_hpa_gene_symbol'] = checked_genes
    work_df['_diagnostic_status'] = statuses

    matched_mask = work_df['_expression_value_numeric'].notna()
    expressed_mask = matched_mask & (work_df['_expression_value_numeric'] >= EXPRESSION_THRESHOLD)
    low_expression_mask = matched_mask & (work_df['_expression_value_numeric'] < EXPRESSION_THRESHOLD)
    unmatched_mask = ~matched_mask

    print(f'   Найдена экспрессия: {matched_mask.sum()}')
    print(f'   Оставлено в итоговой таблице, expression >= {EXPRESSION_THRESHOLD:g}: {expressed_mask.sum()}')
    print(f'   Найдена экспрессия, но ниже порога: {low_expression_mask.sum()}')
    print(f'   Не найдено данных об экспрессии: {unmatched_mask.sum()}')

    print('\n7. Сохранение результатов')

    helper_cols = [
        '_hpa_gene_symbol',
        '_gene_key',
        '_expression_value_numeric',
        '_checked_hpa_gene_symbol',
        '_diagnostic_status'
    ]

    output_df = work_df.loc[expressed_mask].drop(columns=helper_cols, errors='ignore').copy()
    output_df.to_csv(OUTPUT_FILE, sep='\t', index=False, na_rep='')
    print(f'   Итоговая таблица: {OUTPUT_FILE}')

    low_expression_df = work_df.loc[low_expression_mask].drop(columns=helper_cols, errors='ignore').copy()
    low_expression_df.to_csv(LOW_EXPRESSION_FILE, sep='\t', index=False, na_rep='')
    print(f'   Строки с экспрессией ниже порога: {LOW_EXPRESSION_FILE}')

    unmatched_df = work_df.loc[unmatched_mask, [
        gene_col,
        context_col,
        '_checked_hpa_gene_symbol',
        '_diagnostic_status'
    ]].copy()

    unmatched_df = unmatched_df.rename(
        columns={
            gene_col: 'geneSymbol',
            context_col: 'cell_line',
            '_checked_hpa_gene_symbol': 'Gene_symbol_checked_in_HPA',
            '_diagnostic_status': 'reason'
        }
    )

    unmatched_df.to_csv(UNMATCHED_FILE, sep='\t', index=False, na_rep='')
    print(f'   Несовпавшие строки: {UNMATCHED_FILE}')

    diagnostics_df = work_df[[
        gene_col,
        context_col,
        '_checked_hpa_gene_symbol',
        'expression_value',
        'expression_unit',
        'expression_source',
        'matched_hpa_context',
        '_diagnostic_status'
    ]].copy()

    diagnostics_df = diagnostics_df.rename(
        columns={
            gene_col: 'geneSymbol',
            context_col: 'cell_line',
            '_checked_hpa_gene_symbol': 'Gene_symbol_checked_in_HPA',
            '_diagnostic_status': 'diagnostic_status'
        }
    )

    diagnostics_df.to_csv(DIAGNOSTICS_FILE, sep='\t', index=False, na_rep='')
    print(f'   Диагностика: {DIAGNOSTICS_FILE}')

    print('\nГотово.')


if __name__ == '__main__':
    main()