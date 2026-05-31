#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


# ===================== НАСТРОЙКИ =====================

INPUT_FILE = Path('tf_snv_analysis_results/tf_snv_tissue_pioneer_class_family_filtered.tsv')
TF_LIST_FILE = Path('tf_list.tsv')

# nTPM-файлы
RNA_TISSUE_FILE = Path('HPA/rna_tissue_consensus.tsv')
RNA_CANCER_FILE = Path('HPA/rna_cell_line_cancer.tsv')
RNA_CELLINE_FILE = Path('HPA/rna_celline.tsv')

RNA_BRAIN_HPA_FILE = Path('HPA/rna_brain_hpa.tsv')
RNA_IMMUNE_CELL_FILE = Path('HPA/rna_immune_cell.tsv')
RNA_BRAIN_REGION_HPA_FILE = Path('HPA/rna_brain_region_hpa.tsv')
RNA_PFC_BRAIN_HPA_FILE = Path('HPA/rna_pfc_brain_hpa.tsv')

# nCPM-файлы
RNA_SINGLE_CELL_TYPE_FILE = Path('HPA/rna_single_cell_type.tsv')
RNA_SINGLE_CELL_TYPE_GROUP_FILE = Path('HPA/rna_single_cell_type_group.tsv')

# Необязательный файл дополнительных ручных соответствий контекстов.
# Если файла нет, скрипт работает только с ручными соответствиями,
# прописанными ниже в DEFAULT_CONTEXT_ALIASES.
CONTEXT_ALIASES_FILE = Path('context_aliases.tsv')

OUTPUT_FILE = Path('tf_snv_tissue_pioneer_class_family_filtered_with_expression.tsv')
UNMATCHED_FILE = Path('unmatched_hpa_expression.tsv')
DIAGNOSTICS_FILE = Path('hpa_expression_match_diagnostics.tsv')

# ======================================================


# Ручные соответствия контекстов.
# input_tissue — как записано в твоей таблице
# hpa_source — источник HPA в скрипте
# hpa_context — как записано в HPA
DEFAULT_CONTEXT_ALIASES = [
    {
        'input_tissue': 'A3_Jurkat_clone_A3__Childhood_T_acute_lymphoblastic_leukemia',
        'hpa_source': 'cell_line',
        'hpa_context': 'JURKAT'
    },
    {
        'input_tissue': 'HeLa_S3__cervical_adenocarcinoma',
        'hpa_source': 'cell_line',
        'hpa_context': 'HeLa'
    },
    {
        'input_tissue': 'TALL1__T-Acute_Lymphoblastic_Leukemia',
        'hpa_source': 'cell_line',
        'hpa_context': 'TALL-1 [Human adult T-ALL]'
    },
    {
        'input_tissue': 'heart',
        'hpa_source': 'tissue_consensus',
        'hpa_context': 'heart muscle'
    },
    {
        'input_tissue': 'periphery_of_retina',
        'hpa_source': 'tissue_consensus',
        'hpa_context': 'retina'
    }
]


# Ручные соответствия TF.
# input_tf — как TF записан в твоей таблице
# hpa_gene_symbol — как Gene name записан в HPA
DEFAULT_TF_ALIASES = [
    {
        'input_tf': 'AP2D',
        'hpa_gene_symbol': 'TFAP2D'
    },
    {
        'input_tf': 'BMAL1',
        'hpa_gene_symbol': 'BMAL1'
    },
    {
        'input_tf': 'ZF64A',
        'hpa_gene_symbol': 'ZFP64'
    }
]


def read_tsv(path):
    if not path.exists():
        raise FileNotFoundError(f'Не найден файл: {path}')

    df = pd.read_csv(path, sep='\t', dtype=str)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def clean_value(x):
    if pd.isna(x):
        return ''
    return str(x).strip()


def normalize_gene_symbol(x):
    return clean_value(x).upper()


def get_left_part_before_double_underscore(x):
    """
    Для строки:
    A673__Ewing_sarcoma

    возвращает:
    A673
    """
    x = clean_value(x)
    x = re.sub(r'_+$', '', x)

    parts = [p for p in re.split(r'__+', x) if clean_value(p) != '']

    if len(parts) >= 2:
        return clean_value(parts[0])

    return x


def get_right_part_after_double_underscore(x):
    """
    Для строки:
    A673__Ewing_sarcoma

    возвращает:
    Ewing_sarcoma
    """
    x = clean_value(x)
    x = re.sub(r'_+$', '', x)

    parts = [p for p in re.split(r'__+', x) if clean_value(p) != '']

    if len(parts) >= 2:
        return clean_value(parts[-1])

    return x


def normalize_general_context(x):
    """
    Общая нормализация для тканей, типов рака, иммунных клеток,
    регионов мозга и single-cell типов.

    Учитывает разные регистры:
    Heart, heart, HEART -> heart

    Разделители превращаются в пробелы:
    B-cells -> b cells
    adipose_tissue -> adipose tissue
    """
    x = clean_value(x)

    if x == '':
        return ''

    x = re.sub(r'_+$', '', x)

    x = x.replace('_', ' ')
    x = x.replace('-', ' ')
    x = x.replace('/', ' ')
    x = x.replace('\\', ' ')

    x = re.sub(r'[^A-Za-zА-Яа-я0-9]+', ' ', x)

    x = x.lower()
    x = re.sub(r'\s+', ' ', x).strip()

    return x


def normalize_cell_line_context(x):
    """
    Специальная нормализация для клеточных линий.

    Здесь разделители полностью удаляются:
    A673  -> a673
    A-673 -> a673
    A_673 -> a673
    22RV1 -> 22rv1
    22Rv1 -> 22rv1
    MCF7  -> mcf7
    MCF-7 -> mcf7
    """
    x = clean_value(x)

    if x == '':
        return ''

    x = re.sub(r'_+$', '', x)

    # Если передана полная строка вида MCF7__Invasive_ductal_breast_carcinoma,
    # для клеточной линии берём только левую часть.
    x = get_left_part_before_double_underscore(x)

    x = x.lower()

    # Оставляем только буквы и цифры.
    x = re.sub(r'[^a-z0-9]+', '', x)

    return x


def add_unique_variant(variants, value):
    value = clean_value(value)

    if value != '' and value not in variants:
        variants.append(value)


def get_context_variants(x, source):
    """
    Возвращает список ключей для поиска.

    Для cell_line используется специальная нормализация.
    Для остальных источников используется общая нормализация.

    Для single-cell источников добавлены дополнительные варианты:
    trophoblast cell <-> trophoblast cells
    adipocyte <-> adipocytes
    fibroblast -> fibroblasts
    pancreatic islets -> pancreatic islet cells
    """
    variants = []

    if source == 'cell_line':
        base = normalize_cell_line_context(x)

        if base != '':
            add_unique_variant(variants, base)

    else:
        base = normalize_general_context(x)

        if base != '':
            add_unique_variant(variants, base)

        # cell / cells
        if ' cells' in base:
            add_unique_variant(variants, base.replace(' cells', ' cell'))

        if ' cell' in base and ' cells' not in base:
            add_unique_variant(variants, base.replace(' cell', ' cells'))

        # lymphocyte / lymphocytes
        if ' lymphocytes' in base:
            add_unique_variant(variants, base.replace(' lymphocytes', ' lymphocyte'))

        if ' lymphocyte' in base and ' lymphocytes' not in base:
            add_unique_variant(variants, base.replace(' lymphocyte', ' lymphocytes'))

        # Общие варианты для single-cell данных.
        if source in ['single_cell_type', 'single_cell_type_group']:
            if base.endswith(' fibroblast') or ' fibroblast' in base:
                add_unique_variant(variants, 'fibroblasts')

            if base == 'skin fibroblast':
                add_unique_variant(variants, 'fibroblasts')

            if base == 'foreskin fibroblast':
                add_unique_variant(variants, 'fibroblasts')

            if base == 'human secondary fibroblast':
                add_unique_variant(variants, 'fibroblasts')

            if base == 'fibroblast of gingiva':
                add_unique_variant(variants, 'fibroblasts')

            if base == 'fibroblast like synoviocytes':
                add_unique_variant(variants, 'fibroblasts')

            if base == 'adipocyte':
                add_unique_variant(variants, 'adipocytes')

            if base == 'primary white adipocytes':
                add_unique_variant(variants, 'adipocytes')

            if base == 'hepatocyte':
                add_unique_variant(variants, 'hepatocytes')

            if base == 'hepatocytes':
                add_unique_variant(variants, 'hepatocyte')

            if base == 'myocyte':
                add_unique_variant(variants, 'myocytes')

            if base == 'myoblast':
                add_unique_variant(variants, 'myoblasts')

            if base == 'myoblasts':
                add_unique_variant(variants, 'myoblast')

            if base == 'smooth muscle cell':
                add_unique_variant(variants, 'smooth muscle cells')

            if base == 'skeletal muscle cell':
                add_unique_variant(variants, 'skeletal myocytes')
                add_unique_variant(variants, 'myocytes')

            if base == 'pancreatic islets':
                add_unique_variant(variants, 'pancreatic islet cells')

            if base == 'pancreatic islet':
                add_unique_variant(variants, 'pancreatic islet cells')

            if base == 'islet precursor cell':
                add_unique_variant(variants, 'pancreatic islet cells')

            if base == 'kidney tubule cell':
                add_unique_variant(variants, 'proximal tubular cells')
                add_unique_variant(variants, 'distal tubular cells')
                add_unique_variant(variants, 'renal tubular cells')
                add_unique_variant(variants, 'renal nephron cells')

            if base == 'glomerular visceral epithelial cell':
                add_unique_variant(variants, 'podocytes')

            if base == 'trophoblast cell':
                add_unique_variant(variants, 'trophoblast cells')

            if base == 'hematopoietic multipotent progenitor cell':
                add_unique_variant(variants, 'hematopoietic stem cells')

            if base == 'cd34 hematopoietic stem cells':
                add_unique_variant(variants, 'hematopoietic stem cells')

            if base == 'bipolar neuron':
                add_unique_variant(variants, 'retinal bipolar cells')

            if base == 'neural progenitors':
                add_unique_variant(variants, 'neural progenitor cells')

            if base == 'ecto neural progenitor cell':
                add_unique_variant(variants, 'neural progenitor cells')

            if base == 'endodermal cell':
                add_unique_variant(variants, 'endodermal cells')

            if base == 'chondrocytes':
                add_unique_variant(variants, 'chondrocyte')

            if base == 'ips derived chondrocytes':
                add_unique_variant(variants, 'chondrocytes')
                add_unique_variant(variants, 'chondrocyte')

            if base == 'foreskin keratinocyte':
                add_unique_variant(variants, 'keratinocytes')
                add_unique_variant(variants, 'keratinocyte')

            if base == 'erythroid blood cells':
                add_unique_variant(variants, 'erythrocytes')

            if base == 'plasmablasts':
                add_unique_variant(variants, 'plasma cells')

            if base == 'peripheral blood derived mast cells':
                add_unique_variant(variants, 'mast cells')

            if base == 'monocyte derived macrophages from peripheral blood':
                add_unique_variant(variants, 'macrophages')

        # Иммунные варианты.
        if source in ['immune_cell', 'single_cell_type', 'single_cell_type_group']:
            if base == 'b cells':
                add_unique_variant(variants, 'b cells')
                add_unique_variant(variants, 'b-cells')
                add_unique_variant(variants, 'b cell')
                add_unique_variant(variants, 'b-cell')

            if base == 'bulk b cells':
                add_unique_variant(variants, 'b cells')
                add_unique_variant(variants, 'b-cells')

            if base == 'primary b cells':
                add_unique_variant(variants, 'b cells')
                add_unique_variant(variants, 'b-cells')

            if base == 'cd19 b cells':
                add_unique_variant(variants, 'b cells')
                add_unique_variant(variants, 'b-cells')

            if base == 'cd19cd5 b cell':
                add_unique_variant(variants, 'b cells')
                add_unique_variant(variants, 'b-cells')

            if base == 't cells':
                add_unique_variant(variants, 't cells')
                add_unique_variant(variants, 't-cells')
                add_unique_variant(variants, 't cell')
                add_unique_variant(variants, 't-cell')

            if base == 'cd8 t cells':
                add_unique_variant(variants, 'cd8 t cells')
                add_unique_variant(variants, 'naive cd8 t-cell')
                add_unique_variant(variants, 'memory cd8 t-cell')
                add_unique_variant(variants, 'cytotoxic t cells')

            if base == 'cd8 lymphocytes':
                add_unique_variant(variants, 'cd8 t cells')
                add_unique_variant(variants, 'naive cd8 t-cell')
                add_unique_variant(variants, 'memory cd8 t-cell')

            if base == 'cd4 t cells':
                add_unique_variant(variants, 'cd4 t cells')
                add_unique_variant(variants, 'naive cd4 t-cell')
                add_unique_variant(variants, 'memory cd4 t-cell')
                add_unique_variant(variants, 'helper t cells')

            if base == 'primary cd4 t cells':
                add_unique_variant(variants, 'cd4 t cells')
                add_unique_variant(variants, 'naive cd4 t-cell')
                add_unique_variant(variants, 'memory cd4 t-cell')

            if base == 'cd14 monocytes':
                add_unique_variant(variants, 'monocytes')
                add_unique_variant(variants, 'classical monocyte')
                add_unique_variant(variants, 'classical monocytes')

            if base == 'monocytes':
                add_unique_variant(variants, 'monocyte')
                add_unique_variant(variants, 'classical monocyte')
                add_unique_variant(variants, 'classical monocytes')

            if base == 'mdm monocyte derived macrophages':
                add_unique_variant(variants, 'macrophages')

            if base == 'regulatory t cell':
                add_unique_variant(variants, 't reg')
                add_unique_variant(variants, 't-reg')
                add_unique_variant(variants, 'regulatory t cells')

            if base == 'memory t regulatory':
                add_unique_variant(variants, 't-reg')
                add_unique_variant(variants, 'regulatory t cells')

            if base == 'gamma delta t cells':
                add_unique_variant(variants, 'gdt cell')
                add_unique_variant(variants, 'gdt-cell')
                add_unique_variant(variants, 'gamma delta t cell')

            if base == 'natural killer cell':
                add_unique_variant(variants, 'nk cell')
                add_unique_variant(variants, 'nk-cell')
                add_unique_variant(variants, 'nk cells')
                add_unique_variant(variants, 'nk-cells')

            if base == 'mature natural killer':
                add_unique_variant(variants, 'nk cell')
                add_unique_variant(variants, 'nk-cell')
                add_unique_variant(variants, 'nk cells')
                add_unique_variant(variants, 'nk-cells')

    return variants


def find_first_column(df, possible_names, table_name):
    for col in possible_names:
        if col in df.columns:
            return col

    raise ValueError(
        f'В таблице {table_name} не найден ни один из столбцов: {possible_names}. '
        f'Доступные столбцы: {list(df.columns)}'
    )


def load_hpa_expression(path, context_cols, value_col, source_name, expression_unit):
    df = read_tsv(path)

    if isinstance(context_cols, str):
        context_cols = [context_cols]

    context_col = find_first_column(df, context_cols, path.name)

    required_cols = ['Gene name', context_col, value_col]
    missing = [c for c in required_cols if c not in df.columns]

    if len(missing) > 0:
        raise ValueError(
            f'В файле {path} отсутствуют столбцы: {missing}. '
            f'Доступные столбцы: {list(df.columns)}'
        )

    df = df[['Gene name', context_col, value_col]].copy()
    df = df.rename(
        columns={
            context_col: 'HPA_context',
            value_col: 'expression_value'
        }
    )

    df['Gene name'] = df['Gene name'].apply(clean_value)
    df['HPA_context'] = df['HPA_context'].apply(clean_value)
    df['expression_value'] = df['expression_value'].apply(clean_value)

    df = df[
        (df['Gene name'] != '') &
        (df['HPA_context'] != '') &
        (df['expression_value'] != '')
    ].copy()

    df['expression_unit'] = expression_unit
    df['source'] = source_name
    df['gene_key'] = df['Gene name'].apply(normalize_gene_symbol)

    return df


def add_attempt(attempts, seen, source, context):
    context_variants = get_context_variants(context, source)

    for context_key in context_variants:
        key = (source, context_key)

        if key not in seen:
            attempts.append((source, context, context_key))
            seen.add(key)


def load_context_aliases():
    """
    Загружает ручные соответствия контекстов.

    Сначала добавляются соответствия, прописанные прямо в скрипте
    в DEFAULT_CONTEXT_ALIASES.

    Дополнительно, если существует файл context_aliases.tsv,
    будут загружены соответствия из него.

    Формат context_aliases.tsv:

    input_tissue    hpa_source    hpa_context
    heart           tissue_consensus    heart muscle
    skin_fibroblast single_cell_type_group    fibroblasts
    """
    aliases = {}

    def add_alias(input_tissue, hpa_source, hpa_context):
        input_tissue = clean_value(input_tissue)
        hpa_source = clean_value(hpa_source)
        hpa_context = clean_value(hpa_context)

        if input_tissue == '' or hpa_source == '' or hpa_context == '':
            return

        key_exact = input_tissue
        key_norm = normalize_general_context(input_tissue)

        value = {
            'source': hpa_source,
            'context': hpa_context
        }

        if key_exact not in aliases:
            aliases[key_exact] = []

        aliases[key_exact].append(value)

        if key_norm not in aliases:
            aliases[key_norm] = []

        aliases[key_norm].append(value)

    for item in DEFAULT_CONTEXT_ALIASES:
        add_alias(
            item['input_tissue'],
            item['hpa_source'],
            item['hpa_context']
        )

    print(f'   Загружено ручных соответствий контекстов из скрипта: {len(DEFAULT_CONTEXT_ALIASES)}')

    if not CONTEXT_ALIASES_FILE.exists():
        print(f'   Файл дополнительных ручных соответствий не найден: {CONTEXT_ALIASES_FILE}')
        return aliases

    alias_df = read_tsv(CONTEXT_ALIASES_FILE)

    required_cols = ['input_tissue', 'hpa_source', 'hpa_context']
    missing = [c for c in required_cols if c not in alias_df.columns]

    if len(missing) > 0:
        raise ValueError(
            f'В файле {CONTEXT_ALIASES_FILE} отсутствуют столбцы: {missing}. '
            f'Доступные столбцы: {list(alias_df.columns)}'
        )

    alias_df = alias_df[required_cols].copy()

    alias_df['input_tissue'] = alias_df['input_tissue'].apply(clean_value)
    alias_df['hpa_source'] = alias_df['hpa_source'].apply(clean_value)
    alias_df['hpa_context'] = alias_df['hpa_context'].apply(clean_value)

    alias_df = alias_df[
        (alias_df['input_tissue'] != '') &
        (alias_df['hpa_source'] != '') &
        (alias_df['hpa_context'] != '')
    ].copy()

    for _, row in alias_df.iterrows():
        add_alias(
            row['input_tissue'],
            row['hpa_source'],
            row['hpa_context']
        )

    print(f'   Загружено дополнительных ручных соответствий контекстов из файла: {len(alias_df)}')

    return aliases


def add_alias_attempts(attempts, seen, raw_context, aliases):
    raw_context = clean_value(raw_context)
    norm_context = normalize_general_context(raw_context)

    alias_items = []

    if raw_context in aliases:
        alias_items.extend(aliases[raw_context])

    if norm_context in aliases:
        alias_items.extend(aliases[norm_context])

    for item in alias_items:
        add_attempt(attempts, seen, item['source'], item['context'])


def build_match_attempts(raw_context, aliases):
    raw_context = clean_value(raw_context)
    raw_context = re.sub(r'_+$', '', raw_context)

    attempts = []
    seen = set()

    # Сначала ручные соответствия контекстов.
    add_alias_attempts(attempts, seen, raw_context, aliases)

    cell_line_part = get_left_part_before_double_underscore(raw_context)
    cancer_part = get_right_part_after_double_underscore(raw_context)

    has_double_underscore = '__' in raw_context

    if has_double_underscore:
        # Для строк вида MCF7__Invasive_ductal_breast_carcinoma
        # в cell_line проверяется только MCF7.
        add_attempt(attempts, seen, 'cell_line', cell_line_part)

        # Дополнительно пробуем всю строку как cell_line,
        # но normalize_cell_line_context всё равно возьмёт левую часть.
        add_attempt(attempts, seen, 'cell_line', raw_context)

        # Для cancer-таблицы проверяем правую часть.
        add_attempt(attempts, seen, 'cell_line_cancer', cancer_part)

        # Резервные попытки для nTPM-источников.
        add_attempt(attempts, seen, 'tissue_consensus', cancer_part)
        add_attempt(attempts, seen, 'immune_cell', cell_line_part)
        add_attempt(attempts, seen, 'brain_pfc_subregion', cancer_part)
        add_attempt(attempts, seen, 'brain_subregion', cancer_part)
        add_attempt(attempts, seen, 'brain_region', cancer_part)

        # Резервные попытки для nCPM single-cell источников.
        add_attempt(attempts, seen, 'single_cell_type', cell_line_part)
        add_attempt(attempts, seen, 'single_cell_type_group', cell_line_part)
        add_attempt(attempts, seen, 'single_cell_type', cancer_part)
        add_attempt(attempts, seen, 'single_cell_type_group', cancer_part)

    else:
        # Более специфичные nTPM-источники.
        add_attempt(attempts, seen, 'cell_line', raw_context)
        add_attempt(attempts, seen, 'immune_cell', raw_context)
        add_attempt(attempts, seen, 'brain_pfc_subregion', raw_context)
        add_attempt(attempts, seen, 'brain_subregion', raw_context)
        add_attempt(attempts, seen, 'brain_region', raw_context)
        add_attempt(attempts, seen, 'tissue_consensus', raw_context)
        add_attempt(attempts, seen, 'cell_line_cancer', raw_context)

        # nCPM single-cell источники.
        add_attempt(attempts, seen, 'single_cell_type', raw_context)
        add_attempt(attempts, seen, 'single_cell_type_group', raw_context)

    return attempts


def add_expression_to_lookup(expression_lookup, expr_df):
    for _, row in expr_df.iterrows():
        gene_key = row['gene_key']
        source = row['source']
        hpa_context = row['HPA_context']

        context_keys = get_context_variants(hpa_context, source)

        for context_key in context_keys:
            key = (source, gene_key, context_key)

            if key not in expression_lookup:
                expression_lookup[key] = {
                    'expression_value': row['expression_value'],
                    'expression_unit': row['expression_unit'],
                    'expression_source': source,
                    'hpa_context': hpa_context
                }


def is_zero_expression(value):
    value = clean_value(value)

    if value == '':
        return False

    try:
        return float(value) == 0.0
    except ValueError:
        return value in ['0', '0.0', '0.00', '0.000']


def add_manual_tf_aliases(tf_map):
    """
    Добавляет ручные соответствия TF.

    В tf_map:
    ключ — как TF записан в твоей таблице;
    значение — как Gene name записан в HPA.
    """
    for item in DEFAULT_TF_ALIASES:
        input_tf = clean_value(item['input_tf'])
        hpa_gene_symbol = clean_value(item['hpa_gene_symbol'])

        if input_tf != '' and hpa_gene_symbol != '':
            tf_map[input_tf] = hpa_gene_symbol

    return tf_map


def main():
    print('1. Загрузка основной таблицы...')

    main_df = read_tsv(INPUT_FILE)

    tf_col = find_first_column(
        main_df,
        ['TF', 'Gene name'],
        INPUT_FILE.name
    )

    tissue_col = find_first_column(
        main_df,
        ['tissue', 'Tissue'],
        INPUT_FILE.name
    )

    print(f'   Основная таблица: {INPUT_FILE}')
    print(f'   Строк: {len(main_df)}')
    print(f'   Столбец TF: {tf_col}')
    print(f'   Столбец tissue: {tissue_col}')

    print('\n2. Загрузка словаря TF...')

    tf_dict = read_tsv(TF_LIST_FILE)

    required_dict_cols = ['curated:uniprot_id', 'auto:gene_symbol']
    missing_dict_cols = [c for c in required_dict_cols if c not in tf_dict.columns]

    if len(missing_dict_cols) > 0:
        raise ValueError(
            f'В файле {TF_LIST_FILE} отсутствуют столбцы: {missing_dict_cols}. '
            f'Доступные столбцы: {list(tf_dict.columns)}'
        )

    tf_dict = tf_dict[required_dict_cols].copy()
    tf_dict['curated:uniprot_id'] = tf_dict['curated:uniprot_id'].apply(clean_value)
    tf_dict['auto:gene_symbol'] = tf_dict['auto:gene_symbol'].apply(clean_value)

    tf_dict = tf_dict[
        (tf_dict['curated:uniprot_id'] != '') &
        (tf_dict['auto:gene_symbol'] != '')
    ].copy()

    duplicated_tf = tf_dict['curated:uniprot_id'].duplicated(keep=False).sum()

    if duplicated_tf > 0:
        print(f'   Предупреждение: повторяющихся curated-названий в tf_list.tsv: {duplicated_tf}')
        print('   Для повторов будет использовано первое встретившееся соответствие.')

    tf_dict = tf_dict.drop_duplicates(subset=['curated:uniprot_id'], keep='first')

    # Основной словарь из tf_list.tsv:
    # ключ — как TF записан в твоей таблице;
    # значение — как Gene name записан в HPA.
    tf_map = dict(zip(tf_dict['curated:uniprot_id'], tf_dict['auto:gene_symbol']))

    print(f'   Соответствий в словаре tf_list.tsv: {len(tf_map)}')

    tf_map = add_manual_tf_aliases(tf_map)

    print(f'   Ручных соответствий TF добавлено: {len(DEFAULT_TF_ALIASES)}')
    print(f'   Всего соответствий TF после добавления ручных замен: {len(tf_map)}')

    main_df['_tf_key'] = main_df[tf_col].apply(clean_value)
    main_df['_hpa_gene_symbol'] = main_df['_tf_key'].map(tf_map)

    rows_with_tf_map = main_df['_hpa_gene_symbol'].notna().sum()
    rows_without_tf_map = len(main_df) - rows_with_tf_map

    print(f'   Строк с TF, найденным в словаре: {rows_with_tf_map}')
    print(f'   Строк с TF, не найденным в словаре: {rows_without_tf_map}')

    print('\n3. Загрузка ручных соответствий контекстов...')

    aliases = load_context_aliases()

    print('\n4. Загрузка таблиц HPA...')

    tissue_expr = load_hpa_expression(
        RNA_TISSUE_FILE,
        ['Tissue'],
        'nTPM',
        'tissue_consensus',
        'nTPM'
    )

    cancer_expr = load_hpa_expression(
        RNA_CANCER_FILE,
        ['Cancer'],
        'nTPM',
        'cell_line_cancer',
        'nTPM'
    )

    cellline_expr = load_hpa_expression(
        RNA_CELLINE_FILE,
        ['Cell line'],
        'nTPM',
        'cell_line',
        'nTPM'
    )

    brain_expr = load_hpa_expression(
        RNA_BRAIN_HPA_FILE,
        ['Subregion'],
        'nTPM',
        'brain_subregion',
        'nTPM'
    )

    immune_expr = load_hpa_expression(
        RNA_IMMUNE_CELL_FILE,
        ['Immune cell'],
        'nTPM',
        'immune_cell',
        'nTPM'
    )

    brain_region_expr = load_hpa_expression(
        RNA_BRAIN_REGION_HPA_FILE,
        ['Brain region'],
        'nTPM',
        'brain_region',
        'nTPM'
    )

    pfc_brain_expr = load_hpa_expression(
        RNA_PFC_BRAIN_HPA_FILE,
        ['Subregion'],
        'nTPM',
        'brain_pfc_subregion',
        'nTPM'
    )

    single_cell_type_expr = load_hpa_expression(
        RNA_SINGLE_CELL_TYPE_FILE,
        ['Cell type'],
        'nCPM',
        'single_cell_type',
        'nCPM'
    )

    single_cell_type_group_expr = load_hpa_expression(
        RNA_SINGLE_CELL_TYPE_GROUP_FILE,
        ['Cell type group', 'Cell type'],
        'nCPM',
        'single_cell_type_group',
        'nCPM'
    )

    print(f'   {RNA_TISSUE_FILE}: {len(tissue_expr)} строк после фильтрации')
    print(f'   {RNA_CANCER_FILE}: {len(cancer_expr)} строк после фильтрации')
    print(f'   {RNA_CELLINE_FILE}: {len(cellline_expr)} строк после фильтрации')
    print(f'   {RNA_BRAIN_HPA_FILE}: {len(brain_expr)} строк после фильтрации')
    print(f'   {RNA_IMMUNE_CELL_FILE}: {len(immune_expr)} строк после фильтрации')
    print(f'   {RNA_BRAIN_REGION_HPA_FILE}: {len(brain_region_expr)} строк после фильтрации')
    print(f'   {RNA_PFC_BRAIN_HPA_FILE}: {len(pfc_brain_expr)} строк после фильтрации')
    print(f'   {RNA_SINGLE_CELL_TYPE_FILE}: {len(single_cell_type_expr)} строк после фильтрации')
    print(f'   {RNA_SINGLE_CELL_TYPE_GROUP_FILE}: {len(single_cell_type_group_expr)} строк после фильтрации')

    print('\n5. Создание индекса HPA...')

    expression_lookup = {}

    # Порядок важен: более точные источники добавляются раньше.
    # Если один и тот же ключ встретится повторно, сохранится первое значение.
    add_expression_to_lookup(expression_lookup, cellline_expr)
    add_expression_to_lookup(expression_lookup, immune_expr)
    add_expression_to_lookup(expression_lookup, pfc_brain_expr)
    add_expression_to_lookup(expression_lookup, brain_expr)
    add_expression_to_lookup(expression_lookup, brain_region_expr)
    add_expression_to_lookup(expression_lookup, tissue_expr)
    add_expression_to_lookup(expression_lookup, cancer_expr)

    # Single-cell nCPM добавляется после nTPM-источников,
    # чтобы не заменять более прямое bulk/cell line совпадение.
    add_expression_to_lookup(expression_lookup, single_cell_type_expr)
    add_expression_to_lookup(expression_lookup, single_cell_type_group_expr)

    print(f'   Уникальных ключей Gene name × context × source в HPA: {len(expression_lookup)}')

    all_hpa_genes = set()
    all_hpa_contexts_by_source = {}

    for key in expression_lookup:
        source, gene_key, context_key = key

        all_hpa_genes.add(gene_key)

        if source not in all_hpa_contexts_by_source:
            all_hpa_contexts_by_source[source] = set()

        all_hpa_contexts_by_source[source].add(context_key)

    print('\n6. Быстрая проверка нормализации:')

    test_contexts = [
        'A673',
        'A-673',
        'A-172',
        '22RV1',
        '22Rv1',
        'MCF7',
        'Heart',
        'heart',
        'HEART',
        'skin_fibroblast',
        'trophoblast_cell',
        'pancreatic_islets',
        'A3_Jurkat_clone_A3__Childhood_T_acute_lymphoblastic_leukemia',
        'HeLa_S3__cervical_adenocarcinoma',
        'TALL1__T-Acute_Lymphoblastic_Leukemia',
        'periphery_of_retina'
    ]

    for context in test_contexts:
        print(
            f'   {context} -> general: {normalize_general_context(context)} | '
            f'cell_line: {normalize_cell_line_context(context)}'
        )

    print('\n7. Быстрая проверка ручных TF-соответствий:')

    for item in DEFAULT_TF_ALIASES:
        input_tf = item['input_tf']
        hpa_gene_symbol = item['hpa_gene_symbol']
        print(f'   {input_tf} -> {hpa_gene_symbol}')

    print('\n8. Сопоставление основной таблицы с HPA...')

    expression_values = []
    expression_units = []
    expression_sources = []
    matched_hpa_contexts = []
    diagnostic_statuses = []
    checked_gene_symbols = []
    checked_contexts = []

    for _, row in main_df.iterrows():
        hpa_gene = clean_value(row['_hpa_gene_symbol'])
        raw_context = clean_value(row[tissue_col])

        checked_gene_symbols.append(hpa_gene)
        checked_contexts.append(raw_context)

        if hpa_gene == '':
            expression_values.append(pd.NA)
            expression_units.append('')
            expression_sources.append('')
            matched_hpa_contexts.append('')
            diagnostic_statuses.append('TF not found in tf_list.tsv')
            continue

        gene_key = normalize_gene_symbol(hpa_gene)
        attempts = build_match_attempts(raw_context, aliases)

        found = False

        for source, original_context, context_key in attempts:
            lookup_key = (source, gene_key, context_key)

            if lookup_key in expression_lookup:
                hit = expression_lookup[lookup_key]

                expression_values.append(hit['expression_value'])
                expression_units.append(hit['expression_unit'])
                expression_sources.append(hit['expression_source'])
                matched_hpa_contexts.append(hit['hpa_context'])

                if is_zero_expression(hit['expression_value']):
                    diagnostic_statuses.append(
                        f'matched: {hit["expression_unit"]} is zero'
                    )
                else:
                    diagnostic_statuses.append(
                        f'matched: {hit["expression_unit"]} found'
                    )

                found = True
                break

        if not found:
            expression_values.append(pd.NA)
            expression_units.append('')
            expression_sources.append('')
            matched_hpa_contexts.append('')

            gene_exists = gene_key in all_hpa_genes

            context_exists = False

            for source, original_context, context_key in attempts:
                if source in all_hpa_contexts_by_source:
                    if context_key in all_hpa_contexts_by_source[source]:
                        context_exists = True
                        break

            if not gene_exists:
                diagnostic_statuses.append(
                    'Gene symbol from tf_list.tsv not found in HPA expression tables'
                )
            elif not context_exists:
                diagnostic_statuses.append(
                    'Context name not found in checked HPA sources'
                )
            else:
                diagnostic_statuses.append(
                    'Gene exists and context exists, but exact Gene name × context pair was not found'
                )

    main_df['expression_value'] = expression_values
    main_df['expression_unit'] = expression_units
    main_df['expression_source'] = expression_sources

    main_df['_matched_hpa_context'] = matched_hpa_contexts
    main_df['_diagnostic_status'] = diagnostic_statuses
    main_df['_checked_hpa_gene_symbol'] = checked_gene_symbols
    main_df['_checked_context'] = checked_contexts

    matched_mask = main_df['expression_value'].notna()
    n_matched = matched_mask.sum()
    n_unmatched = len(main_df) - n_matched

    print(f'   Совпадений найдено: {n_matched}')
    print(f'   Совпадений не найдено: {n_unmatched}')

    print('\n9. Совпадения по источникам HPA:')

    source_counts = main_df.loc[matched_mask, 'expression_source'].value_counts()

    if len(source_counts) == 0:
        print('   Совпадений по HPA-таблицам нет.')
    else:
        for source, count in source_counts.items():
            print(f'   {source}: {count}')

    print('\n10. Совпадения по единицам экспрессии:')

    unit_counts = main_df.loc[matched_mask, 'expression_unit'].value_counts()

    if len(unit_counts) == 0:
        print('   Совпадений нет.')
    else:
        for unit, count in unit_counts.items():
            print(f'   {unit}: {count}')

    print('\n11. Диагностика сопоставления:')

    diagnostic_counts = main_df['_diagnostic_status'].value_counts(dropna=False)

    for status, count in diagnostic_counts.items():
        print(f'   {status}: {count}')

    print('\n12. Сохранение результата...')

    output_df = main_df.drop(
        columns=[
            '_tf_key',
            '_hpa_gene_symbol',
            '_matched_hpa_context',
            '_diagnostic_status',
            '_checked_hpa_gene_symbol',
            '_checked_context'
        ],
        errors='ignore'
    )

    output_df.to_csv(
        OUTPUT_FILE,
        sep='\t',
        index=False,
        na_rep=''
    )

    print(f'   Итоговая таблица сохранена: {OUTPUT_FILE}')

    unmatched_df = main_df.loc[~matched_mask, [
        tf_col,
        tissue_col,
        '_checked_hpa_gene_symbol',
        '_diagnostic_status'
    ]].copy()

    unmatched_df = unmatched_df.rename(
        columns={
            tf_col: 'TF',
            tissue_col: 'tissue',
            '_checked_hpa_gene_symbol': 'Gene_symbol_for_HPA',
            '_diagnostic_status': 'reason'
        }
    )

    unmatched_df.to_csv(
        UNMATCHED_FILE,
        sep='\t',
        index=False,
        na_rep=''
    )

    print(f'   Таблица несовпавших строк сохранена: {UNMATCHED_FILE}')

    diagnostics_df = main_df[[
        tf_col,
        tissue_col,
        '_checked_hpa_gene_symbol',
        'expression_value',
        'expression_unit',
        'expression_source',
        '_matched_hpa_context',
        '_diagnostic_status'
    ]].copy()

    diagnostics_df = diagnostics_df.rename(
        columns={
            tf_col: 'TF',
            tissue_col: 'tissue',
            '_checked_hpa_gene_symbol': 'Gene_symbol_checked_in_HPA',
            '_matched_hpa_context': 'matched_hpa_context',
            '_diagnostic_status': 'diagnostic_status'
        }
    )

    diagnostics_df.to_csv(
        DIAGNOSTICS_FILE,
        sep='\t',
        index=False,
        na_rep=''
    )

    print(f'   Диагностическая таблица сохранена: {DIAGNOSTICS_FILE}')

    print('\nГотово.')


if __name__ == '__main__':
    main()