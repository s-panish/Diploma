import pandas as pd

# === Пути к файлам ===
MASTER_FILE      = 'tf_masterlist.tsv'
ENRICHMENT_FILE  = "tf_snv_analysis_results/tf_snv_tissue_pioneer_with_KN.tsv"
OUTPUT_FILE      = "tf_snv_analysis_results/tf_snv_tissue_pioneer_class_family.tsv"

print("Загрузка файлов...")
df_master = pd.read_csv(MASTER_FILE, sep='\t', dtype=str)
df_enrich = pd.read_csv(ENRICHMENT_FILE, sep='\t', dtype=str)

# 1. Оставляем только Human записи
df_human = df_master[df_master['curated:uniprot_id'].str.contains('_HUMAN', na=False, regex=False)].copy()
print(f"   Найдено {len(df_human)} записей для Human.")

# 2. Строим словари маппинга из ОБОИХ столбцов
class_dict = {}
family_dict = {}

for _, row in df_human.iterrows():
    # Безопасная очистка от NaN и пробелов
    cls = str(row['tfclass:class']).strip() if pd.notna(row['tfclass:class']) else ''
    fam = str(row['tfclass:family']).strip() if pd.notna(row['tfclass:family']) else ''
    cls_id = str(row['tfclass:id']).strip() if pd.notna(row['tfclass:id']) else ''

    keys_to_add = []

    # Обрабатываем auto:gene_symbol (может содержать несколько символов через / или ;)
    raw_sym = str(row['auto:gene_symbol']).strip()
    if raw_sym.upper() not in ('NAN', 'NA', 'NONE', ''):
        for sym in raw_sym.replace(';', '/').split('/'):
            key = sym.strip().upper()
            if key:
                keys_to_add.append(key)

    # Обрабатываем curated:uniprot_id (убираем суффикс _HUMAN)
    raw_up = str(row['curated:uniprot_id']).strip()
    if raw_up.upper() not in ('NAN', 'NA', 'NONE', ''):
        key = raw_up.replace('_HUMAN', '').strip().upper()
        if key:
            keys_to_add.append(key)

    # Форматируем класс: добавляем ID в фигурных скобках, если он присутствует
    if cls_id and cls_id.upper() not in ('NAN', 'NA', 'NONE', ''):
        formatted_cls = f"{cls} {{{cls_id}}}"
    else:
        formatted_cls = cls

    # Добавляем в словари (сохраняем первое найденное значение)
    for k in keys_to_add:
        if k not in class_dict:
            class_dict[k] = formatted_cls
            family_dict[k] = fam

print(f"   Словарь маппинга создан для {len(class_dict)} уникальных идентификаторов.")

# 3. Нормализуем колонку TF в таблице обогащения
tf_norm = df_enrich['TF'].str.strip().str.upper()

# 4. Применяем маппинг к двум ОТДЕЛЬНЫМ колонкам
df_enrich['tfclass:class']  = tf_norm.map(class_dict)
df_enrich['tfclass:family'] = tf_norm.map(family_dict)

# Заполняем ненайденные TF меткой Not_found (pandas 3.0 совместимо)
df_enrich['tfclass:class'] = df_enrich['tfclass:class'].fillna('Not_found')
df_enrich['tfclass:family'] = df_enrich['tfclass:family'].fillna('Not_found')

# 5. Сохраняем результат
df_enrich.to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f"\nГотово! Файл сохранён: {OUTPUT_FILE}")
not_found = (df_enrich['tfclass:class'] == 'Not_found').sum()
print(f"TF без найденного класса/семейства: {not_found} из {len(df_enrich)}")