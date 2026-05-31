import pandas as pd
import os
import glob

MOTIFBREAKR_FILE = 'motifbreakr_filtered_001_log.tsv'
TF_DIR = 'TF'
OUTPUT_FILE = 'motifbr+adastra_001_log.tsv'

# --- 1. Загрузка таблицы MotifBreakr ---
print(f"Чтение файла: {MOTIFBREAKR_FILE}...")
if not os.path.exists(MOTIFBREAKR_FILE):
    raise FileNotFoundError(f"Файл {MOTIFBREAKR_FILE} не найден.")

df_mb = pd.read_csv(MOTIFBREAKR_FILE, sep='\t', low_memory=False)

required_mb_cols = ['SNP_id', 'geneSymbol']
for col in required_mb_cols:
    if col not in df_mb.columns:
        raise ValueError(f"В файле {MOTIFBREAKR_FILE} отсутствует обязательный столбец: {col}")

# --- 2. Загрузка и обработка таблиц из папки TF ---
print(f"Чтение файлов из директории: {TF_DIR}/...")

if not os.path.exists(TF_DIR):
    raise FileNotFoundError(f"Директория {TF_DIR} не найдена.")

tf_files = glob.glob(os.path.join(TF_DIR, '*_HUMAN.tsv'))

if not tf_files:
    raise FileNotFoundError(f"В директории {TF_DIR} не найдено файлов по маске *_HUMAN.tsv")

df_tf_list = []

for file_path in tf_files:
    filename = os.path.basename(file_path)
    tf_name = filename.replace('_HUMAN.tsv', '')
    
    try:
        df_temp = pd.read_csv(file_path, sep='\t', comment=None, low_memory=False)
        df_temp['geneSymbol'] = tf_name
        rename_map = {
            '#chr': 'tf_chr',
            'start': 'tf_start',
            'end': 'tf_end',
            'ref': 'tf_ref',
            'alt': 'tf_alt',
            'ID': 'SNP_id'
        }
        
        existing_cols_to_rename = {k: v for k, v in rename_map.items() if k in df_temp.columns}
        df_temp.rename(columns=existing_cols_to_rename, inplace=True)
        
        df_tf_list.append(df_temp)
        
    except Exception as e:
        print(f"Ошибка при чтении файла {file_path}: {e}")

if not df_tf_list:
    raise ValueError("Не удалось загрузить ни одного файла из папки TF.")

df_tf = pd.concat(df_tf_list, ignore_index=True)
print(f"Загружено данных TF: {len(df_tf)} строк из {len(tf_files)} файлов.")
print("Выполнение объединения таблиц...")

df_merged = pd.merge(df_mb, df_tf, on=['SNP_id', 'geneSymbol'], how='left')

cols_mb = [
    'seqnames', 'start', 'end', 'width', 'strand', 'SNP_id', 'REF', 'ALT', 
    'varType', 'motifPos', 'geneSymbol', 'dataSource', 'providerName', 
    'providerId', 'seqMatch', 'pctRef', 'pctAlt', 'scoreRef', 'scoreAlt', 
    'Refpvalue', 'Altpvalue', 'snpPos', 'alleleRef', 'alleleAlt', 'effect', 
    'altPos', 'alleleDiff', 'alleleEffectSize', 'motif_score_0.001', 
    'motif_score_0.0005', 'motif_score_0.0001', 'comparison_motif_score_0.001', 
    'comparison_motif_score_0.0005', 'comparison_motif_score_0.0001'
]

cols_tf = [
    'tf_chr', 'tf_start', 'tf_end', 'SNP_id', 'tf_ref', 'tf_alt', 'repeat_type', 
    'mean_BAD', 'mean_SNP_per_segment', 'n_aggregated', 'total_cover', 
    'es_mean_ref', 'es_mean_alt', 'fdrp_bh_ref', 'fdrp_bh_alt', 'motif_log_pref', 
    'motif_log_palt', 'motif_fc', 'motif_pos', 'motif_orient', 'motif_conc', 
    'motif_index'
]

# SNP_id и geneSymbol уже есть в cols_mb, убираем дубли из списка tf для финального порядка
cols_tf_final = [c for c in cols_tf if c not in ['SNP_id', 'geneSymbol']]

final_columns = cols_mb + cols_tf_final

# Проверка, все ли колонки присутствуют в объединенном датафрейме
missing_cols = [c for c in final_columns if c not in df_merged.columns]
if missing_cols:
    print(f"Предупреждение: Следующие столбцы не найдены в данных и будут пропущены: {missing_cols}")
    # Фильтруем список колонок только по тем, что реально есть
    final_columns = [c for c in final_columns if c in df_merged.columns]

df_final = df_merged[final_columns]

# --- 5. Сохранение результата ---
print(f"Сохранение результата в {OUTPUT_FILE}...")
df_final.to_csv(OUTPUT_FILE, sep='\t', index=False)

print("Готово!")
print(f"Итоговое количество строк: {len(df_final)}")
print(f"Итоговое количество столбцов: {len(df_final.columns)}")