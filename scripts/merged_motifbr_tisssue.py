import pandas as pd
import numpy as np
import multiprocessing as mp
from functools import partial
import os

def process_chunk(chunk_df, bed_df):
    merged = pd.merge(chunk_df, bed_df, left_on='SNP_id', right_on='rsid', how='left', sort=False)
    
    if 'rsid' in merged.columns:
        merged = merged.drop(columns=['rsid'])    
    return merged.drop_duplicates()

def main():
    BED_FILE = 'merged_data_udacha_cleaned.bed'
    MOTIF_FILE = 'motifbr_001+adastra_001_significant_analyzed_log.tsv'
    OUTPUT_FILE = 'motifbreakr_adastra_with_cell_lines.tsv'
    N_WORKERS = 40  

    print(" Загрузка справочной таблицы (BED)...")
    bed_df = pd.read_csv(BED_FILE, sep='\t', usecols=['rsid', 'cell_line'])
    bed_df = bed_df.drop_duplicates()  

    print(" Загрузка основной таблицы...")
    motif_df = pd.read_csv(MOTIF_FILE, sep='\t')
    print(f"  Строк в motifbreakr: {len(motif_df):,}")

    print(f" Разбиение на {N_WORKERS} частей...")
    chunks = np.array_split(motif_df, N_WORKERS)
    del motif_df 

    print(f"⚡ Запуск {N_WORKERS} параллельных процессов...")

    process_func = partial(process_chunk, bed_df=bed_df)

    with mp.Pool(processes=N_WORKERS) as pool:
        results = pool.map(process_func, chunks)

    print(" Сборка результатов...")
    final_df = pd.concat(results, ignore_index=True)
    del results 

    print(" Финальное удаление полных дубликатов (глобально)...")
    final_df = final_df.drop_duplicates()

    print(f"Сохранение в {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print(f" Готово! Итоговых уникальных строк: {len(final_df):,}")

if __name__ == '__main__':
    if os.name != 'nt':
        mp.set_start_method('fork', force=True)
    main()