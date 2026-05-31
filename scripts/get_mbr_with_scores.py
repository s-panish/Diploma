#!/usr/bin/env python3
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
import time

N_PROCESSES = 110

print("Загрузка таблицы пороговых значений...", file=sys.stderr)
scores_df = pd.read_csv('motif_scores_hocomoco_merged.tsv', sep='\t', dtype={'PWM': str})
scores_dict = scores_df.set_index('PWM')[['motif_score_0.001', 'motif_score_0.0005', 'motif_score_0.0001']].to_dict('index')

print("Загрузка motifbreakr данных...", file=sys.stderr)
col_count = len(open('motifbreakr_from_bed_001_log.tsv').readline().split('\t'))
df = pd.read_csv('motifbreakr_from_bed_001_log.tsv', sep='\t', header=None, names=range(col_count), engine='python', na_filter=False, skipinitialspace=True, quoting=3)
df = df.loc[:, (df != '').any(axis=0)].replace('', np.nan)
df = df.iloc[:, :28]
df.columns = ['seqnames', 'start', 'end', 'width', 'strand', 'SNP_id', 'REF', 'ALT', 'varType', 'motifPos', 'geneSymbol', 'dataSource', 'providerName', 'providerId', 'seqMatch', 'pctRef', 'pctAlt', 'scoreRef', 'scoreAlt', 'Refpvalue', 'Altpvalue', 'snpPos', 'alleleRef', 'alleleAlt', 'effect', 'altPos', 'alleleDiff', 'alleleEffectSize']

# Преобразуем числовые столбцы для корректной сортировки
df['start'] = pd.to_numeric(df['start'], errors='coerce')
df['end'] = pd.to_numeric(df['end'], errors='coerce')

total_rows = len(df)
chunk_size = max(1, total_rows // N_PROCESSES)
chunks = [(i, df.iloc[i*chunk_size:min((i+1)*chunk_size, total_rows)].copy()) for i in range(N_PROCESSES)]

def process_chunk(args):
    idx, chunk = args
    chunk['motif_score_0.001'] = chunk['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.001', np.nan))
    chunk['motif_score_0.0005'] = chunk['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.0005', np.nan))
    chunk['motif_score_0.0001'] = chunk['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.0001', np.nan))
    return idx, len(chunk)

print(f"\nОбработка {total_rows:,} строк ({N_PROCESSES} процессов)...", file=sys.stderr)
start_time = time.time()

with ProcessPoolExecutor(max_workers=N_PROCESSES) as executor:
    futures = [executor.submit(process_chunk, chunk) for chunk in chunks]
    completed = 0
    processed_rows = 0
    
    for future in as_completed(futures):
        idx, rows = future.result()
        completed += 1
        processed_rows += rows
        progress = completed / N_PROCESSES * 100
        elapsed = time.time() - start_time
        eta = elapsed / completed * (N_PROCESSES - completed) if completed > 0 else 0
        print(f"  {progress:5.1f}% | {completed:4d}/{N_PROCESSES} чанков | {processed_rows:10,}/{total_rows:,} строк | ETA: {eta:5.1f}s", 
              end='\r', file=sys.stderr)

print("\n\nДобавление пороговых значений...", file=sys.stderr)
df['motif_score_0.001'] = df['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.001', np.nan))
df['motif_score_0.0005'] = df['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.0005', np.nan))
df['motif_score_0.0001'] = df['providerId'].map(lambda x: scores_dict.get(x, {}).get('motif_score_0.0001', np.nan))

print("Сортировка по 'seqnames', 'start', 'end'...", file=sys.stderr)
df = df.sort_values(by=['seqnames', 'start', 'end'], ignore_index=True)

final_columns = ['seqnames', 'start', 'end', 'width', 'strand', 'SNP_id', 'REF', 'ALT', 'varType', 'motifPos', 'geneSymbol', 'dataSource', 'providerName', 'providerId', 'seqMatch', 'pctRef', 'pctAlt', 'scoreRef', 'scoreAlt', 'Refpvalue', 'Altpvalue', 'snpPos', 'alleleRef', 'alleleAlt', 'effect', 'altPos', 'alleleDiff', 'alleleEffectSize', 'motif_score_0.001', 'motif_score_0.0005', 'motif_score_0.0001']
df = df[final_columns]
df.to_csv('motifbreakr_with_scores_001_log.tsv', sep='\t', index=False, na_rep='NA', float_format='%.6f')

elapsed_total = time.time() - start_time
matched = df['motif_score_0.001'].notna().sum()
print(f"\n Готово за {elapsed_total:.1f} сек", file=sys.stderr)
print(f"   Обработано: {total_rows:,} строк", file=sys.stderr)
print(f"   Совпадений: {matched:,} ({matched/total_rows*100:.1f}%)", file=sys.stderr)
print(f"   Результат: motifbreakr_with_scores_001_log.tsv (отсортирован по seqnames, start, end)", file=sys.stderr)