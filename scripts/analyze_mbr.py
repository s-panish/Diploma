#!/usr/bin/env python3
"""
УСЛОВИЯ СРАВНЕНИЯ:
Для каждого порога motif_score_(###) (0.001, 0.0005, 0.0001):
- Если scoreRef > motif_score и scoreAlt < motif_score → "-2"
- Если scoreRef < motif_score и scoreAlt < motif_score → "-1"
- Если scoreRef > motif_score и scoreAlt > motif_score → "1"
- Если scoreRef < motif_score и scoreAlt > motif_score → "2"
"""
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from tqdm import tqdm

def process_row(row, threshold_name):
    score_ref = row['scoreRef']
    score_alt = row['scoreAlt']
    threshold = row[threshold_name]
    
    if score_ref > threshold and score_alt < threshold:
        return "-2"
    elif score_ref < threshold and score_alt < threshold:
        return "-1"
    elif score_ref > threshold and score_alt > threshold:
        return "1"
    elif score_ref < threshold and score_alt > threshold:
        return "2"
    else:
        return "0"

def process_chunk(chunk, threshold_names):
    for threshold in threshold_names:
        new_col_name = f"comparison_{threshold}"
        chunk[new_col_name] = chunk.apply(
            lambda row: process_row(row, threshold), axis=1
        )
    return chunk

def main():
    input_file = "motifbreakr_with_scores_001_log.tsv"
    output_file = "motifbreakr_analyzed_001_log.tsv"
    num_threads = 110

    print(f"Чтение данных из {input_file}...")
    df = pd.read_csv(input_file, sep='\t')

    threshold_columns = [
        'motif_score_0.001',
        'motif_score_0.0005',
        'motif_score_0.0001'
    ]

    chunk_size = max(1, len(df) // (num_threads * 4))
    chunks = [df[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
    
    print(f"Обработка {len(df)} строк с использованием {num_threads} потоков...")

    results = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(process_chunk, chunk.copy(), threshold_columns)
            for chunk in chunks
        ]
        
        for future in tqdm(as_completed(futures), total=len(futures)):
            results.append(future.result())

    result_df = pd.concat(results, ignore_index=True)

    print("Сортировка данных...")
    result_df = result_df.sort_values(by=['seqnames', 'start', 'end'])

    print(f"Сохранение результата в {output_file}...")
    result_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Анализ завершен! Результат сохранен в {output_file}")
    print(f"Добавлено {len(threshold_columns)} новых столбцов сравнения")

if __name__ == "__main__":
    main()