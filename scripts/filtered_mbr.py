#!/usr/bin/env python3
"""
УСЛОВИЯ ФИЛЬТРАЦИИ:
Отбираются строки, в которых хотя бы в одном из столбцов:
- comparison_motif_score_0.001
- comparison_motif_score_0.0005  
- comparison_motif_score_0.0001
значение равно "-2" или "2"
"""

import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from tqdm import tqdm
import os

def filter_chunk(chunk, comparison_columns):
    mask = pd.Series(False, index=chunk.index)
    
    for col in comparison_columns:
        if col in chunk.columns:
            col_mask = (chunk[col] == -2) | (chunk[col] == 2)
            mask = mask | col_mask
    
    filtered_chunk = chunk[mask].copy()
    
    return filtered_chunk

def process_file(input_file, output_file, num_threads=110):
    print(f"Чтение данных из {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"Всего строк: {len(df)}")

    comparison_columns = [
        'comparison_motif_score_0.001',
        'comparison_motif_score_0.0005',
        'comparison_motif_score_0.0001'
    ]
    
    available_columns = [col for col in comparison_columns if col in df.columns]
    missing_columns = [col for col in comparison_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Внимание: отсутствуют столбцы: {missing_columns}")
    
    if not available_columns:
        print("Ошибка: ни один из столбцов сравнения не найден!")
        return
    
    print(f"Доступные столбцы для фильтрации: {available_columns}")

    chunk_size = max(1, len(df) // (num_threads * 4))
    chunks = [df[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
    
    print(f"Разделение на {len(chunks)} чанков для обработки...")

    filtered_chunks = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        futures = [
            executor.submit(filter_chunk, chunk, available_columns)
            for chunk in chunks
        ]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Фильтрация"):
            result = future.result()
            if len(result) > 0:
                filtered_chunks.append(result)

    if filtered_chunks:
        result_df = pd.concat(filtered_chunks, ignore_index=True)
    else:
        result_df = pd.DataFrame(columns=df.columns)
    
    print(f"После фильтрации: {len(result_df)} строк")

    print("Сортировка данных...")
    result_df_sorted = result_df.sort_values(
        by=['seqnames', 'start', 'end']
    )

    print(f"Сохранение результата в {output_file}...")
    result_df_sorted.to_csv(output_file, sep='\t', index=False)

    print("\n=== СТАТИСТИКА ===")
    for col in available_columns:
        if col in result_df_sorted.columns:
            value_counts = result_df_sorted[col].value_counts()
            print(f"\nСтолбец {col}:")
            for value, count in value_counts.items():
                print(f"  {value}: {count} строк")
    
    print(f"\nАнализ завершен! Результат сохранен в {output_file}")
    print(f"Итоговое количество строк: {len(result_df_sorted)}")

def main():
    input_file = "motifbreakr_analyzed_001_log.tsv"
    output_file = "motifbreakr_filtered_001_log.tsv"
    num_threads = 110

    process_file(input_file, output_file, num_threads)

if __name__ == "__main__":
    main()