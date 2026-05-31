import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

INPUT_FILE = "/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/merged_data_udacha.bed"
OUTPUT_FILE = "/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/merged_data_udacha_with_intervals.bed"
N_WORKERS = 12
INTERVAL_SIZE = 100

def process_chunk(args):
    """
    Обрабатывает кусок DataFrame: вычисляет новые интервалы.
    args: (chunk_df, chunk_id)
    start_interval 0-based
    end_interval 1-based
    """
    df, chunk_id = args
    
    # Вычисляем новые колонки
    # start_interval = start - 100 (не меньше 0)
    df['start_interval'] = df['start'].apply(lambda x: max(0, x - INTERVAL_SIZE))
    # end_interval = end + 100
    df['end_interval'] = df['end'] + INTERVAL_SIZE
    
    return df

def main():
    # Проверка существования файла
    if not os.path.exists(INPUT_FILE):
        print(f"Ошибка: Файл {INPUT_FILE} не найден.")
        return

    print(f"Чтение файла {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE, sep='\t')
    
    total_rows = len(df)
    print(f"Всего строк: {total_rows}")
    
    # Разбиваем DataFrame на части для параллельной обработки
    # np.array_split равномерно делит данные даже если они не делятся нацело
    chunks = np.array_split(df, N_WORKERS)
    
    print(f"Разбито на {len(chunks)} частей. Запуск обработки на {N_WORKERS} воркерах...")
    
    processed_chunks = []
    
    # Используем ProcessPoolExecutor для параллелизации
    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        tasks = [(chunk, i) for i, chunk in enumerate(chunks)]
        
        futures = {executor.submit(process_chunk, task): task[1] for task in tasks}
        
        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result_df = future.result()
                processed_chunks.append(result_df)
            except Exception as e:
                print(f"Ошибка в воркере {chunk_id}: {e}")
    
    if not processed_chunks:
        print("Ошибка: не удалось обработать ни одного чанка.")
        return

    print("Объединение результатов...")
    final_df = pd.concat(processed_chunks, ignore_index=True)
    
    print(f"Строк до удаления дубликатов: {len(final_df)}")
    
    # Фильтрация уникальных строк по первым трем колонкам
    final_df = final_df.drop_duplicates()
    
    print(f"Строк после удаления дубликатов: {len(final_df)}")
    
    expected_cols = ['chr', 'start', 'end', 'rsid', 'ref', 'alt', 
                     'ref_fdr_comb_pval', 'alt_fdr_comb_pval', 'cell_line',
                     'start_interval', 'end_interval']

    if all(col in final_df.columns for col in expected_cols):
        final_df = final_df[expected_cols]

    print(f"Сохранение в {OUTPUT_FILE}...")
    final_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print("Готово!")

if __name__ == "__main__":
    main()