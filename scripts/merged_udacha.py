import os
import glob
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  # Для отображения прогресса (опционально, если установлена)

# ================= КОНФИГУРАЦИЯ =================
# Пути к директориям
DIR_ATAC = "/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/atac"
DIR_DNASE = "/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/dnase"

# Выходной файл
OUTPUT_FILE = "/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/merged_data_udacha.bed"

# Количество воркеров (потоков/процессов)
N_WORKERS = 12

# Названия колонок согласно описанию
COLUMN_NAMES = [
    'chr', 
    'start', 
    'end', 
    'rsid', 
    'ref', 
    'alt', 
    'ref_fdr_comb_pval', 
    'alt_fdr_comb_pval'
]
# =================================================

def get_bed_files(directories):
    """Собирает список всех .bed файлов из указанных директорий."""
    files = []
    for directory in directories:
        if not os.path.exists(directory):
            print(f"Предупреждение: Директория не найдена: {directory}")
            continue
        # Ищем все файлы с расширением .bed
        pattern = os.path.join(directory, "*.bed")
        files.extend(glob.glob(pattern))
    return files

def process_file(filepath):
    """
    Читает один bed файл, добавляет колонку cell_line и возвращает DataFrame.
    Эта функция будет выполняться в отдельном процессе.
    """
    try:
        # Извлекаем имя клеточной линии из имени файла
        filename = os.path.basename(filepath)
        cell_line = filename.replace('.bed', '')
        
        # Читаем файл (отсутствие заголовка и разделитель табуляция)
        df = pd.read_csv(
            filepath, 
            sep='\t', 
            header=None, 
            names=COLUMN_NAMES,
            dtype={'chr': str, 'start': int, 'end': int} # Оптимизация типов
        )
        
        df['cell_line'] = cell_line
        
        return df
    except Exception as e:
        print(f"Ошибка при обработке файла {filepath}: {e}")
        return None

def main():
    print("Поиск файлов...")
    bed_files = get_bed_files([DIR_ATAC, DIR_DNASE])
    
    if not bed_files:
        print("Файлы .bed не найдены в указанных директориях.")
        return

    print(f"Найдено файлов: {len(bed_files)}. Запуск обработки на {N_WORKERS} воркерах...")
    
    all_dfs = []
    
    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        future_to_file = {executor.submit(process_file, f): f for f in bed_files}

        iterator = as_completed(future_to_file)
        try:
            from tqdm import tqdm
            iterator = tqdm(iterator, total=len(bed_files), desc="Обработка файлов")
        except ImportError:
            pass # Если tqdm не установлен, просто итерируемся

        for future in iterator:
            df = future.result()
            if df is not None:
                all_dfs.append(df)
    
    if not all_dfs:
        print("Не удалось прочитать ни одного файла.")
        return

    print("Объединение таблиц...")
    merged_df = pd.concat(all_dfs, ignore_index=True)
    
    print(f"Всего строк до фильтрации: {len(merged_df)}")
    
    # Фильтрация уникальных строк по первым трем колонкам (chr, start, end)
    # keep='first' оставляет первую встретившуюся строку, остальные дубли удаляет
    print("Удаление дубликатов ")
    unique_df = merged_df.drop_duplicates()
    
    print(f"Всего строк после фильтрации: {len(unique_df)}")
    print(f"Сохранение результата в {OUTPUT_FILE}...")
    unique_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
    
    print("Готово!")

if __name__ == "__main__":
    main()