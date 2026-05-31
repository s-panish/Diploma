import csv
import os

input_file = 'motifbr+adastra_001_log.tsv'
output_file = 'motifbr+adastra_clean_001_log.tsv'

target_columns = ['es_mean_ref', 'es_mean_alt', 'fdrp_bh_ref', 'fdrp_bh_alt']

def filter_tsv():
    if not os.path.exists(input_file):
        print(f"Ошибка: Файл '{input_file}' не найден.")
        return

    try:
        with open(input_file, 'r', encoding='utf-8', newline='') as f_in, \
             open(output_file, 'w', encoding='utf-8', newline='') as f_out:
            
            reader = csv.DictReader(f_in, delimiter='\t')

            if reader.fieldnames is None:
                print("Ошибка: Не удалось прочитать заголовки файла.")
                return

            missing_cols = [col for col in target_columns if col not in reader.fieldnames]
            if missing_cols:
                print(f"Ошибка: В файле отсутствуют следующие столбцы: {missing_cols}")
                return

            writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
            writer.writeheader()

            count_total = 0
            count_kept = 0

            for row in reader:
                count_total += 1
                keep_row = True

                for col in target_columns:
                    val = row.get(col, '')
                    if val is None or not str(val).strip():
                        keep_row = False
                        break
                
                if keep_row:
                    writer.writerow(row)
                    count_kept += 1

            print(f"Обработано строк: {count_total}")
            print(f"Сохранено строк: {count_kept}")
            print(f"Результат записан в: {output_file}")

    except Exception as e:
        print(f"Произошла ошибка при обработке файла: {e}")

if __name__ == '__main__':
    filter_tsv()