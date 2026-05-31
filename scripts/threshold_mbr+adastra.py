import csv
import os

input_file = 'motifbr+adastra_clean_001_log.tsv'
output_file = 'motifbr_001+adastra_001_significant_log.tsv'

# Пороговое значение
threshold = 0.01

def filter_significant():
    try:
        with open(input_file, 'r', encoding='utf-8', newline='') as f_in, \
             open(output_file, 'w', encoding='utf-8', newline='') as f_out:
            
            reader = csv.DictReader(f_in, delimiter='\t')
            
            if reader.fieldnames is None:
                print("Ошибка: Не удалось прочитать заголовки файла.")
                return

            # Проверка наличия нужных столбцов
            required_cols = ['fdrp_bh_ref', 'fdrp_bh_alt']
            missing_cols = [col for col in required_cols if col not in reader.fieldnames]
            if missing_cols:
                print(f"Ошибка: Отсутствуют столбцы: {missing_cols}")
                return

            writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
            writer.writeheader()

            count_total = 0
            count_kept = 0
            count_ref = 0
            count_alt = 0
            count_both = 0

            for row in reader:
                count_total += 1
                
                try:
                    # Получаем значения и преобразуем к float
                    fdrp_ref = row.get('fdrp_bh_ref', '').strip()
                    fdrp_alt = row.get('fdrp_bh_alt', '').strip()
                    
                    # Пропускаем строки с пустыми значениями
                    if not fdrp_ref or not fdrp_alt:
                        continue
                    
                    fdrp_ref_val = float(fdrp_ref)
                    fdrp_alt_val = float(fdrp_alt)
                    
                    # Проверяем условие
                    ref_sig = fdrp_ref_val < threshold
                    alt_sig = fdrp_alt_val < threshold
                    
                    if ref_sig or alt_sig:
                        writer.writerow(row)
                        count_kept += 1
                        
                        if ref_sig and alt_sig:
                            count_both += 1
                        elif ref_sig:
                            count_ref += 1
                        elif alt_sig:
                            count_alt += 1
                            
                except ValueError as e:
                    # Пропускаем строки с некорректными числовыми значениями
                    continue

            print(f"Всего строк обработано: {count_total}")
            print(f"Строк прошло фильтр: {count_kept}")
            print(f"  - только fdrp_bh_ref < {threshold}: {count_ref}")
            print(f"  - только fdrp_bh_alt < {threshold}: {count_alt}")
            print(f"  - оба значения < {threshold}: {count_both}")
            print(f"Результат записан в: {output_file}")

    except Exception as e:
        print(f"Произошла ошибка при обработке файла: {e}")

if __name__ == '__main__':
    filter_significant()