import csv
import os
from collections import defaultdict

# Имена файлов
input_file = 'motifbr_001+adastra_001_significant_log.tsv'
output_file = 'motifbr+adastra_tf_top_001_log.tsv'

# Количество топ элементов (0 = все)
top_n = 0

def create_tf_top():
    if not os.path.exists(input_file):
        print(f"Ошибка: Файл '{input_file}' не найден.")
        return

    try:
        # Словарь: geneSymbol -> set уникальных SNP_id
        tf_snps = defaultdict(set)
        # Словарь: geneSymbol -> общее количество строк
        tf_counts = defaultdict(int)
        
        with open(input_file, 'r', encoding='utf-8', newline='') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            
            if reader.fieldnames is None:
                print("Ошибка: Не удалось прочитать заголовки файла.")
                return
            
            required_cols = ['geneSymbol', 'SNP_id']
            missing_cols = [col for col in required_cols if col not in reader.fieldnames]
            if missing_cols:
                print(f"Ошибка: Отсутствуют столбцы: {missing_cols}")
                return
            
            for row in reader:
                gene_symbol = row.get('geneSymbol', '').strip()
                snp_id = row.get('SNP_id', '').strip()
                
                if gene_symbol and snp_id:
                    tf_snps[gene_symbol].add(snp_id)  # Уникальные SNV
                    tf_counts[gene_symbol] += 1       # Всего строк
        
        # Сортируем по количеству уникальных SNV (по убыванию)
        sorted_tfs = sorted(tf_snps.items(), key=lambda x: len(x[1]), reverse=True)
        
        # Берем топ-N
        if top_n > 0:
            top_tfs = sorted_tfs[:top_n]
        else:
            top_tfs = sorted_tfs
        
        # Выводим результаты в консоль
        print(f"\n{'='*80}")
        print(f"ТОП {len(top_tfs)} TF по количеству уникальных SNV")
        print(f"{'='*80}")
        print(f"{'#':<5} {'geneSymbol':<15} {'Уникальных SNV':<18} {'Всего строк':<15}")
        print(f"{'-'*80}")
        
        for i, (gene_symbol, snps) in enumerate(top_tfs, 1):
            print(f"{i:<5} {gene_symbol:<15} {len(snps):<18} {tf_counts[gene_symbol]:<15}")
        
        print(f"{'='*80}")
        print(f"Всего уникальных TF: {len(tf_snps)}")
        print(f"Всего строк в файле: {sum(tf_counts.values())}")
        
        # Сохраняем результат в файл
        with open(output_file, 'w', encoding='utf-8', newline='') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            writer.writerow(['Rank', 'geneSymbol', 'Unique_SNV_Count'])
            for i, (gene_symbol, snps) in enumerate(sorted_tfs, 1):
                writer.writerow([i, gene_symbol, len(snps)])
        
        print(f"\nПолный результат сохранён в: {output_file}")

    except Exception as e:
        print(f"Произошла ошибка при обработке файла: {e}")

if __name__ == '__main__':
    create_tf_top()