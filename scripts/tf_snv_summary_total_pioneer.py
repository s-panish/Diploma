#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

# Пути к файлам
pioneer_tf_path = '/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/PioneerTF_list_apr2026.xlsx'
snv_summary_path = '/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/tf_snv_analysis_results/tf_snv_summary.tsv'
output_path = '/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/tf_snv_analysis_results/tf_snv_tissue_pioneer.tsv'

# 1. Читаем таблицу с pioneer TF
pioneer_df = pd.read_excel(pioneer_tf_path, sheet_name='Table_S10', skiprows=1)
pioneer_tf_list = set(pioneer_df['TF'].dropna().astype(str).str.strip())

print(f"Найдено {len(pioneer_tf_list)} pioneer TF в таблице Pioneer_TF_list.xlsx")
print(f"Примеры pioneer TF: {list(pioneer_tf_list)[:10]}")

# 2. Читаем таблицу с SNV summary
snv_df = pd.read_csv(snv_summary_path, sep='\t')
print(f"\nЗагружено {len(snv_df)} записей из tf_snv_summary_total.tsv")

# 3. Добавляем столбец is_pioneer
snv_df['is_pioneer'] = snv_df['TF'].astype(str).str.strip().isin(pioneer_tf_list)

# 4. Сохраняем результат
snv_df.to_csv(output_path, sep='\t', index=False)
print(f"\nРезультат сохранён в: {output_path}")

# 5. Выводим статистику
print(f"\nСтатистика:")
print(f"  - Всего TF: {len(snv_df)}")
print(f"  - Pioneer TF: {snv_df['is_pioneer'].sum()}")
print(f"  - Не pioneer TF: {(~snv_df['is_pioneer']).sum()}")

# 6. Показываем pioneer TF из вашего файла
pioneer_found = snv_df[snv_df['is_pioneer']]['TF'].tolist()
if pioneer_found:
    print(f"\nPioneer TF найдены в tf_snv_summary_total.tsv: {pioneer_found}")
else:
    print("\nPioneer TF не найдены в tf_snv_summary_total.tsv")