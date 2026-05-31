import pandas as pd
import os
import numpy as np

file_path = '/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/tf_snv_analysis_results/tf_snv_summary_by_tissue.tsv'

print("Чтение файла...")
df = pd.read_csv(file_path, sep='\t', low_memory=False)

print("Преобразование столбцов к числовому типу...")
df['N(tf+t)'] = pd.to_numeric(df['N(tf+t)'], errors='coerce')
df['K(tf+t)'] = pd.to_numeric(df['K(tf+t)'], errors='coerce')

print("Вычисление столбца K/N...")
df['K/N'] = df['K(tf+t)'] / df['N(tf+t)']
df['K/N'] = df['K/N'].replace([np.inf, -np.inf], np.nan)
df.loc[df['K(tf+t)'] == 0, 'K/N'] = np.nan
df['K/N'] = df['K/N'].round(4)

print(f"Строк до удаления: {len(df)}")
# df = df.dropna(subset=['K/N'])  # Удаляем только где K/N пустой
print(f"Строк после удаления: {len(df)}")
print("Сохранение файла...")
df.to_csv(file_path, sep='\t', index=False)
print("Готово!")