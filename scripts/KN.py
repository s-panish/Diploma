#!/usr/bin/env python3

import pandas as pd
import numpy as np

input_file = "tf_snv_analysis_results/tf_snv_tissue_pioneer.tsv"
output_file = "tf_snv_analysis_results/tf_snv_tissue_pioneer_with_KN.tsv"

# Чтение файла
df = pd.read_csv(input_file, sep="\t", engine="python")

# Очистка имен колонок от лишних пробелов
df.columns = df.columns.astype(str).str.strip()

print("Исходные колонки:", df.columns.tolist())
print(f"Исходное количество строк: {len(df)}")

# Проверка наличия необходимых колонок
required_cols = ["TF", "N(tf+t)", "K(tf+t)", "tissue", "is_pioneer"]
missing_cols = [col for col in required_cols if col not in df.columns]

if missing_cols:
    raise ValueError(f"Отсутствуют необходимые колонки: {missing_cols}")

# Приведение к числовым типам (NaN остаются как есть)
df["N(tf+t)"] = pd.to_numeric(df["N(tf+t)"], errors="coerce")
df["K(tf+t)"] = pd.to_numeric(df["K(tf+t)"], errors="coerce")

# Очистка строковых колонок
df["TF"] = df["TF"].astype(str).str.strip()
df["tissue"] = df["tissue"].astype(str).str.strip()
df["is_pioneer"] = df["is_pioneer"].astype(str).str.strip().str.lower()

# Расчет K/N
# Избегаем деления на ноль, но сохраняем все строки
df["K/N"] = np.where(
    df["N(tf+t)"] > 0,
    df["K(tf+t)"] / df["N(tf+t)"],
    0  # Если N=0 или NaN, то K/N = 0
)

# Дополнительно: можно добавить проверку, что K не может быть больше N
invalid_kn = df[df["K(tf+t)"] > df["N(tf+t)"]]
if len(invalid_kn) > 0:
    print(f"\nПредупреждение: В {len(invalid_kn)} строках K > N:")
    print(invalid_kn[["TF", "N(tf+t)", "K(tf+t)", "tissue"]].head(10))

# НЕ удаляем строки с NaN - сохраняем все строки
print(f"\nВсего строк после расчета K/N: {len(df)}")
print(f"Из них строк с NaN в N(tf+t): {df['N(tf+t)'].isna().sum()}")
print(f"Из них строк с NaN в K(tf+t): {df['K(tf+t)'].isna().sum()}")

# Статистика по K/N (только для валидных значений)
valid_kn = df[df["K/N"].notna() & (df["N(tf+t)"] > 0)]
if len(valid_kn) > 0:
    print("\n=== Статистика по K/N (только для валидных значений) ===")
    print(f"Mean K/N: {valid_kn['K/N'].mean():.4f}")
    print(f"Median K/N: {valid_kn['K/N'].median():.4f}")
    print(f"Std K/N: {valid_kn['K/N'].std():.4f}")
    print(f"Min K/N: {valid_kn['K/N'].min():.4f}")
    print(f"Max K/N: {valid_kn['K/N'].max():.4f}")

# Статистика отдельно для пионерных и непионерных TF
print("\n=== Статистика по K/N для пионерных vs непионерных TF ===")
for pioneer_val in [True, False, "true", "false"]:
    mask = df["is_pioneer"] == str(pioneer_val).lower()
    valid_mask = mask & (df["K/N"].notna()) & (df["N(tf+t)"] > 0)
    if valid_mask.sum() > 0:
        group_name = "Pioneer" if str(pioneer_val).lower() == "true" else "Non-pioneer"
        print(f"\n{group_name} (n={valid_mask.sum()}):")
        print(f"  Mean K/N: {df.loc[valid_mask, 'K/N'].mean():.4f}")
        print(f"  Median K/N: {df.loc[valid_mask, 'K/N'].median():.4f}")
        print(f"  Std K/N: {df.loc[valid_mask, 'K/N'].std():.4f}")

# Переупорядочиваем колонки для удобства
column_order = ["TF", "tissue", "N(tf+t)", "K(tf+t)", "K/N", "is_pioneer"]
existing_cols = [col for col in column_order if col in df.columns]
other_cols = [col for col in df.columns if col not in column_order]
df = df[existing_cols + other_cols]

# Сохраняем результат (со всеми строками, включая пустые)
df.to_csv(output_file, sep="\t", index=False, na_rep="NaN")

print(f"\nГотово! Результат сохранен в файл: {output_file}")
print(f"Всего строк в выходном файле: {len(df)} (включая строки с NaN)")

# Показываем первые строки результата
print("\n=== Первые 10 строк результата ===")
print(df.head(10).to_string(index=False, na_rep="NaN"))

# Дополнительно: статистика по тканям (только для валидных значений)
print("\n=== Статистика по тканям (первые 10) ===")
tissue_stats = df[df["K/N"].notna() & (df["N(tf+t)"] > 0)].groupby("tissue").agg({
    "K/N": ["count", "mean", "median", "std"],
    "TF": "nunique"
}).round(4)
if len(tissue_stats) > 0:
    tissue_stats.columns = ["n_records", "mean_KN", "median_KN", "std_KN", "n_unique_TF"]
    print(tissue_stats.head(10).to_string())
    
    # Сохраняем также агрегированную статистику по тканям
    tissue_stats_file = "tf_snv_analysis_results/tissue_KN_statistics.tsv"
    tissue_stats.to_csv(tissue_stats_file, sep="\t")
    print(f"\nСтатистика по тканям сохранена в: {tissue_stats_file}")
else:
    print("Нет валидных данных для статистики по тканям")