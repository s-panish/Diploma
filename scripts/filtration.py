import pandas as pd

def clean_tissue_name(value):
    if pd.isna(value):
        return value
    return str(value).rstrip('_')

df = pd.read_csv('tf_snv_tissue_pioneer_class_family.tsv', sep='\t')

# чистим названия тканей: убираем нижнее подчёркивание в конце
df['tissue'] = df['tissue'].apply(clean_tissue_name)

n_unique_before = df['tissue'].nunique()

filtered_df = df[(df['N(tf+t)'] >= 100) & (df['K(tf+t)'] >= 10)]

n_unique_after = filtered_df['tissue'].nunique()

print(f'Уникальных тканей в исходном файле: {n_unique_before}')
print(f'Уникальных тканей после фильтрации: {n_unique_after}')
print(f'Уникальных тканей отсеяно: {n_unique_before - n_unique_after}')

filtered_df.to_csv('tf_snv_tissue_pioneer_class_family_filtered.tsv', sep='\t', index=False)