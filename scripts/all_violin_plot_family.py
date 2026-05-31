#!/usr/bin/env python3

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import PolyCollection

input_file = 'tf_snv_analysis_results/total_tf_snv_pioneer_class_family.tsv'
output_dir = 'violin_plots_by_family'
os.makedirs(output_dir, exist_ok=True)

plt.rcParams['ps.fonttype'] = 42

df = pd.read_csv(input_file, sep='\t', engine='python')

# если есть мусорный индекс-столбец, удаляем
first_col = str(df.columns[0]).strip()
if first_col.startswith('Unnamed') or first_col == '':
    df = df.drop(columns=df.columns[0])

# чистим имена столбцов
df.columns = [str(c).strip() for c in df.columns]

required_cols = ['TF', 'is_pioneer', 'K/N', 'TF_family']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(
        f'В файле отсутствуют обязательные столбцы: {missing}\n'
        f'Найденные столбцы: {list(df.columns)}'
    )

df['K/N'] = pd.to_numeric(df['K/N'], errors='coerce')
df['TF_family'] = df['TF_family'].astype(str).str.strip()
df['TF'] = df['TF'].astype(str).str.strip()

def pioneer_label(x):
    x = str(x).strip().lower()
    if x == 'true':
        return 'ПТФ'
    return 'ТФ'

df['TF_type'] = df['is_pioneer'].apply(pioneer_label)

df = df.dropna(subset=['K/N', 'TF_family', 'TF_type', 'TF']).copy()
df = df[df['TF_family'] != ''].copy()
df = df[df['TF'] != ''].copy()

type_order = ['ПТФ', 'ТФ']

palette = {
    'ПТФ': '#D7AEF6',
    'ТФ': '#70E0E3'
}

strip_palette = {
    'ПТФ': 'gray',
    'ТФ': 'gray'
}

sns.set_style('white')

families = sorted(df['TF_family'].unique())
print(f'Всего семейств: {len(families)}')

valid_families = []
skipped_count = 0

for family in families:
    sub = df[df['TF_family'] == family].copy()

    if sub.empty:
        skipped_count += 1
        continue

    counts = sub.groupby('TF_type')['TF'].nunique()
    n_ptf = counts.get('ПТФ', 0)
    n_tf = counts.get('ТФ', 0)

    if n_ptf < 3 or n_tf < 3:
        print(
            f"Пропуск семейства '{family}': "
            f'ПТФ = {n_ptf}, ТФ = {n_tf} (нужно минимум 3 и 3)'
        )
        skipped_count += 1
        continue

    valid_families.append(family)

print(f'Подходящих семейств: {len(valid_families)}')
print(f'Пропущено семейств: {skipped_count}')

if len(valid_families) == 0:
    raise ValueError('Нет ни одного семейства, где есть минимум 3 ПТФ и минимум 3 ТФ.')

plot_df = df[df['TF_family'].isin(valid_families)].copy()
family_order = valid_families

# таблица с количеством уникальных TF для каждой пары семейство × группа
count_table = (
    plot_df
    .groupby(['TF_family', 'TF_type'])['TF']
    .nunique()
    .unstack(fill_value=0)
)

fig_width = max(18, len(family_order) * 2.8)
fig_height = 10.5

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# violin
sns.violinplot(
    data=plot_df,
    x='TF_family',
    y='K/N',
    hue='TF_type',
    order=family_order,
    hue_order=type_order,
    palette=palette,
    cut=0,
    inner='box',
    linewidth=1.3,
    dodge=True,
    ax=ax
)

# points
sns.stripplot(
    data=plot_df,
    x='TF_family',
    y='K/N',
    hue='TF_type',
    order=family_order,
    hue_order=type_order,
    dodge=True,
    palette=strip_palette,
    alpha=0.55,
    size=3.5,
    jitter=0.08,
    ax=ax
)

# убираем дубли легенды
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles[:2],
    labels[:2],
    title='',
    frameon=False,
    fontsize=13,
    loc='upper right'
)

# оформление
ax.set_title('The fraction of sites with SNV by TF family', fontsize=18, pad=18, color='black')
ax.set_xlabel('')
ax.set_ylabel('The fraction of sites with SNV', fontsize=15, color='black')

plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12, color='black')

ax.spines['bottom'].set_color('gray')
ax.spines['left'].set_color('gray')
ax.spines['top'].set_color('gray')
ax.spines['right'].set_color('gray')

ax.tick_params(axis='x', colors='gray', labelcolor='black', labelsize=12, pad=10)
ax.tick_params(axis='y', colors='gray', labelcolor='black', labelsize=13)

ax.xaxis.label.set_color('black')
ax.yaxis.label.set_color('black')

ax.grid(axis='y', color='lightgray', linestyle='--', linewidth=0.7, alpha=0.7)

# диапазон по Y + место под подписи количества TF
y_min = plot_df['K/N'].min()
y_max = plot_df['K/N'].max()
y_range = y_max - y_min if y_max > y_min else 1.0

ax.set_ylim(y_min, y_max + y_range * 0.16)

# чуть-чуть поля справа, чтобы ничего не прилипало
ax.set_xlim(-0.5, len(family_order) - 0.5 + 0.08)

# -------- точное определение центров violin --------
# берем только сами тела violin из PolyCollection
violin_bodies = [c for c in ax.collections if isinstance(c, PolyCollection)]

# для каждого violin вычисляем средний x по его вершинам
violin_centers = []
for body in violin_bodies:
    paths = body.get_paths()
    if not paths:
        continue
    verts = paths[0].vertices
    x_center = verts[:, 0].mean()
    violin_centers.append(x_center)

# сортируем слева направо
violin_centers = sorted(violin_centers)

expected_n = len(family_order) * 2
if len(violin_centers) < expected_n:
    raise ValueError(
        f'Найдено только {len(violin_centers)} violin bodies, ожидалось минимум {expected_n}.'
    )

# если seaborn создал лишние объекты, берем первые нужные после сортировки
violin_centers = violin_centers[:expected_n]

label_y = y_max + y_range * 0.025

# на каждое семейство приходится 2 violin: ПТФ и ТФ
for i, family in enumerate(family_order):
    x_ptf = violin_centers[2 * i]
    x_tf = violin_centers[2 * i + 1]

    n_ptf = int(count_table.loc[family, 'ПТФ']) if 'ПТФ' in count_table.columns else 0
    n_tf = int(count_table.loc[family, 'ТФ']) if 'ТФ' in count_table.columns else 0

    ax.text(
        x_ptf,
        label_y,
        f'n={n_ptf}',
        ha='center',
        va='bottom',
        fontsize=11,
        color='black'
    )

    ax.text(
        x_tf,
        label_y,
        f'n={n_tf}',
        ha='center',
        va='bottom',
        fontsize=11,
        color='black'
    )

plt.tight_layout()

png_path = os.path.join(output_dir, 'all_families_one_plot_with_tf_counts.png')
pdf_path = os.path.join(output_dir, 'all_families_one_plot_with_tf_counts.pdf')

plt.savefig(png_path, dpi=300, bbox_inches='tight')
plt.savefig(pdf_path, bbox_inches='tight')
plt.close(fig)

print('Готово.')
print(f'Сохранено семейств на общем графике: {len(valid_families)}')
print(f'PNG: {png_path}')
print(f'PDF: {pdf_path}')