import pandas as pd
import numpy as np
from scipy.stats import hypergeom

# ================= НАСТРОЙКИ =================
INPUT_FILE = 'motifbreakr_with_cell_lines_unique_broken_expression_ge1.tsv'
OUTPUT_FILE = 'top5_frequency_enrichment_motifbr_expression.tsv'
TOP5_PAIRS_FILE = 'top5_tf_cell_lines_expression.tsv'
# =============================================

print("1. Загрузка данных...")
df = pd.read_csv(INPUT_FILE, sep='\t')
print(f"   Всего строк: {len(df)} | Уникальных TF: {df['providerId'].nunique()}")

# 2. Определение "сломанных" мотивов
score_cols = ['comparison_motif_score_0.001', 'comparison_motif_score_0.0005', 'comparison_motif_score_0.0001']
df[score_cols] = df[score_cols].apply(pd.to_numeric, errors='coerce')
df['is_broken'] = df[score_cols].abs().eq(2).any(axis=1)

# Извлекаем базовое имя ТФ из providerId (например, AHRR.H13CORE... -> AHRR)
df['TF'] = df['providerId'].str.split('.').str[0]

# 3. Агрегация по каждой паре TF × cell_line
print("2. Расчёт k, K и соотношения для каждой пары TF×ткань...")
tf_cell_agg = df.groupby(['TF', 'cell_line']).agg(
    k=('is_broken', 'sum'),
    K=('is_broken', 'size')
).reset_index()

# ==================== НОВЫЙ ФИЛЬТР ====================
# Подсчёт уникальных тканей ДО фильтрации
unique_tissues_before = tf_cell_agg['cell_line'].nunique()
print(f"   Уникальных тканей ДО фильтрации: {unique_tissues_before}")

print("   Применяем фильтр: K >= 100 и k >= 10...")
tf_cell_agg = tf_cell_agg[(tf_cell_agg['K'] >= 100) & (tf_cell_agg['k'] >= 10)].copy()

# Подсчёт уникальных тканей ПОСЛЕ фильтрации
unique_tissues_after = tf_cell_agg['cell_line'].nunique()
print(f"   Уникальных тканей ПОСЛЕ фильтрации: {unique_tissues_after}")
print(f"   Осталось пар TF×ткань после фильтрации: {len(tf_cell_agg)}")
# =====================================================

if len(tf_cell_agg) == 0:
    raise ValueError("После фильтрации не осталось данных для анализа. Проверьте входные данные или пороги K/k.")

tf_cell_agg['ratio'] = tf_cell_agg['k'] / tf_cell_agg['K']

# 4. Определение глобального порога топ-5%
print("3. Расчёт глобального топ-5% по соотношению k/K...")
threshold_95 = np.percentile(tf_cell_agg['ratio'], 95)
tf_cell_agg['in_top5'] = tf_cell_agg['ratio'] >= threshold_95

n = tf_cell_agg['in_top5'].sum()  # n: всего пар в топ-5% (после фильтра)
N = len(tf_cell_agg)              # N: всего пар в датасете (после фильтра)

print(f"   Порог 95-го перцентиля: {threshold_95:.4f}")
print(f"   n (в топ-5%): {n}")
print(f"   N (всего после фильтра): {N}")

# === СОХРАНЕНИЕ ВСЕХ ПАР TF × ТКАНЬ, ПОПАВШИХ В ТОП-5% ===
print("4. Сохранение всех пар TF×ткань, попавших в топ-5% (ratio >= 95-й перцентиль)...")
top5_percent_pairs = tf_cell_agg[tf_cell_agg['in_top5']].copy()
top5_percent_pairs = top5_percent_pairs.sort_values('ratio', ascending=False).reset_index(drop=True)

top5_cols = ['TF', 'cell_line', 'k', 'K', 'ratio', 'in_top5']
top5_percent_pairs[top5_cols].to_csv(TOP5_PAIRS_FILE, sep='\t', index=False)
print(f"Сохранено {len(top5_percent_pairs)} пар в {TOP5_PAIRS_FILE}")
print(top5_percent_pairs[['TF', 'cell_line', 'ratio']].head(10).to_string(index=False))
# ==========================================

# 5. Агрегация по каждому ТФ: сколько его тканей попало в топ-5%
print("5. Подсчёт частоты вхождения каждого ТФ в топ-5%...")
tf_agg = tf_cell_agg.groupby('TF')['in_top5'].agg(
    k=('sum'),          # k: число тканей, где ТФ в топ-5%
    K=('size')          # K: число тканей, где ТФ НЕ в топ-5% (в рамках фильтра)
).reset_index()

tf_agg['n'] = n
tf_agg['N'] = N

# 6. Гипергеометрический тест: P(X >= k)
print("6. Расчёт p-value...")
tf_agg['p_value'] = tf_agg.apply(
    lambda row: hypergeom.sf(row['k'] - 1, row['N'], row['n'], row['k'] + row['K']), axis=1
)

# 7. Fold Enrichment
obs_freq = tf_agg['k'] / (tf_agg['k'] + tf_agg['K'])
bg_freq = n / N
tf_agg['fold_enrichment'] = np.where(
    (tf_agg['k'] + tf_agg['K']) > 0, 
    obs_freq / bg_freq, 
    0.0
)

# 8. FDR коррекция (Benjamini-Hochberg)
def fdr_bh(pvals):
    pvals = np.asarray(pvals, dtype=float)
    m = len(pvals)
    if m == 0: return pvals
    order = pvals.argsort()
    sorted_pvals = pvals[order]
    qvals = np.zeros(m)
    qvals[-1] = sorted_pvals[-1]
    for i in range(m - 2, -1, -1):
        qvals[i] = min(sorted_pvals[i] * m / (i + 1), qvals[i + 1])
    return qvals[order.argsort()]

tf_agg['p_value_fdr'] = fdr_bh(tf_agg['p_value'].values)

# 9. Сортировка и сохранение основной таблицы
print("7. Сортировка и сохранение основной таблицы...")
tf_agg = tf_agg.sort_values('p_value_fdr', ascending=True).reset_index(drop=True)

out_cols = ['TF', 'k', 'K', 'n', 'N', 'p_value', 'p_value_fdr', 'fold_enrichment']
tf_agg[out_cols].to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f"\n Готово! Результат сохранён в {OUTPUT_FILE}")
print(f"  ТФ в таблице: {len(tf_agg)}")
print(f"  Значимых при FDR < 0.05: {(tf_agg['p_value_fdr'] < 0.05).sum()}")
print(f" Медианный fold_enrichment: {tf_agg['fold_enrichment'].median():.3f}")