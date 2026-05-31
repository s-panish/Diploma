#!/bin/bash

# Настройка путей
SNV_FILE="/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/merged_data_udacha_with_intervals.bed"
TF_DIR="/shared_data/rscf_scatac-seq_2024-2025/hocomoco_tracks"
OUT_DIR="./tf_snv_analysis_results"
TEMP_DIR=$(mktemp -d)

mkdir -p "$OUT_DIR"
mkdir -p "$TEMP_DIR"

echo "=== Подготовка данных SNV ==="

# 1. Точные позиции SNV (chr, start, end, tissue)
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $2, $3, $9}' "$SNV_FILE" | sort -k1,1 -k2,2n > "$TEMP_DIR/snv_exact.bed"

# 2. Интервалы вокруг SNV (chr, start_interval, end_interval, tissue)
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 {print $1, $10, $11, $9}' "$SNV_FILE" | sort -k1,1 -k2,2n > "$TEMP_DIR/snv_interval.bed"

echo "Строк в snv_exact.bed: $(wc -l < "$TEMP_DIR/snv_exact.bed")"
echo "Строк в snv_interval.bed: $(wc -l < "$TEMP_DIR/snv_interval.bed")"

echo "=== Запуск параллельной обработки TF (12 потоков) ==="

process_tf() {
    local tf_file="$1"
    local temp_dir="$2"
    local snv_exact="$3"
    local snv_interval="$4"
    
    local base_name=$(basename "$tf_file" .bed)
    local tf_norm="$temp_dir/${base_name}_norm.bed"
    local out_k="$temp_dir/${base_name}_K.tmp"
    local out_n="$temp_dir/${base_name}_N.tmp"
    
    # 1. Нормализация хромосом
    grep -v "^track" "$tf_file" 2>/dev/null | grep -v "^browser" | grep -v "^#" | \
    awk 'BEGIN{OFS="\t"} {
        if ($1 !~ /^chr/) {
            print "chr"$1, $2, $3, $4, $5, $6
        } else {
            print $1, $2, $3, $4, $5, $6
        }
    }' | sort -k1,1 -k2,2n > "$tf_norm"
    
    if [ ! -s "$tf_norm" ]; then
        return 0
    fi
    
    # 2. Извлечение имени ТФ
    local tf_name=$(head -n1 "$tf_norm" | cut -f5 | awk -F'_HUMAN' '{print $1}')
    [ -z "$tf_name" ] && tf_name="$base_name"
    
    # 3. Расчет K (точное попадание SNV в мотив)
    ### ИЗМЕНЕНО ###
    # Выводим: tissue \t motif_key (без подсчета, подсчет будет в финале)
    bedtools intersect -a "$tf_norm" -b "$snv_exact" -wa -wb 2>/dev/null | \
    awk 'BEGIN{OFS="\t"} {
        motif_key = $1":"$2":"$3;
        tissue = $10;
        print tissue, motif_key;
    }' > "$out_k"
    
    # 4. Расчет N (попадание мотива в интервал SNV)
    ### ИЗМЕНЕНО ###
    bedtools intersect -a "$tf_norm" -b "$snv_interval" -wa -wb 2>/dev/null | \
    awk 'BEGIN{OFS="\t"} {
        motif_key = $1":"$2":"$3;
        tissue = $10;
        print tissue, motif_key;
    }' > "$out_n"
    
    rm -f "$tf_norm"
}

export -f process_tf
export TEMP_DIR

find "$TF_DIR" -name "*.bed" -type f | \
parallel -j 12 --halt now,fail=1 process_tf {} "$TEMP_DIR" "$TEMP_DIR/snv_exact.bed" "$TEMP_DIR/snv_interval.bed"

echo "=== Агрегация результатов ==="

K_FILES=$(find "$TEMP_DIR" -name "*_K.tmp" -type f | wc -l)
N_FILES=$(find "$TEMP_DIR" -name "*_N.tmp" -type f | wc -l)
echo "Файлов с K: $K_FILES, с N: $N_FILES"

# Собираем все пересечения в один файл
cat "$TEMP_DIR"/*_K.tmp 2>/dev/null > "$TEMP_DIR/all_K.tmp"
cat "$TEMP_DIR"/*_N.tmp 2>/dev/null > "$TEMP_DIR/all_N.tmp"

# Подсчет УНИКАЛЬНЫХ мотивов по тканям (глобально по всем TF)
# Сортируем и убираем дубликаты (один мотив в одной ткани считается 1 раз)
echo "Подсчет уникальных мотивов для K..."
sort -u "$TEMP_DIR/all_K.tmp" | \
awk -F'\t' '{count[$1]++} END {for (t in count) print t, count[t], "K"}' OFS='\t' > "$TEMP_DIR/K_by_tissue.tmp"

echo "Подсчет уникальных мотивов для N..."
sort -u "$TEMP_DIR/all_N.tmp" | \
awk -F'\t' '{count[$1]++} END {for (t in count) print t, count[t], "N"}' OFS='\t' > "$TEMP_DIR/N_by_tissue.tmp"

# Финальная агрегация
cat "$TEMP_DIR/K_by_tissue.tmp" "$TEMP_DIR/N_by_tissue.tmp" > "$TEMP_DIR/all_counts.tmp"

echo "Пример all_counts.tmp:"
head -n 5 "$TEMP_DIR/all_counts.tmp"

awk 'BEGIN{OFS="\t"} {
    tissue = $1;
    count = $2;
    type = $3;
    
    if (type == "K") {
        K[tissue] = count;
    } else if (type == "N") {
        N[tissue] = count;
    }
    TISSUE[tissue] = tissue;
}
END {
    for (t in TISSUE) {
        n_val = (t in N) ? N[t] : 0;
        k_val = (t in K) ? K[t] : 0;
        if (n_val > 0 || k_val > 0) {
            print TISSUE[t], n_val, k_val;
        }
    }
}' "$TEMP_DIR/all_counts.tmp" | sort -t$'\t' -k2,2nr -k3,3nr > "$TEMP_DIR/data_sorted.tmp"

# Добавляем заголовок
echo -e "tissue\tN(tf+t)\tK(tf+t)" > "$OUT_DIR/tf_snv_summary_by_tissue.tsv"
cat "$TEMP_DIR/data_sorted.tmp" >> "$OUT_DIR/tf_snv_summary_by_tissue.tsv"

# Проверка: K должно быть <= N
echo "=== Проверка корректности (K <= N) ==="
awk -F'\t' 'NR>1 && $3 > $2 {print "ОШИБКА: K > N для", $1, "(K="$3", N="$2")"}' "$OUT_DIR/tf_snv_summary_by_tissue.tsv"
ERROR_COUNT=$(awk -F'\t' 'NR>1 && $3 > $2 {count++} END {print count+0}' "$OUT_DIR/tf_snv_summary_by_tissue.tsv")
if [ "$ERROR_COUNT" -eq 0 ]; then
    echo "✅ Все корректно: K <= N для всех строк"
else
    echo "⚠️ Найдено $ERROR_COUNT строк где K > N"
fi

rm -rf "$TEMP_DIR"

echo "=== Готово ==="
echo "Результат: $OUT_DIR/tf_snv_summary_by_tissue.tsv"
echo "Строк (включая заголовок): $(wc -l < "$OUT_DIR/tf_snv_summary_by_tissue.tsv")"
echo "Пример (первые 10 строк):"
head -n 10 "$OUT_DIR/tf_snv_summary_by_tissue.tsv"