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
    ### ИСПРАВЛЕНО: выводим TF, tissue, motif_key без подсчета
    ### Уникальность будем определять после объединения всех файлов
    bedtools intersect -a "$tf_norm" -b "$snv_exact" -wa -wb 2>/dev/null | \
    awk -v tf="$tf_name" 'BEGIN{OFS="\t"} {
        motif_key = $1":"$2":"$3;
        tissue = $10;
        print tf, tissue, motif_key;
    }' > "$out_k"
    
    # 4. Расчет N (попадание мотива в интервал SNV)
    ### ИСПРАВЛЕНО: выводим TF, tissue, motif_key без подсчета
    bedtools intersect -a "$tf_norm" -b "$snv_interval" -wa -wb 2>/dev/null | \
    awk -v tf="$tf_name" 'BEGIN{OFS="\t"} {
        motif_key = $1":"$2":"$3;
        tissue = $10;
        print tf, tissue, motif_key;
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

# === КЛЮЧЕВОЕ ИСПРАВЛЕНИЕ ===
# 1. Сначала объединяем все файлы
# 2. sort -u убирает дубликаты (один мотив в одной ткани для одного ТФ = 1 раз)
# 3. Только потом считаем количество уникальных мотивов

echo "Удаление дубликатов мотивов и подсчет..."

# Для K: убираем дубликаты TF+tissue+motif, затем считаем по TF+tissue
cat "$TEMP_DIR"/*_K.tmp 2>/dev/null | sort -u | \
awk -F'\t' '{
    key = $1 "\t" $2;
    count[key]++;
}
END {
    for (k in count) {
        print k, count[k], "K";
    }
}' OFS='\t' > "$TEMP_DIR/K_final.tmp"

# Для N: та же логика
cat "$TEMP_DIR"/*_N.tmp 2>/dev/null | sort -u | \
awk -F'\t' '{
    key = $1 "\t" $2;
    count[key]++;
}
END {
    for (k in count) {
        print k, count[k], "N";
    }
}' OFS='\t' > "$TEMP_DIR/N_final.tmp"

# Объединяем K и N
cat "$TEMP_DIR/K_final.tmp" "$TEMP_DIR/N_final.tmp" > "$TEMP_DIR/all_counts.tmp"

echo "Пример all_counts.tmp:"
head -n 5 "$TEMP_DIR/all_counts.tmp"

# Финальная агрегация
awk 'BEGIN{OFS="\t"} {
    tf = $1;
    tissue = $2;
    count = $3;
    type = $4;
    
    key = tf "\t" tissue;
    
    if (type == "K") {
        K[key] = count;
    } else if (type == "N") {
        N[key] = count;
    }
    TF[key] = tf;
    TISSUE[key] = tissue;
}
END {
    print "TF", "N(tf+t)", "K(tf+t)", "tissue";
    for (k in TF) {
        n_val = (k in N) ? N[k] : 0;
        k_val = (k in K) ? K[k] : 0;
        if (n_val > 0 || k_val > 0) {
            print TF[k], n_val, k_val, TISSUE[k];
        }
    }
}' "$TEMP_DIR/all_counts.tmp" | sort -t$'\t' -k1,1 -k4,4 > "$TEMP_DIR/data_sorted.tmp"

# Добавляем заголовок
echo -e "TF\tN(tf+t)\tK(tf+t)\ttissue" > "$OUT_DIR/tf_snv_summary.tsv"
cat "$TEMP_DIR/data_sorted.tmp" >> "$OUT_DIR/tf_snv_summary.tsv"

# === ПРОВЕРКА: K должно быть <= N ===
echo "=== Проверка корректности (K <= N) ==="
awk -F'\t' 'NR>1 && $3 > $2 {print "ОШИБКА: K > N для", $1, $4, "(K="$3", N="$2")"}' "$OUT_DIR/tf_snv_summary.tsv"
ERROR_COUNT=$(awk -F'\t' 'NR>1 && $3 > $2 {count++} END {print count+0}' "$OUT_DIR/tf_snv_summary.tsv")
if [ "$ERROR_COUNT" -eq 0 ]; then
    echo "✅ Все корректно: K <= N для всех строк"
else
    echo "⚠️ Найдено $ERROR_COUNT строк где K > N"
fi

rm -rf "$TEMP_DIR"

echo "=== Готово ==="
echo "Результат: $OUT_DIR/tf_snv_summary.tsv"
echo "Строк: $(wc -l < "$OUT_DIR/tf_snv_summary.tsv")"
echo "Пример (первые 10 строк):"
head -n 10 "$OUT_DIR/tf_snv_summary.tsv"