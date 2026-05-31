#!/bin/bash
set -euo pipefail

# Новые входные директории
dir_udacha_dnase="/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/dnase"
dir_udacha_atac="/shared_data/rscf_scatac-seq_2024-2025/exp_wo_treat/atac"

# Выходной файл
OUTFILE="merged_coords_001.bed"
TMPFILE="tmp_snps_unsorted.bed"

# Очистка временных и итоговых файлов
> "$TMPFILE"
> "$OUTFILE"

# Функция обработки одного файла (только данные УДАЧИ)
process_file() {
    local file="$1"
    awk -F'\t' 'BEGIN { OFS="\t" }
    # Фильтрация: rsID в колонке 4, все необходимые поля не пустые
    $4 ~ /^rs/ && $1 != "" && $2 != "" && $3 != "" && $5 != "" && $6 != "" {
        # Формат вывода: chr start end chr:end:ref:alt 0 +
        print $1, $2, $3, $1":"$3":"$5":"$6, 0, "+"
    }' "$file" >> "$TMPFILE"
}

# Обработка всех BED-файлов из обеих директорий
for dir in "$dir_udacha_atac" "$dir_udacha_dnase"; do
    if [ -d "$dir" ]; then
        for f in "$dir"/*.bed; do
            [ -e "$f" ] || continue
            echo "Processing $f"
            process_file "$f"
        done
    else
        echo "Warning: directory not found: $dir" >&2
    fi
done

# Удаление дубликатов и сортировка по хромосоме и позиции
if [ -s "$TMPFILE" ]; then
    sort -k1,1 -k2,2n "$TMPFILE" | uniq > "$OUTFILE"
else
    echo "Warning: no SNPs extracted, output file will be empty" >&2
fi

# Удаление временного файла
rm -f "$TMPFILE"

echo "Готово: результат сохранён в $OUTFILE"