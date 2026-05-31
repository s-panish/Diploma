#!/usr/bin/env Rscript
# создаются неуникальныее строки
# === Загрузка библиотек ===
suppressPackageStartupMessages({
  library(motifbreakR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(MotifDb)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
  library(BiocParallel)
})

# === Параметры ===
bed_file <- "merged_coords_001.bed"
output_file <- "motifbreakr_from_bed_001_log.tsv"
motif_file <- "motifdb_hocomocov13_motifs.rds"
chunk_size <- 4000
start_chunk <- 1

# === 1. Загрузка BED ===
bed_data <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
colnames(bed_data) <- c("chr", "start", "end", "name", "score", "strand")

# === 2. Сохраняем как временный BED-файл ===
temp_bed <- tempfile(fileext = ".bed")
write.table(bed_data, file = temp_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# === 3. Импорт SNP с rsID из dbSNP ===
variants <- snps.from.file(
  file = temp_bed,
  format = "bed",
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  check.unnamed.for.rsid = TRUE
)
cat("Импортировано SNP с rsID (если доступны):", length(variants), "\n")

# === 4. Загрузка PWM ===
motifs <- readRDS(motif_file)
cat("Загружено PWM из .rds:", length(motifs), "\n")

# === 5. Разбиение на чанки и поэтапная обработка ===
genome(variants) <- "hg38"
total <- length(variants)
n_chunks <- ceiling(total / chunk_size)

for (i in seq(from = start_chunk, to = n_chunks)) {
  cat(sprintf("Обработка чанка %d из %d...\n", i, n_chunks))

  idx_start <- ((i - 1) * chunk_size) + 1
  idx_end <- min(i * chunk_size, total)
  chunk <- variants[idx_start:idx_end]

  chunk_results <- motifbreakR(
    snpList = chunk,
    pwmList = motifs,
    threshold = 0.9995,
    method = "log",
    filterp = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = 110)
  )

  if (length(chunk_results) == 0) {
    cat("Нет результатов в чанке", i, "\n")
    next
  }

  df <- as.data.frame(chunk_results, row.names = NULL)
  df[] <- lapply(df, function(col) if (is.list(col)) sapply(col, toString) else col)

  write.table(df, file = output_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE, append = TRUE)

  cat(sprintf("Чанк %d сохранён: %d строк\n", i, nrow(df)))
}

cat("Обработка завершена. Результаты добавлены в", output_file, "\n")
