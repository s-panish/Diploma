#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path


# ================= НАСТРОЙКИ =================

INPUT_FILE = "total_tf_snv_pioneer_class_family.tsv"

OUTPUT_SUMMARY = "permutation_test_pioneer_vs_nonpioneer_summary.tsv"
OUTPUT_PERMUTATIONS = "permutation_test_pioneer_vs_nonpioneer_null_distribution.tsv"

N_PERMUTATIONS = 100000
RANDOM_SEED = 42

# Варианты:
# "two-sided"  — проверяет любое отличие между группами
# "greater"    — проверяет, что у пионерных ТФ K/N выше
# "less"       — проверяет, что у пионерных ТФ K/N ниже
ALTERNATIVE = "two-sided"

# =============================================


def parse_bool_column(series: pd.Series) -> pd.Series:
    """
    Приводит колонку is_pioneer к булевому типу.
    Поддерживает значения True/False, TRUE/FALSE, 1/0, yes/no.
    """
    if series.dtype == bool:
        return series

    mapping = {
        "true": True,
        "false": False,
        "1": True,
        "0": False,
        "yes": True,
        "no": False,
        "y": True,
        "n": False,
    }

    parsed = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
        .map(mapping)
    )

    if parsed.isna().any():
        bad_values = sorted(series[parsed.isna()].astype(str).unique())
        raise ValueError(
            "В колонке is_pioneer есть значения, которые не удалось распознать: "
            + ", ".join(bad_values)
        )

    return parsed.astype(bool)


def calculate_p_value(observed_stat: float, permuted_stats: np.ndarray, alternative: str) -> float:
    """
    Расчёт p-value по эмпирическому нулевому распределению.
    Используется поправка +1, чтобы p-value не было равно нулю.
    """
    if alternative == "two-sided":
        extreme_count = np.sum(np.abs(permuted_stats) >= abs(observed_stat))
    elif alternative == "greater":
        extreme_count = np.sum(permuted_stats >= observed_stat)
    elif alternative == "less":
        extreme_count = np.sum(permuted_stats <= observed_stat)
    else:
        raise ValueError("ALTERNATIVE должен быть: 'two-sided', 'greater' или 'less'")

    p_value = (extreme_count + 1) / (len(permuted_stats) + 1)
    return p_value


def main():
    print("1. Загрузка данных...")

    input_path = Path(INPUT_FILE)

    if not input_path.exists():
        raise FileNotFoundError(f"Файл не найден: {INPUT_FILE}")

    df = pd.read_csv(input_path, sep="\t", engine="python")

    # Удаляем возможный индексный столбец, если он появился при сохранении таблицы
    first_col = str(df.columns[0]).strip()
    if first_col.startswith("Unnamed") or first_col == "":
        df = df.drop(columns=df.columns[0])

    # Чистим имена колонок
    df.columns = [str(c).strip() for c in df.columns]

    required_cols = ["TF", "N(tf)", "K(tf)", "is_pioneer", "K/N", "TF_class", "TF_family"]
    missing_cols = [c for c in required_cols if c not in df.columns]

    if missing_cols:
        raise ValueError(
            "В таблице отсутствуют обязательные колонки: "
            + ", ".join(missing_cols)
        )

    print(f"   Всего строк: {len(df)}")
    print(f"   Уникальных TF: {df['TF'].nunique()}")

    print("2. Подготовка колонок...")

    df["is_pioneer"] = parse_bool_column(df["is_pioneer"])

    df["K/N"] = pd.to_numeric(df["K/N"], errors="coerce")
    df["N(tf)"] = pd.to_numeric(df["N(tf)"], errors="coerce")
    df["K(tf)"] = pd.to_numeric(df["K(tf)"], errors="coerce")

    before_filter = len(df)

    df = df.dropna(subset=["TF", "is_pioneer", "K/N", "N(tf)", "K(tf)"]).copy()

    # На всякий случай исключаем строки с некорректными знаменателями
    df = df[df["N(tf)"] > 0].copy()

    after_filter = len(df)

    print(f"   Строк до фильтрации NA/некорректных N(tf): {before_filter}")
    print(f"   Строк после фильтрации: {after_filter}")

    # Если K/N мог быть сохранён округлённым, можно пересчитать его из K(tf)/N(tf).
    # Это полезно, если нужна максимальная точность.
    df["K_over_N_recalculated"] = df["K(tf)"] / df["N(tf)"]

    # Основной анализ проводим по пересчитанной доле.
    # Если нужно использовать именно готовую колонку K/N, замените analysis_col на "K/N".
    analysis_col = "K_over_N_recalculated"

    n_pioneer = int(df["is_pioneer"].sum())
    n_nonpioneer = int((~df["is_pioneer"]).sum())

    if n_pioneer == 0:
        raise ValueError("В таблице нет пионерных ТФ: is_pioneer == True")

    if n_nonpioneer == 0:
        raise ValueError("В таблице нет не-пионерных ТФ: is_pioneer == False")

    print("3. Описание групп...")

    pioneer_values = df.loc[df["is_pioneer"], analysis_col].to_numpy()
    nonpioneer_values = df.loc[~df["is_pioneer"], analysis_col].to_numpy()

    observed_median_pioneer = float(np.median(pioneer_values))
    observed_median_nonpioneer = float(np.median(nonpioneer_values))
    observed_mean_pioneer = float(np.mean(pioneer_values))
    observed_mean_nonpioneer = float(np.mean(nonpioneer_values))

    observed_diff_median = observed_median_pioneer - observed_median_nonpioneer
    observed_diff_mean = observed_mean_pioneer - observed_mean_nonpioneer

    print(f"   Пионерные ТФ: {n_pioneer}")
    print(f"   Не-пионерные ТФ: {n_nonpioneer}")
    print(f"   Медиана K/N у пионерных ТФ: {observed_median_pioneer:.6f}")
    print(f"   Медиана K/N у не-пионерных ТФ: {observed_median_nonpioneer:.6f}")
    print(f"   Разность медиан: {observed_diff_median:.6f}")
    print(f"   Среднее K/N у пионерных ТФ: {observed_mean_pioneer:.6f}")
    print(f"   Среднее K/N у не-пионерных ТФ: {observed_mean_nonpioneer:.6f}")
    print(f"   Разность средних: {observed_diff_mean:.6f}")

    print("4. Пермутационный тест...")

    rng = np.random.default_rng(RANDOM_SEED)

    all_values = df[analysis_col].to_numpy()
    n_total = len(all_values)

    permuted_diff_median = np.empty(N_PERMUTATIONS, dtype=float)
    permuted_diff_mean = np.empty(N_PERMUTATIONS, dtype=float)

    for i in range(N_PERMUTATIONS):
        permuted_indices = rng.permutation(n_total)

        random_pioneer_indices = permuted_indices[:n_pioneer]
        random_nonpioneer_indices = permuted_indices[n_pioneer:]

        random_pioneer_values = all_values[random_pioneer_indices]
        random_nonpioneer_values = all_values[random_nonpioneer_indices]

        permuted_diff_median[i] = (
            np.median(random_pioneer_values)
            - np.median(random_nonpioneer_values)
        )

        permuted_diff_mean[i] = (
            np.mean(random_pioneer_values)
            - np.mean(random_nonpioneer_values)
        )

        if (i + 1) % 10000 == 0:
            print(f"   Выполнено перестановок: {i + 1}/{N_PERMUTATIONS}")

    p_value_median = calculate_p_value(
        observed_stat=observed_diff_median,
        permuted_stats=permuted_diff_median,
        alternative=ALTERNATIVE
    )

    p_value_mean = calculate_p_value(
        observed_stat=observed_diff_mean,
        permuted_stats=permuted_diff_mean,
        alternative=ALTERNATIVE
    )

    print("5. Сохранение результатов...")

    summary = pd.DataFrame([
        {
            "test": "permutation_test",
            "comparison": "pioneer_vs_nonpioneer",
            "alternative": ALTERNATIVE,
            "n_permutations": N_PERMUTATIONS,
            "random_seed": RANDOM_SEED,
            "analysis_column": analysis_col,
            "n_total_tf": n_total,
            "n_pioneer_tf": n_pioneer,
            "n_nonpioneer_tf": n_nonpioneer,

            "median_pioneer": observed_median_pioneer,
            "median_nonpioneer": observed_median_nonpioneer,
            "observed_diff_median": observed_diff_median,
            "p_value_median": p_value_median,

            "mean_pioneer": observed_mean_pioneer,
            "mean_nonpioneer": observed_mean_nonpioneer,
            "observed_diff_mean": observed_diff_mean,
            "p_value_mean": p_value_mean,

            "null_median_diff_mean": float(np.mean(permuted_diff_median)),
            "null_median_diff_sd": float(np.std(permuted_diff_median, ddof=1)),
            "null_median_diff_2.5_percentile": float(np.percentile(permuted_diff_median, 2.5)),
            "null_median_diff_97.5_percentile": float(np.percentile(permuted_diff_median, 97.5)),

            "null_mean_diff_mean": float(np.mean(permuted_diff_mean)),
            "null_mean_diff_sd": float(np.std(permuted_diff_mean, ddof=1)),
            "null_mean_diff_2.5_percentile": float(np.percentile(permuted_diff_mean, 2.5)),
            "null_mean_diff_97.5_percentile": float(np.percentile(permuted_diff_mean, 97.5)),
        }
    ])

    null_distribution = pd.DataFrame({
        "permutation_id": np.arange(1, N_PERMUTATIONS + 1),
        "permuted_diff_median": permuted_diff_median,
        "permuted_diff_mean": permuted_diff_mean,
    })

    summary.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)
    null_distribution.to_csv(OUTPUT_PERMUTATIONS, sep="\t", index=False)

    print(f"   Сводная таблица сохранена: {OUTPUT_SUMMARY}")
    print(f"   Нулевое распределение сохранено: {OUTPUT_PERMUTATIONS}")

    print("6. Итог:")
    print(f"   Разность медиан K/N: {observed_diff_median:.6f}")
    print(f"   p-value для разности медиан: {p_value_median:.6g}")
    print(f"   Разность средних K/N: {observed_diff_mean:.6f}")
    print(f"   p-value для разности средних: {p_value_mean:.6g}")


if __name__ == "__main__":
    main()