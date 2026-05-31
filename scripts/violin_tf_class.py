#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

input_file = "tf_snv_analysis_results/tf_snv_pioneer_class_family.tsv"
output_dir = "violin_plots_by_class"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_file, sep="\t", engine="python")

# если есть мусорный индекс-столбец, удаляем
first_col = str(df.columns[0]).strip()
if first_col.startswith("Unnamed") or first_col == "":
    df = df.drop(columns=df.columns[0])

# чистим имена столбцов
df.columns = [str(c).strip() for c in df.columns]

required_cols = ["TF", "is_pioneer", "K/N", "TF_class"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(
        f"В файле отсутствуют обязательные столбцы: {missing}\n"
        f"Найденные столбцы: {list(df.columns)}"
    )

df["K/N"] = pd.to_numeric(df["K/N"], errors="coerce")
df["TF_class"] = df["TF_class"].astype(str).str.strip()

def pioneer_label(x):
    x = str(x).strip().lower()
    if x == "true":
        return "ПТФ"
    return "ТФ"

df["TF_type"] = df["is_pioneer"].apply(pioneer_label)

df = df.dropna(subset=["K/N", "TF_class", "TF_type"]).copy()
df = df[df["TF_class"] != ""].copy()

# порядок категорий
order = ["ПТФ", "ТФ"]

# цвета
palette = {
    "ПТФ": "#D7AEF6",
    "ТФ": "#70E0E3"
}

sns.set_style("white")

def safe_filename(text):
    text = str(text).strip()
    text = re.sub(r"[\\/*?:\"<>|]", "_", text)
    text = re.sub(r"\s+", "_", text)
    return text

classes = sorted(df["TF_class"].unique())

print(f"Всего классов: {len(classes)}")

saved_count = 0
skipped_count = 0

for i, tf_class in enumerate(classes, start=1):
    sub = df[df["TF_class"] == tf_class].copy()

    if sub.empty:
        skipped_count += 1
        continue

    # считаем число TF в каждой группе
    counts = sub["TF_type"].value_counts()
    n_ptf = counts.get("ПТФ", 0)
    n_tf = counts.get("ТФ", 0)

    # строим график только если есть минимум 3 ПТФ и минимум 3 ТФ
    if n_ptf < 3 or n_tf < 3:
        print(
            f"Пропуск класса '{tf_class}': "
            f"ПТФ = {n_ptf}, ТФ = {n_tf} (нужно минимум 3 и 3)"
        )
        skipped_count += 1
        continue

    present_order = [x for x in order if x in set(sub["TF_type"])]

    fig, ax = plt.subplots(figsize=(5.5, 6.0))

    sns.violinplot(
        data=sub,
        x="TF_type",
        y="K/N",
        hue="TF_type",
        order=present_order,
        hue_order=present_order,
        palette=palette,
        cut=0,
        inner="box",
        linewidth=1.2,
        dodge=False,
        legend=False,
        ax=ax
    )

    sns.stripplot(
        data=sub,
        x="TF_type",
        y="K/N",
        order=present_order,
        color="gray",
        alpha=0.6,
        size=4,
        jitter=0.18,
        ax=ax
    )

    ax.set_title(tf_class, fontsize=13, pad=12, color="black")
    ax.set_xlabel("")
    ax.set_ylabel("Доля сайтов с ОНВ", fontsize=11, color="black")

    # серые линии осей
    ax.spines["bottom"].set_color("gray")
    ax.spines["left"].set_color("gray")
    ax.spines["top"].set_color("gray")
    ax.spines["right"].set_color("gray")

    # деления серые, подписи черные
    ax.tick_params(axis="x", colors="gray", labelcolor="black", labelsize=11)
    ax.tick_params(axis="y", colors="gray", labelcolor="black", labelsize=10)

    ax.xaxis.label.set_color("black")
    ax.yaxis.label.set_color("black")

    ax.grid(axis="y", color="lightgray", linestyle="--", linewidth=0.7, alpha=0.7)

    plt.tight_layout()

    out_name = safe_filename(tf_class)
    png_path = os.path.join(output_dir, f"{out_name}.png")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")

    plt.close(fig)
    saved_count += 1

    if i % 20 == 0 or i == len(classes):
        print(f"Обработано {i} из {len(classes)}")

print(f"Готово. Сохранено графиков: {saved_count}")
print(f"Пропущено классов: {skipped_count}")
print(f"Графики лежат в папке: {output_dir}")