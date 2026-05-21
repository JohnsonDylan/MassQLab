import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_area_heatmap(
    df,
    value_suffix,
    output_file,
    row_label_col="name",
    query_col="query",
    log_transform=True,
    center_method="mean",
    figsize=None,
    cmap="vlag",
    abs_limit=8,
    na_color="#d3d3d3",
):
    value_cols = [c for c in df.columns if c.endswith(value_suffix)]
    if not value_cols:
        raise ValueError(f"No columns found ending with {value_suffix!r}")

    missing_label_cols = [c for c in [row_label_col, query_col] if c not in df.columns]
    if missing_label_cols:
        raise KeyError(f"Missing required label columns: {missing_label_cols}")

    if center_method not in ["mean", "median", None]:
        raise ValueError("center_method must be 'mean', 'median', or None")

    plot_df = df[[row_label_col, query_col] + value_cols].copy()

    row_labels = (plot_df[row_label_col].astype(str))

    matrix_df = plot_df[value_cols].copy()

    renamed_cols = [c[: -len(value_suffix)] for c in value_cols]
    matrix_df.columns = renamed_cols
    matrix_df.index = row_labels

    if log_transform:
        matrix_df = np.log2(matrix_df + 1)

    if center_method == "mean":
        matrix_df = matrix_df.sub(matrix_df.mean(axis=1, skipna=True), axis=0)
        cbar_label = "row mean-centered log2(area + 1)"
    elif center_method == "median":
        matrix_df = matrix_df.sub(matrix_df.median(axis=1, skipna=True), axis=0)
        cbar_label = "row median-centered log2(area + 1)"
    else:
        cbar_label = "log2(area + 1)" if log_transform else "area"

    if figsize is None:
        width = max(6, 0.45 * len(matrix_df.columns) + 2)
        height = max(6, 0.22 * len(matrix_df.index) + 2)
        figsize = (width, height)

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color=na_color)

    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        matrix_df,
        cmap=cmap_obj,
        center=0,
        vmin=-abs_limit,
        vmax=abs_limit,
        mask=matrix_df.isna(),
        linewidths=0,
        linecolor=None,
        cbar_kws={"label": cbar_label},
    )

    ax.set_xlabel("Source file")
    ax.set_ylabel("Name | Query")
    ax.set_title(os.path.splitext(os.path.basename(output_file))[0])

    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def plot_area_and_peak_area_heatmaps(
    area_df,
    peak_area_df,
    area_output_file="area_heatmap.pdf",
    peak_area_output_file="peak_area_heatmap.pdf",
    row_label_col="name",
    query_col="query",
    log_transform=True,
    center_method="mean",
    abs_limit=8,
    na_color="#d3d3d3",
):
    plot_area_heatmap(
        df=area_df,
        value_suffix=" area",
        output_file=area_output_file,
        row_label_col=row_label_col,
        query_col=query_col,
        log_transform=log_transform,
        center_method=center_method,
        abs_limit=abs_limit,
        na_color=na_color,
    )

    plot_area_heatmap(
        df=peak_area_df,
        value_suffix=" peak_area",
        output_file=peak_area_output_file,
        row_label_col=row_label_col,
        query_col=query_col,
        log_transform=log_transform,
        center_method=center_method,
        abs_limit=abs_limit,
        na_color=na_color,
    )