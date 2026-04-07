import re
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_peak_gaussians(
    peak_gaussian_df: pd.DataFrame,
    raw_ms1_df: pd.DataFrame,
    output_pdf: str,
    x_col: str = "rt",
    y_col: str = "i",
):
    required_peak_cols = {"name", "query"}
    required_raw_cols = {"source_file", "name", "query", x_col, y_col}

    missing_peak = required_peak_cols - set(peak_gaussian_df.columns)
    missing_raw = required_raw_cols - set(raw_ms1_df.columns)

    if missing_peak:
        raise KeyError(f"peak_gaussian_df is missing required columns: {sorted(missing_peak)}")
    if missing_raw:
        raise KeyError(f"raw_ms1_df is missing required columns: {sorted(missing_raw)}")

    gaussian_cols = [c for c in peak_gaussian_df.columns if c.endswith(" gaussian")]
    if not gaussian_cols:
        raise ValueError("peak_gaussian_df has no columns ending with ' gaussian'")

    raw_plot_df = raw_ms1_df[["source_file", "name", "query", x_col, y_col]].copy()
    raw_plot_df = raw_plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col])

    trace_dict = {}
    for (source_file, name, query), grouped_df in raw_plot_df.groupby(
        ["source_file", "name", "query"],
        sort=False
    ):
        trace_df = (
            grouped_df.groupby(x_col, as_index=False)[y_col]
            .mean()
            .sort_values(x_col)
            .reset_index(drop=True)
        )
        trace_dict[(source_file, name, query)] = trace_df

    number_pattern = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
    gaussian_pattern = re.compile(
        rf"^\s*{number_pattern}\s*\+\s*{number_pattern}\s*\*\s*exp\("
        rf"-0\.5\s*\*\s*\(\(x\s*-\s*{number_pattern}\)\s*/\s*{number_pattern}\)\^2\)\s*$"
    )

    rtmin_pattern = re.compile(r"rtmin\s*[:=]\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)
    rtmax_pattern = re.compile(r"rtmax\s*[:=]\s*([0-9]*\.?[0-9]+)", re.IGNORECASE)

    plot_count = 0

    with mpl.rc_context({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "lines.linewidth": 1.2,
        "savefig.bbox": "tight",
    }):
        with PdfPages(output_pdf) as pdf:
            for _, row in peak_gaussian_df.iterrows():
                name = row["name"]
                query = row["query"]

                rtmin = np.nan
                rtmax = np.nan
                if isinstance(query, str):
                    rtmin_match = rtmin_pattern.search(query)
                    rtmax_match = rtmax_pattern.search(query)
                    if rtmin_match:
                        rtmin = float(rtmin_match.group(1))
                    if rtmax_match:
                        rtmax = float(rtmax_match.group(1))

                for gaussian_col in gaussian_cols:
                    formula = row[gaussian_col]

                    if pd.isna(formula):
                        continue

                    match = gaussian_pattern.match(str(formula))
                    if match is None:
                        continue

                    baseline, amplitude, center_rt, sigma = map(float, match.groups())
                    source_file = gaussian_col.rsplit(" gaussian", 1)[0]

                    trace_df = trace_dict.get((source_file, name, query))
                    if trace_df is None or trace_df.empty:
                        continue

                    x = trace_df[x_col].to_numpy(dtype=float)
                    y = trace_df[y_col].to_numpy(dtype=float)

                    if len(x) == 0 or not np.isfinite(sigma) or sigma <= 0:
                        continue

                    raw_x_min = float(np.min(x))
                    raw_x_max = float(np.max(x))

                    if pd.notna(rtmin) and pd.notna(rtmax) and rtmax > rtmin:
                        range_rt = rtmax - rtmin
                        x_plot_min = rtmin - 0.5 * range_rt
                        x_plot_max = rtmax + 0.5 * range_rt
                    else:
                        x_plot_min = center_rt - 4 * sigma
                        x_plot_max = center_rt + 4 * sigma

                    if x_plot_max <= x_plot_min:
                        x_plot_min = raw_x_min
                        x_plot_max = raw_x_max

                    mask = (x >= x_plot_min) & (x <= x_plot_max)
                    if np.any(mask):
                        x_plot = x[mask]
                        y_plot = y[mask]
                    else:
                        x_plot = x
                        y_plot = y

                    x_fit = np.linspace(x_plot_min, x_plot_max, 500)
                    y_fit = baseline + amplitude * np.exp(-0.5 * ((x_fit - center_rt) / sigma) ** 2)

                    fig, ax = plt.subplots(figsize=(4.8, 3.6))

                    ax.plot(
                        x_plot,
                        y_plot,
                        color="black",
                        linewidth=0.9,
                        marker="o",
                        markersize=3.2,
                        label="Raw trace",
                    )
                    ax.plot(
                        x_fit,
                        y_fit,
                        color="red",
                        linewidth=1.4,
                        label="Gaussian",
                    )

                    if pd.notna(rtmin):
                        ax.axvline(
                            rtmin,
                            color="black",
                            linestyle="--",
                            linewidth=0.9,
                        )
                    if pd.notna(rtmax):
                        ax.axvline(
                            rtmax,
                            color="black",
                            linestyle="--",
                            linewidth=0.9,
                        )

                    ax.set_xlim(x_plot_min, x_plot_max)
                    ax.set_title(f"{name} | {source_file}")
                    ax.set_xlabel("Retention time")
                    ax.set_ylabel("Intensity")

                    ax.spines["top"].set_visible(False)
                    ax.spines["right"].set_visible(False)
                    ax.margins(x=0.02, y=0.08)

                    text_lines = [
                        f"center={center_rt:.4f}",
                        f"sigma={sigma:.4f}",
                        f"area={amplitude * sigma * np.sqrt(2 * np.pi):.4g}",
                    ]
                    if pd.notna(rtmin):
                        text_lines.append(f"RTMIN={rtmin:.4f}")
                    if pd.notna(rtmax):
                        text_lines.append(f"RTMAX={rtmax:.4f}")

                    ax.text(
                        0.02,
                        0.98,
                        "\n".join(text_lines),
                        transform=ax.transAxes,
                        ha="left",
                        va="top",
                        fontsize=7,
                    )

                    ax.legend(frameon=False, loc="best")
                    fig.tight_layout(pad=0.6)
                    pdf.savefig(fig)
                    plt.close(fig)

                    plot_count += 1

    return plot_count