from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

# Safer for command-line / non-notebook execution.
matplotlib.use("Agg")

from . import core_dev
from . import ms1_analysis_dev
from . import ms1_plots
from . import ms1_traces


QUERY_EXTENSIONS = {".json", ".csv", ".xlsx", ".xls"}


def resolve_query_file(query_path: str | Path) -> Path:
    """
    Accept either:
    - a direct query file path, or
    - a directory containing exactly one query file.
    """
    query_path = Path(query_path).expanduser().resolve()

    if query_path.is_file():
        if query_path.suffix.lower() not in QUERY_EXTENSIONS:
            raise ValueError(
                f"Unsupported query file type: {query_path.suffix}. "
                f"Expected one of: {sorted(QUERY_EXTENSIONS)}"
            )
        return query_path

    if query_path.is_dir():
        matches = []
        for suffix in QUERY_EXTENSIONS:
            matches.extend(query_path.glob(f"*{suffix}"))

        matches = sorted(
            p for p in matches
            if p.is_file() and not p.name.startswith("~$")
        )

        if len(matches) == 1:
            return matches[0]

        if len(matches) == 0:
            raise FileNotFoundError(
                f"No query files found in {query_path}. "
                f"Expected one of: {sorted(QUERY_EXTENSIONS)}"
            )

        raise ValueError(
            "Multiple query files found. Please pass the exact query file path:\n"
            + "\n".join(str(p) for p in matches)
        )

    raise FileNotFoundError(f"Query path does not exist: {query_path}")


def run_workflow(
    query_file_or_directory: str | Path,
    mzml_directory: str | Path,
    output_directory: str | Path | None = None,
    min_points: int = 5,
    x_col: str = "rt",
    y_col: str = "i",
    center_method: str | None = "median",
    log_transform: bool = True,
    abs_limit: float = 1,
) -> Path:
    query_file = resolve_query_file(query_file_or_directory)
    mzml_directory = Path(mzml_directory).expanduser().resolve()

    if not mzml_directory.is_dir():
        raise NotADirectoryError(f"mzML directory does not exist: {mzml_directory}")

    mzml_files = sorted(mzml_directory.glob("*.mzML"))
    if not mzml_files:
        raise FileNotFoundError(f"No .mzML files found in: {mzml_directory}")

    print(f"Query file: {query_file}")
    print(f"mzML directory: {mzml_directory}")

    raw_ms1_df, raw_ms2_df, area_df, run_output_dir = core_dev.massqlab_main(
        data_directory=mzml_directory,
        query_file=query_file,
        output_directory=output_directory,
    )

    print(f"Base outputs written to: {run_output_dir}")

    if raw_ms1_df.empty:
        print("No MS1 results found. Skipping MS1 downstream analysis and plots.")
        return run_output_dir

    peak_area_df, peak_rt_df, peak_gaussian_df = ms1_analysis_dev.raw_ms1_df_analysis(
        raw_ms1_df,
        min_points=min_points,
        x_col=x_col,
        y_col=y_col,
        output_path=run_output_dir / "ms1_analysis_export_bundle.xlsx",
    )

    print(f"MS1 analysis workbook written to: {run_output_dir / 'ms1_analysis_export_bundle.xlsx'}")

    gaussian_cols = [
        c for c in peak_gaussian_df.columns
        if str(c).endswith(" gaussian")
    ]

    if gaussian_cols:
        n_plots = ms1_traces.plot_peak_gaussians(
            peak_gaussian_df=peak_gaussian_df,
            raw_ms1_df=raw_ms1_df,
            output_pdf=run_output_dir / "peak_gaussian_overlays.pdf",
            x_col=x_col,
            y_col=y_col,
        )

        print(f"Wrote {n_plots} Gaussian overlay plots to:")
        print(run_output_dir / "peak_gaussian_overlays.pdf")
    else:
        print("No Gaussian columns found. Skipping Gaussian overlay PDF.")

    area_cols = [c for c in area_df.columns if str(c).endswith(" area")]
    peak_area_cols = [
        c for c in peak_area_df.columns
        if str(c).endswith(" peak_area")
    ]

    if area_cols and peak_area_cols:
        ms1_plots.plot_area_and_peak_area_heatmaps(
            area_df=area_df,
            peak_area_df=peak_area_df,
            area_output_file=run_output_dir / "area_heatmap.pdf",
            peak_area_output_file=run_output_dir / "peak_area_heatmap.pdf",
            center_method=center_method,
            log_transform=log_transform,
            abs_limit=abs_limit,
        )

        print("Wrote heatmaps to:")
        print(run_output_dir / "area_heatmap.pdf")
        print(run_output_dir / "peak_area_heatmap.pdf")
    else:
        print("Missing area columns. Skipping heatmaps.")

    return run_output_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="massqlab",
        description="Run the MassQLab notebook_dev workflow from the command line.",
    )

    parser.add_argument(
        "query_file_directory",
        help="Path to a query file, or to a directory containing exactly one .json/.csv/.xlsx/.xls query file.",
    )

    parser.add_argument(
        "mzml_file_directory",
        help="Directory containing .mzML files.",
    )

    parser.add_argument(
        "--output-directory",
        default=None,
        help="Optional base output directory. Defaults to the mzML data directory.",
    )

    parser.add_argument(
        "--min-points",
        type=int,
        default=5,
        help="Minimum number of points required for MS1 Gaussian fitting.",
    )

    parser.add_argument(
        "--x-col",
        default="rt",
        help="Retention-time column name.",
    )

    parser.add_argument(
        "--y-col",
        default="i",
        help="Intensity column name.",
    )

    parser.add_argument(
        "--center-method",
        default="median",
        choices=["mean", "median", "none"],
        help="Heatmap row-centering method.",
    )

    parser.add_argument(
        "--log-transform",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Apply log2(x + 1) transform in heatmaps. Enabled by default; use --no-log-transform to disable.",
    )

    parser.add_argument(
        "--abs-limit",
        type=float,
        default=1,
        help="Absolute heatmap color scale limit.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    center_method = None if args.center_method == "none" else args.center_method

    try:
        run_output_dir = run_workflow(
            query_file_or_directory=args.query_file_directory,
            mzml_directory=args.mzml_file_directory,
            output_directory=args.output_directory,
            min_points=args.min_points,
            x_col=args.x_col,
            y_col=args.y_col,
            center_method=center_method,
            log_transform=args.log_transform,
            abs_limit=args.abs_limit,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"MassQLab workflow complete: {run_output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())