
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence
import numpy as np
import pandas as pd
from pandas.api.types import is_scalar
import pyopenms as oms
from massql import msql_engine
import contextlib
import io


# ---------------------------
# Query normalization helpers
# ---------------------------

def _validate_single_query_source(
    query_file: str | Path | None = None,
    query_dict: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    query_list: Sequence[str] | None = None,
) -> str:
    provided = [
        name
        for name, value in {
            "query_file": query_file,
            "query_dict": query_dict,
            "query_list": query_list,
        }.items()
        if value is not None
    ]

    if len(provided) != 1:
        raise ValueError(
            "Exactly one of query_file, query_dict, or query_list must be provided."
        )

    return provided[0]


def infer_mslevel(query: str) -> int:
    q = query.upper()
    return 2 if "MS2DATA" in q else 1


def _load_queries_from_file(query_file: str | Path) -> list[dict[str, Any]]:
    query_path = Path(query_file)
    suffix = query_path.suffix.lower()

    if suffix == ".json":
        with open(query_path, "r", encoding="utf-8") as fh:
            queries = json.load(fh)
    elif suffix == ".csv":
        queries = pd.read_csv(query_path).to_dict("records")
    elif suffix in {".xlsx", ".xls"}:
        queries = pd.read_excel(query_path).to_dict("records")
    else:
        raise ValueError(f"Unsupported query file type: {suffix}")

    if not isinstance(queries, list):
        raise ValueError("Query file must contain a list-like collection of queries.")

    return queries


def _normalize_query_entry(entry: Any, idx: int) -> dict[str, Any]:
    """
    Normalize a single query entry while preserving any extra metadata fields.
    """
    if isinstance(entry, str):
        query_text = entry.strip()
        if not query_text:
            return {}

        return {
            "name": f"query_{idx}",
            "query": query_text,
            "mslevel": infer_mslevel(query_text),
        }

    if isinstance(entry, Mapping):
        normalized = dict(entry)

        query_text = str(normalized.get("query", "")).strip()
        if not query_text:
            return {}

        raw_name = normalized.get("name")
        if raw_name is None or str(raw_name).strip() == "" or pd.isna(raw_name):
            name = f"query_{idx}"
        else:
            name = str(raw_name).strip()

        raw_mslevel = normalized.get("mslevel")
        if raw_mslevel is None or str(raw_mslevel).strip() == "" or pd.isna(raw_mslevel):
            mslevel = infer_mslevel(query_text)
        else:
            mslevel = int(raw_mslevel)

        normalized["name"] = name
        normalized["query"] = query_text
        normalized["mslevel"] = mslevel
        return normalized

    raise TypeError(f"Unsupported query entry type: {type(entry)}")


def massql_query_resolver(
    query_file: str | Path | None = None,
    query_dict: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    query_list: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """
    Normalize query input into a list of dicts and preserve all query metadata fields.
    Enforces that exactly one source is provided.
    """
    source_name = _validate_single_query_source(query_file, query_dict, query_list)

    if source_name == "query_file":
        raw_queries = _load_queries_from_file(query_file)
    elif source_name == "query_dict":
        raw_queries = [query_dict] if isinstance(query_dict, Mapping) else list(query_dict)
    else:
        raw_queries = list(query_list)

    normalized_queries: list[dict[str, Any]] = []
    for idx, entry in enumerate(raw_queries, start=1):
        normalized = _normalize_query_entry(entry, idx)
        if normalized:
            normalized_queries.append(normalized)

    if not normalized_queries:
        raise ValueError("No valid queries were found.")

    names = [q["name"] for q in normalized_queries]
    duplicate_names = pd.Series(names)[pd.Series(names).duplicated()].unique().tolist()
    if duplicate_names:
        raise ValueError(
            "Query names must be unique. Duplicate query names found: "
            + ", ".join(map(str, duplicate_names))
        )

    return normalized_queries


def build_query_metadata_df(queries: Sequence[dict[str, Any]]) -> pd.DataFrame:
    """
    Build a metadata table from normalized queries.
    This is the single source of truth for query metadata throughout the project.
    """
    if not queries:
        return pd.DataFrame(columns=["name", "query", "mslevel"])

    metadata_df = pd.DataFrame(queries).rename(columns={"name": "name"})

    ordered_cols = ["name", "query", "mslevel"] + [
        c for c in metadata_df.columns if c not in {"name", "query", "mslevel"}
    ]
    return metadata_df.loc[:, ordered_cols]


def build_result_metadata(query: Mapping[str, Any]) -> dict[str, Any]:
    """
    Convert normalized query metadata into result metadata.
    Preserves all extra fields from massql_query_resolver.
    """
    metadata = {
        "name": query["name"],
        "query": query["query"],
        "mslevel": query["mslevel"],
    }

    for key, value in query.items():
        if key not in {"name", "query", "mslevel"}:
            metadata[key] = value

    return metadata


# ---------------------------
# DataFrame helpers
# ---------------------------

def combine_frames(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """
    Efficiently combine a sequence of DataFrames.
    - 0 frames: empty DataFrame
    - 1 frame : return that frame with reset index
    - many    : single pd.concat call

    This is the best in-memory strategy here because massql_engine already
    returns DataFrames. The main improvement is reducing the number of concat inputs.
    """
    if not frames:
        return pd.DataFrame()
    if len(frames) == 1:
        return frames[0].reset_index(drop=True)
    return pd.concat(frames, ignore_index=True)


def attach_metadata_columns(
    df: pd.DataFrame,
    metadata: Mapping[str, Any],
    source_file: str,
    source_path: str,
) -> pd.DataFrame:
    """
    Attach metadata to a query result DataFrame.
    """
    out = df.copy()

    out["source_file"] = source_file
    out["source_path"] = source_path

    for col, value in metadata.items():
        if is_scalar(value):
            out[col] = value
        else:
            out[col] = [value] * len(out)

    return out


# ---------------------------
# mzML loading and execution
# ---------------------------

def _load_massql_data(filename: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    exp = oms.MSExperiment()
    mzml_file = oms.MzMLFile()
    mzml_file.setLogType(oms.LogType.CMD)
    mzml_file.load(str(filename), exp)

    ms1_df, ms2_df = exp.get_massql_df()
    return ms1_df, ms2_df


def run_massql_query(
    query_text: str,
    filename: str,
    ms1_df: pd.DataFrame,
    ms2_df: pd.DataFrame,
    cache: Any = None,
    suppress_output: bool = True,
) -> pd.DataFrame:
    """
    Run a MassQL query with optional stdout/stderr suppression.
    """
    kwargs = {
        "cache": cache,
        "ms1_df": ms1_df,
        "ms2_df": ms2_df,
    }

    if suppress_output:
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()

        with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
            return msql_engine.process_query(query_text, filename, **kwargs)

    return msql_engine.process_query(query_text, filename, **kwargs)

    
def _run_queries_on_file(
    filename: str | Path,
    queries: Sequence[dict[str, Any]],
    cache: Any = None,
    suppress_query_output: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load one mzML once, run all queries against it, and combine results once per file.
    """
    filename = Path(filename)
    source_file = filename.name
    source_path = str(filename)

    ms1_df, ms2_df = _load_massql_data(filename)

    file_ms1_frames: list[pd.DataFrame] = []
    file_ms2_frames: list[pd.DataFrame] = []

    for query in queries:
        results_df = run_massql_query(
            query_text=query["query"],
            filename=str(filename),
            ms1_df=ms1_df,
            ms2_df=ms2_df,
            cache=cache,
            suppress_output=suppress_query_output,
        )

        if results_df is None or results_df.empty:
            continue

        metadata = build_result_metadata(query)
        results_df = attach_metadata_columns(
            results_df,
            metadata=metadata,
            source_file=source_file,
            source_path=source_path,
        )

        if query["mslevel"] == 1:
            file_ms1_frames.append(results_df)
        elif query["mslevel"] == 2:
            file_ms2_frames.append(results_df)
        else:
            raise ValueError(
                f"Unsupported mslevel={query['mslevel']} for query '{query['name']}'"
            )

    return combine_frames(file_ms1_frames), combine_frames(file_ms2_frames)


def massqlab_execute(
    filenames: Iterable[str | Path],
    queries: Sequence[dict[str, Any]],
    cache: Any = None,
    suppress_query_output: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Execute all queries across all mzML files.
    """
    all_ms1_file_frames: list[pd.DataFrame] = []
    all_ms2_file_frames: list[pd.DataFrame] = []

    for filename in filenames:
        file_ms1, file_ms2 = _run_queries_on_file(
            filename=filename,
            queries=queries,
            cache=cache,
            suppress_query_output=suppress_query_output,
        )

        if not file_ms1.empty:
            all_ms1_file_frames.append(file_ms1)
        if not file_ms2.empty:
            all_ms2_file_frames.append(file_ms2)

    raw_df_ms1 = combine_frames(all_ms1_file_frames)
    raw_df_ms2 = combine_frames(all_ms2_file_frames)

    return raw_df_ms1, raw_df_ms2


# ---------------------------
# Area summarization
# ---------------------------

def massqlab_area(
    massqlab_raw_ms1: pd.DataFrame,
    queries: Sequence[dict[str, Any]],
) -> pd.DataFrame:
    """
    Collapse raw MS1 results to one row per name and one column per source_file,
    summing the 'i' column.

    Keeps only:
    - name
    - query metadata from build_query_metadata_df
    - one column per source_file containing summed intensity

    Includes only name values present in massqlab_raw_ms1.
    """
    metadata_df = build_query_metadata_df(queries)
    metadata_cols = ["name"] + [c for c in metadata_df.columns if c != "name"]

    if massqlab_raw_ms1.empty:
        return pd.DataFrame(columns=metadata_cols)

    required_cols = {"source_file", "name", "i"}
    missing = required_cols - set(massqlab_raw_ms1.columns)
    if missing:
        raise KeyError(
            f"massqlab_raw_ms1 is missing required columns for area calculation: {sorted(missing)}"
        )

    area_df = (
        massqlab_raw_ms1
        .groupby(["name", "source_file"], sort=False, as_index=False)["i"]
        .sum()
        .pivot(index="name", columns="source_file", values="i")
        .reset_index()
        .merge(
            metadata_df,
            on="name",
            how="left",
            validate="one_to_one",
        )
    )

    area_df.columns.name = None

    source_file_cols = [c for c in area_df.columns if c not in metadata_cols]
    area_df = area_df.rename(columns={c: f"{c} area" for c in source_file_cols})

    renamed_source_file_cols = [f"{c} area" for c in source_file_cols]
    area_df[renamed_source_file_cols] = area_df[renamed_source_file_cols].fillna(0)
    area_df[renamed_source_file_cols] = np.log2(area_df[renamed_source_file_cols] + 1)

    return area_df.loc[:, metadata_cols + renamed_source_file_cols]


# ---------------------------
# Output helpers
# ---------------------------

def prepare_output_directory(
    data_directory: str | Path,
    output_directory: str | Path | None = None,
) -> Path:
    base_dir = Path(output_directory) if output_directory is not None else Path(data_directory)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    run_dir = base_dir / f"MassQLab_Output" / f"massqlab_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def write_dataframe(
    df: pd.DataFrame,
    path_stem: Path,
    write_parquet: bool = False,
) -> Path:
    """
    Write CSV by default.
    If write_parquet=True, write only parquet.

    Returns the written file path.
    """
    if write_parquet:
        output_path = path_stem.with_suffix(".parquet")
        df.to_parquet(output_path, index=False)
    else:
        output_path = path_stem.with_suffix(".csv")
        df.to_csv(output_path, index=False)

    return output_path


def write_queries_json(queries: Sequence[dict[str, Any]], output_path: Path) -> Path:
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(list(queries), fh, indent=2, ensure_ascii=False, default=str)
    return output_path


def _queries_to_dataframe(queries: Any) -> pd.DataFrame:
    """
    Convert queries into a tabular form suitable for an Excel sheet.
    """
    if isinstance(queries, pd.DataFrame):
        return queries.copy()

    if isinstance(queries, dict):
        return pd.json_normalize(queries, sep=".")

    if isinstance(queries, list):
        if len(queries) == 0:
            return pd.DataFrame()

        if all(isinstance(x, dict) for x in queries):
            return pd.json_normalize(queries, sep=".")

        return pd.DataFrame({"value": queries})

    return pd.DataFrame({"value": [queries]})


def write_excel_bundle(
    queries: Any,
    raw_ms1_df: pd.DataFrame,
    raw_ms2_df: pd.DataFrame,
    area_df: pd.DataFrame,
    output_path: Path,
) -> Path:
    """
    Write all outputs into a single Excel workbook with multiple sheets.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    queries_df = _queries_to_dataframe(queries)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        queries_df.to_excel(writer, sheet_name="normalized_queries", index=False)
        raw_ms1_df.to_excel(writer, sheet_name="massqlab_raw_ms1", index=False)
        raw_ms2_df.to_excel(writer, sheet_name="massqlab_raw_ms2", index=False)
        area_df.to_excel(writer, sheet_name="massqlab_area", index=False)

    return output_path

    
# ---------------------------
# Main
# ---------------------------

def massqlab_main(
    data_directory: str | Path,
    query_file: str | Path | None = None,
    query_dict: Mapping[str, Any] | Sequence[Mapping[str, Any]] | None = None,
    query_list: Sequence[str] | None = None,
    output_directory: str | Path | None = None,
    write_parquet: bool = False,
    cache: Any = None,
    suppress_query_output: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Path]:
    """
    Run MassQLab workflow and write outputs to a timestamped subdirectory.
    """
    _validate_single_query_source(query_file, query_dict, query_list)

    data_directory = Path(data_directory)
    filenames = sorted(data_directory.glob("*.mzML"))
    if not filenames:
        raise FileNotFoundError(f"No .mzML files found in: {data_directory}")

    queries = massql_query_resolver(
        query_file=query_file,
        query_dict=query_dict,
        query_list=query_list,
    )

    raw_ms1_df, raw_ms2_df = massqlab_execute(
        filenames=filenames,
        queries=queries,
        cache=cache,
        suppress_query_output=suppress_query_output,
    )

    area_df = massqlab_area(raw_ms1_df, queries)

    run_output_dir = prepare_output_directory(
        data_directory=data_directory,
        output_directory=output_directory,
    )

    excel_path = run_output_dir / "massqlab_export_bundle.xlsx"

    write_excel_bundle(
        queries=queries,
        raw_ms1_df=raw_ms1_df,
        raw_ms2_df=raw_ms2_df,
        area_df=area_df,
        output_path=excel_path,
    )
    
    return raw_ms1_df, raw_ms2_df, area_df, run_output_dir
