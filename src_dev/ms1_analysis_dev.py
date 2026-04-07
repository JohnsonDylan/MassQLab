import re
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


import re
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths


def fit_single_trace(grouped_df, min_points=5, x_col="rt", y_col="i"):
    source_file = grouped_df["source_file"].iloc[0]
    name = grouped_df["name"].iloc[0]
    query = grouped_df["query"].iloc[0]

    result = {
        "source_file": source_file,
        "name": name,
        "query": query,
        "datapoint_count": 0,
        "fit_success": False,
        "fit_status": None,
        "baseline": np.nan,
        "amplitude": np.nan,
        "center_rt": np.nan,
        "sigma": np.nan,
        "fwhm": np.nan,
        "peak_area": np.nan,
        "rtmin_query": np.nan,
        "rtmax_query": np.nan,
        "rt_error": np.nan,
        "gaussian_formula": np.nan,
    }

    df = grouped_df[[x_col, y_col]].copy()
    df = df.replace([np.inf, -np.inf], np.nan).dropna()

    if df.empty:
        result["fit_status"] = "no_data"
        return result

    df = (
        df.groupby(x_col, as_index=False)[y_col]
        .mean()
        .sort_values(x_col)
        .reset_index(drop=True)
    )

    x = df[x_col].to_numpy(dtype=float)
    y = df[y_col].to_numpy(dtype=float)
    result["datapoint_count"] = len(x)

    query_lower = query.lower() if isinstance(query, str) else ""
    rtmin_match = re.search(r"rtmin\s*[:=]\s*([0-9]*\.?[0-9]+)", query_lower)
    rtmax_match = re.search(r"rtmax\s*[:=]\s*([0-9]*\.?[0-9]+)", query_lower)

    if rtmin_match:
        result["rtmin_query"] = float(rtmin_match.group(1))
    if rtmax_match:
        result["rtmax_query"] = float(rtmax_match.group(1))

    rtmin = result["rtmin_query"]
    rtmax = result["rtmax_query"]

    if len(x) < min_points:
        result["fit_status"] = "insufficient_points"
        return result

    baseline = float(np.min(y))
    y_corr = y - baseline
    y_corr[y_corr < 0] = 0
    result["baseline"] = baseline

    if np.max(y_corr) <= 0:
        result["fit_status"] = "flat_signal"
        return result

    prominence = max(np.max(y_corr) / 10.0, 1e-12)
    height = max(np.mean(y_corr) * 1.1, 1e-12)

    try:
        peaks, properties = find_peaks(
            y_corr,
            height=height,
            distance=5,
            prominence=prominence,
        )
    except Exception:
        result["fit_status"] = "peak_detection_failed"
        return result

    if len(peaks) == 0:
        result["fit_status"] = "no_peak_found"
        return result

    if pd.notna(rtmin) and pd.notna(rtmax):
        range_rt = rtmax - rtmin
        filtered_peaks = [
            p for p in peaks
            if (rtmin - range_rt) < x[p] < (rtmax + range_rt)
        ]
    else:
        filtered_peaks = list(peaks)

    if len(filtered_peaks) == 0:
        result["fit_status"] = "no_peak_in_window"
        return result

    try:
        widths, _, l_ips, r_ips = peak_widths(y_corr, filtered_peaks, rel_height=0.5)
    except Exception:
        result["fit_status"] = "width_estimation_failed"
        return result

    base_idx = np.arange(len(x), dtype=float)
    left_x = np.interp(l_ips, base_idx, x)
    right_x = np.interp(r_ips, base_idx, x)
    real_width = right_x - left_x

    if pd.notna(rtmin) and pd.notna(rtmax):
        range_rt = rtmax - rtmin
        final_peaks = [
            p for i, p in enumerate(filtered_peaks)
            if np.isfinite(real_width[i]) and real_width[i] > 0 and real_width[i] < range_rt
        ]
    else:
        final_peaks = [
            p for i, p in enumerate(filtered_peaks)
            if np.isfinite(real_width[i]) and real_width[i] > 0
        ]

    if len(final_peaks) == 0:
        result["fit_status"] = "no_valid_peak_width"
        return result

    peak = max(final_peaks, key=lambda p: y_corr[p])

    try:
        w, _, l, r = peak_widths(y_corr, [peak], rel_height=0.5)
    except Exception:
        result["fit_status"] = "width_estimation_failed"
        return result

    left = np.interp(l[0], base_idx, x)
    right = np.interp(r[0], base_idx, x)
    fwhm = float(right - left)

    if not np.isfinite(fwhm) or fwhm <= 0:
        result["fit_status"] = "invalid_fwhm"
        return result

    center_rt = float((left + right) / 2.0)
    sigma = float(fwhm / (2 * np.sqrt(2 * np.log(2))))
    amplitude = float(y[peak] - baseline)

    if not np.isfinite(amplitude) or amplitude <= 0:
        result["fit_status"] = "invalid_amplitude"
        return result

    result["fit_success"] = True
    result["fit_status"] = "ok"
    result["amplitude"] = amplitude
    result["center_rt"] = center_rt
    result["sigma"] = sigma
    result["fwhm"] = fwhm
    result["peak_area"] = amplitude * sigma * np.sqrt(2 * np.pi)
    result["gaussian_formula"] = (
        f"{baseline:.12g} + {amplitude:.12g} * exp(-0.5 * ((x - {center_rt:.12g}) / {sigma:.12g})^2)"
    )

    if pd.notna(rtmin) and pd.notna(rtmax):
        if rtmin <= center_rt <= rtmax:
            result["rt_error"] = 0.0
        elif center_rt < rtmin:
            result["rt_error"] = center_rt - rtmin
        else:
            result["rt_error"] = center_rt - rtmax

    return result


def raw_ms1_df_analysis(raw_ms1_df, min_points=5, x_col="rt", y_col="i"):
    long_results = []

    for (source_file, name, query), grouped_df in raw_ms1_df.groupby(
        ["source_file", "name", "query"],
        sort=False
    ):
        result = {
            "source_file": source_file,
            "name": name,
            "query": query,
            "peak_area": np.nan,
            "center_rt": np.nan,
            "rt_error": np.nan,
            "gaussian_formula": np.nan,
            "fit_success": False,
            "fit_status": None,
            "datapoint_count": 0,
        }

        df = grouped_df[[x_col, y_col]].copy()
        df = df.replace([np.inf, -np.inf], np.nan).dropna()

        if df.empty:
            result["fit_status"] = "no_data"
            long_results.append(result)
            continue

        df = (
            df.groupby(x_col, as_index=False)[y_col]
            .mean()
            .sort_values(x_col)
            .reset_index(drop=True)
        )

        x = df[x_col].to_numpy(dtype=float)
        y = df[y_col].to_numpy(dtype=float)
        result["datapoint_count"] = len(x)

        if len(x) < min_points:
            result["fit_status"] = "insufficient_points"
            long_results.append(result)
            continue

        query_lower = query.lower() if isinstance(query, str) else ""
        rtmin_match = re.search(r"rtmin\s*[:=]\s*([0-9]*\.?[0-9]+)", query_lower)
        rtmax_match = re.search(r"rtmax\s*[:=]\s*([0-9]*\.?[0-9]+)", query_lower)

        rtmin = float(rtmin_match.group(1)) if rtmin_match else np.nan
        rtmax = float(rtmax_match.group(1)) if rtmax_match else np.nan

        baseline = float(np.min(y))
        y_corr = y - baseline
        y_corr[y_corr < 0] = 0

        if np.max(y_corr) <= 0:
            result["fit_status"] = "flat_signal"
            long_results.append(result)
            continue

        try:
            prominence = max(np.max(y_corr) / 10.0, 1e-12)
            height = max(np.mean(y_corr) * 1.1, 1e-12)

            peaks, _ = find_peaks(
                y_corr,
                height=height,
                distance=5,
                prominence=prominence,
            )
        except Exception:
            result["fit_status"] = "peak_detection_failed"
            long_results.append(result)
            continue

        if len(peaks) == 0:
            result["fit_status"] = "no_peak_found"
            long_results.append(result)
            continue

        if pd.notna(rtmin) and pd.notna(rtmax):
            range_rt = rtmax - rtmin
            filtered_peaks = [
                p for p in peaks
                if (rtmin - range_rt) < x[p] < (rtmax + range_rt)
            ]
        else:
            filtered_peaks = list(peaks)

        if len(filtered_peaks) == 0:
            result["fit_status"] = "no_peak_in_window"
            long_results.append(result)
            continue

        try:
            widths, _, l_ips, r_ips = peak_widths(y_corr, filtered_peaks, rel_height=0.5)
        except Exception:
            result["fit_status"] = "width_estimation_failed"
            long_results.append(result)
            continue

        base_idx = np.arange(len(x), dtype=float)
        left_x = np.interp(l_ips, base_idx, x)
        right_x = np.interp(r_ips, base_idx, x)
        real_width = right_x - left_x

        if pd.notna(rtmin) and pd.notna(rtmax):
            range_rt = rtmax - rtmin
            final_peaks = [
                p for i, p in enumerate(filtered_peaks)
                if np.isfinite(real_width[i]) and real_width[i] > 0 and real_width[i] < range_rt
            ]
        else:
            final_peaks = [
                p for i, p in enumerate(filtered_peaks)
                if np.isfinite(real_width[i]) and real_width[i] > 0
            ]

        if len(final_peaks) == 0:
            result["fit_status"] = "no_valid_peak_width"
            long_results.append(result)
            continue

        peak = max(final_peaks, key=lambda p: y_corr[p])

        try:
            w, _, l, r = peak_widths(y_corr, [peak], rel_height=0.5)
        except Exception:
            result["fit_status"] = "width_estimation_failed"
            long_results.append(result)
            continue

        left = np.interp(l[0], base_idx, x)
        right = np.interp(r[0], base_idx, x)
        fwhm = float(right - left)

        if not np.isfinite(fwhm) or fwhm <= 0:
            result["fit_status"] = "invalid_fwhm"
            long_results.append(result)
            continue

        center_rt = float((left + right) / 2.0)
        sigma = float(fwhm / (2 * np.sqrt(2 * np.log(2))))
        amplitude = float(y[peak] - baseline)

        if not np.isfinite(amplitude) or amplitude <= 0:
            result["fit_status"] = "invalid_amplitude"
            long_results.append(result)
            continue

        peak_area = amplitude * sigma * np.sqrt(2 * np.pi)

        rt_error = np.nan
        if pd.notna(rtmin) and pd.notna(rtmax):
            if rtmin <= center_rt <= rtmax:
                rt_error = 0.0
            elif center_rt < rtmin:
                rt_error = center_rt - rtmin
            else:
                rt_error = center_rt - rtmax

        result["peak_area"] = peak_area
        result["center_rt"] = center_rt
        result["rt_error"] = rt_error
        result["gaussian_formula"] = (
            f"{baseline:.12g} + {amplitude:.12g} * exp(-0.5 * ((x - {center_rt:.12g}) / {sigma:.12g})^2)"
        )
        result["fit_success"] = True
        result["fit_status"] = "ok"

        long_results.append(result)

    long_df = pd.DataFrame(long_results)

    if long_df.empty:
        peak_area_df = pd.DataFrame(columns=["name", "query"])
        peak_rt_df = pd.DataFrame(columns=["name", "query"])
        peak_gaussian_df = pd.DataFrame(columns=["name", "query"])
        return peak_area_df, peak_rt_df, peak_gaussian_df

    peak_area_df = (
        long_df.pivot(index=["name", "query"], columns="source_file", values="peak_area")
        .reset_index()
    )
    peak_area_df.columns = [
        col if col in ["name", "query"] else f"{col} peak_area"
        for col in peak_area_df.columns
    ]

    peak_gaussian_df = (
        long_df.pivot(index=["name", "query"], columns="source_file", values="gaussian_formula")
        .reset_index()
    )
    peak_gaussian_df.columns = [
        col if col in ["name", "query"] else f"{col} gaussian"
        for col in peak_gaussian_df.columns
    ]

    peak_rt_value_df = (
        long_df.pivot(index=["name", "query"], columns="source_file", values="center_rt")
        .reset_index()
    )
    peak_rt_value_df.columns = [
        col if col in ["name", "query"] else f"{col} peak_rt"
        for col in peak_rt_value_df.columns
    ]

    peak_rt_error_df = (
        long_df.pivot(index=["name", "query"], columns="source_file", values="rt_error")
        .reset_index()
    )
    peak_rt_error_df.columns = [
        col if col in ["name", "query"] else f"{col} rt_error"
        for col in peak_rt_error_df.columns
    ]

    peak_rt_df = peak_rt_value_df.merge(
        peak_rt_error_df,
        on=["name", "query"],
        how="outer"
    )
    
    source_files = list(long_df["source_file"].drop_duplicates())
    
    peak_rt_cols = [
        f"{source_file} peak_rt"
        for source_file in source_files
        if f"{source_file} peak_rt" in peak_rt_df.columns
    ]
    
    rt_error_cols = [
        f"{source_file} rt_error"
        for source_file in source_files
        if f"{source_file} rt_error" in peak_rt_df.columns
    ]
    
    peak_rt_df = peak_rt_df.loc[:, ["name", "query"] + peak_rt_cols + rt_error_cols]

    return peak_area_df, peak_rt_df, peak_gaussian_df