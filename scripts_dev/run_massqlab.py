from pathlib import Path
import sys

project_root = Path(__file__).resolve().parents[1]
src_path = project_root / "src_dev"
if str(src_path) not in sys.path:
    sys.path.insert(0, str(src_path))

import core_dev

DATA_DIRECTORY = Path(r"C:/Users/johnsondyj/Documents/Projects/DylanJohnson/MassQLab/data")
QUERY_FILE = Path(r"C:/Users/johnsondyj/Documents/Projects/DylanJohnson/MassQLab/MassQL_Queries.csv")

if __name__ == "__main__":
    raw_ms1_df, raw_ms2_df, area_df, run_output_dir = core_dev.massqlab_main(
        data_directory=DATA_DIRECTORY,
        query_file=QUERY_FILE,
    )

    print(f"MS1 results shape:  {raw_ms1_df.shape}")
    print(f"MS2 results shape:  {raw_ms2_df.shape}")
    print(f"Area results shape: {area_df.shape}")
    print(f"Outputs written to: {run_output_dir}")