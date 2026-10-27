import argparse
from pathlib import Path

from mudata import read_h5mu

from create_zarr_stores import create_zarr_for_masks
INPUT_FILE_NAMES = ["objects.h5mu"]

def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5mu_file in INPUT_FILE_NAMES:
        # Raises if the file is missing, rather than leaving an empty output directory.
        mdata = read_h5mu(input_dir / h5mu_file)
        zarr_path = output_dir / (Path(h5mu_file).stem + ".zarr")
        create_zarr_for_masks(mdata, zarr_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Transform Segmentation Mudata into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing MuData .h5mu files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where seg-muData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
