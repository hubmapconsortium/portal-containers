import argparse
from pathlib import Path
from os import path

from mudata import read_h5mu

VAR_CHUNK_SIZE = 10
INPUT_FILE_NAMES = ["secondary_analysis.h5mu"]

def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5mu_file in INPUT_FILE_NAMES:
        # Check if input file exists, skip it if it doesn't exist
        input_path = path.join(input_dir, h5mu_file)
        if not path.exists(input_path):
            print(f"Input file {h5mu_file} does not exist in input directory.")
            continue
        mdata = read_h5mu(input_dir / h5mu_file)

        # It is now possible for adata.X to be empty and have shape (0, 0)
        # so we need to check for that here, otherwise there will
        # be a division by zero error during adata.write_zarr
        # Reference: https://github.com/hubmapconsortium/salmon-rnaseq/blob/dfb0e2a/bin/analysis/scvelo_analysis.py#L69
        chunks = (mdata.shape[0], VAR_CHUNK_SIZE) if mdata.shape[1] >= VAR_CHUNK_SIZE else None
        
        zarr_path = output_dir / (Path(h5mu_file).stem + ".zarr")
        mdata.write_zarr(zarr_path, chunks=chunks)

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
