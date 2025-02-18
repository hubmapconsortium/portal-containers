import argparse
from pathlib import Path
from os import path

import zarr
from scipy.io import mmread
import pandas as pd
from vitessce import SnapWrapper

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
SECONDARY_ANALYSIS = "secondary_analysis.h5ad"
SCVELO_ANNOTATED = "scvelo_annotated.h5ad"

def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    mtx = mmread(path.join(input_dir, 'filtered_cell_by_bin.mtx'))
    barcodes_df = pd.read_csv(path.join(input_dir, 'barcodes.txt'), header=None)
    bins_df = pd.read_csv(path.join(input_dir, 'bins.txt'), header=None)
    clusters_df = pd.read_csv(path.join(input_dir, 'umap_coords_clusters.csv'), index_col=0)
    zarr_filepath = path.join(output_dir, 'hubmap-ui.snap.multires.zarr')
    w = SnapWrapper(mtx, barcodes_df, bins_df, clusters_df)
    # In theory, it would be nice to create an AnnData store instead of json
    # We could then attach things to it, like clusters in `uns`
    # w.create_anndata()
    w.create_genomic_multivec_zarr(zarr_filepath)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Transform Snap into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing HuBMAP SnapATAC data",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where (AnnData) zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
