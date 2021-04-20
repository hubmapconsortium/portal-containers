import argparse
from glob import glob
from pathlib import Path
from os import makedirs
import json

import anndata
import zarr
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

CELL_X_ANTIGENS_FILE_SUFFIX = ".ome.tiff-cell_channel_mean.csv"
CLUSTER_FILE_SUFFIX = ".ome.tiff-cell_cluster.csv"
POLYGON_FILE_SUFFUX = ".ome.tiff-cell_polygons_spatial.csv"
NUM_VERTICES = 10


def downsample_shape(shape):
    idx = np.round(np.linspace(0, len(shape[0]) - 1, NUM_VERTICES + 1)).astype(int)
    shape_downsample = np.asarray(shape[0])[idx]
    poly_downsample = Polygon(shape_downsample).exterior.coords
    return list(poly_downsample)


def get_centroid(shape):
    poly = Polygon(np.asarray(shape[0]))
    # poly.centroid.coords is an object that should be converted to a list, but when we get the 0th
    # coordinate i.e the centroid, it is a tuple which needs to be a list for further processing.
    return list(list(poly.centroid.coords)[0])


def read_sprm_to_pandas(input_file, converters={}):
    df = pd.read_csv(input_file, converters=converters).set_index("ID")
    return df


def sprm_to_anndata(img_name, input_dir, output_dir):
    cell_x_antigen_file = Path(input_dir) / (img_name + CELL_X_ANTIGENS_FILE_SUFFIX)
    df_cell_x_antigen = read_sprm_to_pandas(cell_x_antigen_file)

    cluster_file = Path(input_dir) / (img_name + CLUSTER_FILE_SUFFIX)
    # Pandas reads in as int64 by default which can't be interpreted by Javascript.
    df_cluster = read_sprm_to_pandas(cluster_file).astype("uint8")

    polygon_file = Path(input_dir) / (img_name + POLYGON_FILE_SUFFUX)
    df_spatial = read_sprm_to_pandas(
        # It seems like the list of points in the "Shape" column is read in as a string
        # so it needs be evaluated.
        polygon_file,
        converters={"Shape": json.loads},
    )
    # Get centroid and downsampled polygon coordinates.
    df_poly = df_spatial.apply(downsample_shape, axis=1).to_frame(name="Shape")
    df_xy = df_spatial.apply(get_centroid, axis=1).to_frame(name="Shape")

    adata = anndata.AnnData(
        X=df_cell_x_antigen.to_numpy(),
        obsm={
            "poly": np.array(df_poly.values.tolist()),
            "xy": np.array([x[0] for x in df_xy.values.tolist()]),
        },
        obs=df_cluster,
        var=pd.DataFrame(index=df_cell_x_antigen.columns),
    )
    adata.write_zarr(str(output_dir / (img_name + "-anndata.zarr")))


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for input_file in input_dir.glob(f"*{POLYGON_FILE_SUFFUX}"):
        img_name = Path(input_file).name.split(".")[0]
        sprm_to_anndata(img_name, input_dir, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Transform SPRM into AnnData.
        """
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing SPRM .csv files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where AnnData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
