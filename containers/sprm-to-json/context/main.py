import argparse
from glob import glob
from pathlib import Path
from os import makedirs
import json

from shapely.geometry import Polygon
import pandas as pd
import numpy as np


def sprm_to_items(input_file):
    df = pd.read_csv(input_file)
    df = df.set_index("Unnamed: 0")
    df_items = df.T.to_dict().items()
    return df_items


def write_genes_or_factors_to_cells(input_file, cells={}, genes=True):
    vitessce_key = "genes" if genes else "factors"
    df_items = sprm_to_items(input_file)
    for (cell_key, value) in df_items:
        cells[cell_key][vitessce_key] = {}
        for value_key in value:
            cells[cell_key][vitessce_key][value_key] = value[value_key]
    return cells


def write_polyon_bounds(input_file, cells):
    df_items = sprm_to_items(input_file)
    for (k, v) in df_items:
        shape = eval(v["Shape"])
        poly = Polygon(shape)
        idx = np.round(np.linspace(0, len(shape) - 1, 8)).astype(int)
        shape_downsample = np.asarray(shape)[idx]
        poly_downsample = Polygon(shape_downsample)
        cells[k] = {
            "xy": list(poly.centroid.coords),
            "poly": list(poly_downsample.exterior.coords),
        }

    return cells


def create_cells(tile_str, input_dir, output_dir):
    cells = {}
    # Get shapes.
    polygon_input = Path(input_dir) / Path(
        tile_str + ".ome.tiff-cell_polygons_spatial.csv"
    )
    cells = write_polyon_bounds(polygon_input, cells)
    genes_input = Path(input_dir) / Path(tile_str + ".ome.tiff-cell_channel_total.csv")
    cells = write_genes_or_factors_to_cells(
        input_file=genes_input, cells=cells, genes=True
    )
    cluster_input = Path(input_dir) / Path(tile_str + ".ome.tiff-cell_cluster.csv")
    cells = write_genes_or_factors_to_cells(
        input_file=cluster_input, cells=cells, genes=False
    )
    with open(Path(output_dir) / Path(tile_str).with_suffix(".cells.json"), "w") as f:
        f.write(json.dumps(cells))


def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input_file in glob(input_dir + "/*.ome.tiff-cell_polygons_spatial.csv"):
        tile_str = Path(input_file).name.split(".")[0]
        create_cells(tile_str, input_dir, output_dir)
        # create_factors(tile_str, input_dir, output_dir)
        # create_genes(tile_str, input_dir, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Transform SPRM into JSON.
        """
    )
    parser.add_argument("--input_dir", required=True, help="directory containing SPRM")
    parser.add_argument(
        "--output_dir",
        required=True,
        help="directory where arrow files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
