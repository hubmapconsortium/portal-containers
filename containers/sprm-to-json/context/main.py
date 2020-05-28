import argparse
from glob import glob
from pathlib import Path
from os import makedirs
import json

from shapely.geometry import Polygon
import pandas as pd
import numpy as np

# Number of vertices for the polygon boundary data.
NUM_VERTICES = 20
GENES_FILE_SUFFIX = ".ome.tiff-cell_channel_mean.csv"
CLUSTER_FILE_SUFFIX = ".ome.tiff-cell_cluster.csv"
POLYGON_FILE_SUFFUX = ".ome.tiff-cell_polygons_spatial.csv"


def sprm_to_items(input_file):
    df = pd.read_csv(input_file)
    df = df.set_index("ID")
    df_items = df.T.to_dict().items()
    return df_items


def write_genes_or_factors_to_cells(input_file, cells={}, genes=True):
    vitessce_key = "genes" if genes else "factors"
    df_items = sprm_to_items(input_file)
    for (cell_key, value) in df_items:
        cells[cell_key][vitessce_key] = {}
        for value_key in value:
            cells[cell_key][vitessce_key][value_key] = (
                value[value_key] if genes else str(value[value_key])
            )
    return cells


def write_polyon_bounds(input_file, cells):
    df_items = sprm_to_items(input_file)
    for (k, v) in df_items:
        shape = eval(v["Shape"])
        poly = Polygon(shape)
        idx = np.round(np.linspace(0, len(shape) - 1, NUM_VERTICES + 1)).astype(int)
        shape_downsample = np.asarray(shape)[idx]
        poly_downsample = Polygon(shape_downsample)
        cells[k] = {
            "xy": list(poly.centroid.coords)[0],
            "poly": list(poly_downsample.exterior.coords),
        }

    return cells


def create_cells(tile_str, input_dir, output_dir):
    cells = {}
    # Get shapes.
    polygon_input = Path(input_dir) / Path(tile_str + POLYGON_FILE_SUFFUX)
    cells = write_polyon_bounds(polygon_input, cells)
    # Write genes to cells file.
    genes_input = Path(input_dir) / Path(tile_str + GENES_FILE_SUFFIX)
    cells = write_genes_or_factors_to_cells(
        input_file=genes_input, cells=cells, genes=True
    )
    # Write factors to cells file.
    cluster_input = Path(input_dir) / Path(tile_str + CLUSTER_FILE_SUFFIX)
    cells = write_genes_or_factors_to_cells(
        input_file=cluster_input, cells=cells, genes=False
    )
    with open(Path(output_dir) / Path(tile_str).with_suffix(".cells.json"), "w") as f:
        f.write(json.dumps(cells))


def create_factors(tile_str, input_dir, output_dir):
    factors = {}
    cluster_input = Path(input_dir) / Path(tile_str + CLUSTER_FILE_SUFFIX)
    df = pd.read_csv(cluster_input).set_index("ID")
    cluster_types = df.columns.values
    df_items = sprm_to_items(cluster_input)
    for cluster_type in cluster_types:
        cluster_names = sorted(df[cluster_type].unique().astype("uint8"))
        factors[cluster_type] = {
            "map": [str(cluster) for cluster in cluster_names],
            "cells": {k: v[cluster_type] for (k, v) in df_items},
        }
    with open(Path(output_dir) / Path(tile_str).with_suffix(".factors.json"), "w") as f:
        f.write(json.dumps(factors))


def create_genes(tile_str, input_dir, output_dir):
    genes = {}
    genes_input = Path(input_dir) / Path(tile_str + GENES_FILE_SUFFIX)
    df = pd.read_csv(genes_input).set_index("ID")
    gene_types = df.columns.values
    df_items = sprm_to_items(genes_input)
    for gene_type in gene_types:
        genes[gene_type] = {
            "max": max(df[gene_type]),
            "cells": {k: v[gene_type] for (k, v) in df_items},
        }
    with open(Path(output_dir) / Path(tile_str).with_suffix(".genes.json"), "w") as f:
        f.write(json.dumps(genes))


def create_clusters(tile_str, input_dir, output_dir):
    clusters = {}
    genes_input = Path(input_dir) / Path(tile_str + GENES_FILE_SUFFIX)
    df = pd.read_csv(genes_input).set_index("ID")
    gene_types = df.columns.values
    cell_ids = df.index.values
    clusters["rows"] = [str(gene_type) for gene_type in gene_types]
    clusters["cols"] = [str(cell_id) for cell_id in cell_ids]
    clusters["matrix"] = []
    for gene_type in gene_types:
        gene_col = df[gene_type]
        max_min_diff = max(gene_col) - min(gene_col)
        # If max_min difference is 0, just output a column of 0's.
        normalized = (
            (gene_col - min(gene_col)) / max_min_diff
            if max_min_diff != 0
            else (gene_col - gene_col)
        )
        clusters["matrix"].append(list(normalized.values))
    with open(
        Path(output_dir) / Path(tile_str).with_suffix(".clusters.json"), "w"
    ) as f:
        f.write(json.dumps(clusters))


def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input_file in glob(input_dir + "/*.ome.tiff-cell_polygons_spatial.csv"):
        tile_str = Path(input_file).name.split(".")[0]
        # Create cells schema.
        create_cells(tile_str, input_dir, output_dir)
        # Create factors schema.
        create_factors(tile_str, input_dir, output_dir)
        # Create genes schema.
        create_genes(tile_str, input_dir, output_dir)
        # Create clusters schema.
        create_clusters(tile_str, input_dir, output_dir)


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
