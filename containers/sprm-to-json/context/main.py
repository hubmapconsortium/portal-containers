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
    df = pd.read_csv(input_file).set_index("ID")
    df_items = df.T.to_dict().items()
    return df_items


def write_genes_or_factors_to_cells(input_file, cells={}, is_genes=True):
    vitessce_key = "genes" if is_genes else "factors"
    df_items = sprm_to_items(input_file)
    for (cell_key, value) in df_items:
        cells[cell_key][vitessce_key] = {}
        for value_key in value:
            # Genes look for a string while factors do not.
            cells[cell_key][vitessce_key][value_key] = (
                value[value_key] if is_genes else str(value[value_key])
            )
    return cells


def create_factors_or_genes(input_dir, tile_str, is_genes=True):
    suffix = GENES_FILE_SUFFIX if is_genes else CLUSTER_FILE_SUFFIX
    input_path = Path(input_dir) / f"{tile_str}{suffix}"
    df = pd.read_csv(input_path).set_index("ID")
    types = df.columns.values
    df_items = sprm_to_items(input_path)
    return (df_items, df, types)


def write_polyon_bounds(input_file, cells):
    df_items = sprm_to_items(input_file)
    for (k, v) in df_items:
        shape = json.loads(v["Shape"])
        poly = Polygon(shape)
        idx = np.round(np.linspace(0, len(shape) - 1, NUM_VERTICES + 1)).astype(int)
        shape_downsample = np.asarray(shape)[idx]
        poly_downsample = Polygon(shape_downsample).exterior.coords
        cells[k] = {
            "xy": list(poly.centroid.coords)[0],
            "poly": list(poly_downsample.exterior.coords),
        }

    return cells


def create_cells(tile_str, input_dir, output_dir):
    cells = {}
    # Get shapes.
    polygon_input = Path(input_dir) / f"{tile_str}{POLYGON_FILE_SUFFUX}"
    cells = write_polyon_bounds(polygon_input, cells)
    # Write genes to cells file.
    genes_input = Path(input_dir) / Path(tile_str + GENES_FILE_SUFFIX)
    cells = write_genes_or_factors_to_cells(
        input_file=genes_input, cells=cells, is_genes=True
    )
    # Write factors to cells file.
    cluster_input = Path(input_dir) / Path(tile_str + CLUSTER_FILE_SUFFIX)
    cells = write_genes_or_factors_to_cells(
        input_file=cluster_input, cells=cells, is_genes=False
    )
    with open(Path(output_dir) / f"{tile_str}.cells.json", "w") as f:
        f.write(json.dumps(cells, indent=4))


def create_factors(tile_str, input_dir, output_dir):
    factors = {}
    (df_items, df, cluster_types) = create_factors_or_genes(
        input_dir, tile_str, is_genes=False
    )
    for cluster_type in cluster_types:
        cluster_names = sorted(df[cluster_type].unique().astype("uint8"))
        factors[cluster_type] = {
            "map": [str(cluster) for cluster in cluster_names],
            "cells": {k: v[cluster_type] for (k, v) in df_items},
        }
    with open(Path(output_dir) / f"{tile_str}.factors.json", "w") as f:
        f.write(json.dumps(factors, indent=4))


def create_cell_sets(tile_str, input_dir, output_dir):
    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
    cell_sets = {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": []
    }

    (df_items, df, cluster_types) = create_factors_or_genes(
        input_dir, tile_str, is_genes=False
    )

    for cluster_type in cluster_types:
        cluster_names = sorted(df[cluster_type].unique().astype("uint8"))
        cell_sets["tree"].append({
            "name": cluster_type,
            "children": [
                {
                    "name": f"Cluster {cluster}",
                    "set": df.loc[df[cluster_type] == cluster].index.astype(str).values.tolist()
                }
                for cluster in cluster_names
            ]
        })
    with open(Path(output_dir) / f"{tile_str}.cell-sets.json", "w") as f:
        f.write(json.dumps(cell_sets, indent=4))


def create_genes(tile_str, input_dir, output_dir):
    genes = {}
    (df_items, df, gene_types) = create_factors_or_genes(
        input_dir, tile_str, is_genes=True
    )
    for gene_type in gene_types:
        genes[gene_type] = {
            "max": max(df[gene_type]),
            "cells": {k: v[gene_type] for (k, v) in df_items},
        }
    with open(Path(output_dir) / f"{tile_str}.genes.json", "w") as f:
        f.write(json.dumps(genes, indent=4))


def create_clusters(tile_str, input_dir, output_dir):
    clusters = {}
    genes_input = Path(input_dir) / f"{tile_str}{GENES_FILE_SUFFIX}"
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
            else 0 * gene_col
        )
        clusters["matrix"].append(list(normalized.values))
    with open(Path(output_dir) / f"{tile_str}.clusters.json", "w") as f:
        f.write(json.dumps(clusters, indent=4))


def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input_file in glob(input_dir + "/*.ome.tiff-cell_polygons_spatial.csv"):
        tile_str = Path(input_file).name.split(".")[0]
        # Create cells schema.
        create_cells(tile_str, input_dir, output_dir)
        # Create factors schema.
        create_factors(tile_str, input_dir, output_dir)
        # Create cell sets schema.
        create_cell_sets(tile_str, input_dir, output_dir)
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
