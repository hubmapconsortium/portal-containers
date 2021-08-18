import argparse
from glob import glob
from pathlib import Path
from os import makedirs
import json
from functools import reduce

import anndata
import zarr
import pandas as pd
import numpy as np
from shapely.geometry import Polygon

QUANTIFICATION_TYPES = ["cell", "nuclei", "cell_boundaries", "nucleus_boundaries"]
AGG_MODES = ['mean', 'total']
TYPE_X_ANTIGEN_FILE_SUFFIX = ".ome.tiff-TYPE_channel_AGG_MODE.csv"
TYPE_X_ANTIGEN_FILE_SUFFIX = ".ome.tiff-TYPE_channel_AGG_MODE.csv"
CLUSTER_FILE_SUFFIX = ".ome.tiff-TYPE_cluster.csv"
POLYGON_FILE_SUFFUX = ".ome.tiff-cell_polygons_spatial.csv"
TSNE_FILE_SUFFIX = ".ome.tiff-tSNE_allfeatures.csv"
NUM_VERTICES = 8


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

def get_xy(img_name, input_dir):
    polygon_file = Path(input_dir) / (img_name + POLYGON_FILE_SUFFUX)
    df_spatial = read_sprm_to_pandas(
        # It seems like the list of points in the "Shape" column is read in as a string
        # so it needs be evaluated.
        polygon_file,
        converters={"Shape": json.loads},
    )
    df_xy = df_spatial.apply(get_centroid, axis=1).to_frame(name="Shape")
    return np.array([x[0] for x in df_xy.values.tolist()])

def get_type_x_antigen_df(img_name, input_dir, quantification_type, agg_mode):
    quantification_path = Path(input_dir) / (img_name + TYPE_X_ANTIGEN_FILE_SUFFIX.replace("TYPE", quantification_type).replace("AGG_MODE", agg_mode))
    quantification_df = read_sprm_to_pandas(quantification_path)
    return quantification_df

def get_type_x_antigen_dict(img_name, input_dir):
    quantification_type_dict = {}
    for quantification_type in QUANTIFICATION_TYPES:
        for agg_mode in AGG_MODES:
            quantification_df = get_type_x_antigen_df(img_name, input_dir, quantification_type, agg_mode)
            quantification_type_dict[f"{quantification_type}_x_antigen_{agg_mode}"] = quantification_df.to_numpy()
    return quantification_type_dict

def get_cluster_df(img_name, input_dir):
    df_list = []
    for quantification_type in QUANTIFICATION_TYPES:
        cluster_file = Path(input_dir) / (img_name + CLUSTER_FILE_SUFFIX.replace("TYPE", quantification_type))
        # Pandas reads in as int64 by default which can't be interpreted by Javascript.
        df_cluster = read_sprm_to_pandas(cluster_file).astype("uint8")
        prefix = quantification_type.replace('_', ' ').title() + ' '
        df_list += [df_cluster.add_prefix(prefix)]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['ID']), df_list)
    return df_merged

def get_antigen_labels(img_name, input_dir):
    # Does not matter which file from which we pull the labels.
    return get_type_x_antigen_df(img_name, input_dir, "cell", "mean").columns

def get_tsne(img_name, input_dir):
    tsne_file = Path(input_dir) / (img_name + TSNE_FILE_SUFFIX)
    df_tsne = read_sprm_to_pandas(tsne_file)
    return np.array(df_tsne.values.tolist())

def sprm_to_anndata(img_name, input_dir, output_dir):
    df_cluster = get_cluster_df(img_name, input_dir)
    type_x_antigen_dict = get_type_x_antigen_dict(img_name, input_dir)

    adata = anndata.AnnData(
        X=type_x_antigen_dict["cell_x_antigen_mean"],
        layers=type_x_antigen_dict,
        obsm={
            "xy": get_xy(img_name, input_dir),
            "tsne": get_tsne(img_name, input_dir)
        },
        obs=df_cluster,
        var=pd.DataFrame(index=get_antigen_labels(img_name, input_dir)),
    )
    adata.write_zarr(str(output_dir / (img_name + "-anndata.zarr")))


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for input_file in input_dir.glob(f"*.ome.tiff-cell_polygons_spatial.csv"):
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
