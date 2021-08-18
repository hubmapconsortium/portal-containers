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

from utils import read_csv_to_pandas, get_centroid, get_type_x_antigen_df

SEGMENTATION_TYPES = ["cell", "nuclei", "cell_boundaries", "nucleus_boundaries"]
AGG_TYPES = ["mean", "total"]
SEGMENTATION_X_ANTIGEN_FILE_SUFFIX = ".ome.tiff-SEGMENTATION_TYPE_channel_AGG_TYPE.csv"
CLUSTER_FILE_SUFFIX = ".ome.tiff-SEGMENTATION_TYPE_cluster.csv"
POLYGON_FILE_SUFFUX = ".ome.tiff-cell_polygons_spatial.csv"
TSNE_FILE_SUFFIX = ".ome.tiff-tSNE_allfeatures.csv"


def get_xy(img_name, input_dir):
    """Main method for converting an input image to a numpy array of centroid coordinates

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :rtype: <class 'numpy.ndarray'>
    """
    polygon_file = Path(input_dir) / (img_name + POLYGON_FILE_SUFFUX)
    df_spatial = read_csv_to_pandas(
        # It seems like the list of points in the "Shape" column is read in as a string
        # so it needs be evaluated.
        polygon_file,
        converters={"Shape": json.loads},
    )
    df_xy = df_spatial.apply(get_centroid, axis=1).to_frame(name="Shape")
    # This is the return type that AnnData seems to like, just like the tsne method below.
    return np.array([x[0] for x in df_xy.values.tolist()])


def get_type_x_antigen_dict(img_name, input_dir):
    """Main method for converting an input image to a dict of numpy arrays, each of which will be a layer in the AnnData store,
    one per aggregation type (mean/total) and segmentation type (i.e cell_boundaries, nuclei etc.)

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :rtype: dict
    """
    segmentation_type_dict = {}
    for segmentation_type in SEGMENTATION_TYPES:
        for agg_type in AGG_TYPES:
            segmentation_quantification_df = get_type_x_antigen_df(
                img_name, input_dir, segmentation_type, agg_type
            )
            segmentation_type_dict[
                f"{segmentation_type}_x_antigen_{agg_type}"
            ] = segmentation_quantification_df.to_numpy()
    return segmentation_type_dict


def get_cluster_df(img_name, input_dir):
    """Main method for converting an input image to a merged dataframe of clusterings,
    where each column is a clustering type for a given segmentation and aggregation type
    i.e a column for "total nuclei K-Means [Covariance]"

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :rtype: <class 'pandas.core.frame.DataFrame'>
    """
    df_list = []
    for segmentation_type in SEGMENTATION_TYPES:
        cluster_file = Path(input_dir) / (
            img_name
            + CLUSTER_FILE_SUFFIX.replace("SEGMENTATION_TYPE", segmentation_type)
        )
        # Pandas reads in as int64 by default which can't be interpreted by Javascript.
        df_cluster = read_csv_to_pandas(cluster_file).astype("uint8")
        prefix = segmentation_type.replace("_", " ").title() + " "
        df_list += [df_cluster.add_prefix(prefix)]
    df_merged = reduce(lambda left, right: pd.merge(left, right, on=["ID"]), df_list)
    return df_merged


def get_antigen_labels(img_name, input_dir):
    """Main method for converting an input image to a empty pandas dataframe indexed by antigen label

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :rtype: <class 'pandas.core.frame.DataFrame'>
    """
    # Does not matter which file from which we pull the labels.
    return pd.DataFrame(
        index=get_type_x_antigen_df(img_name, input_dir, "cell", "mean").columns
    )


def get_tsne(img_name, input_dir):
    """Main method for converting an input image to a numpy array of tsne coordinates

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :rtype: <class 'numpy.ndarray'>
    """
    tsne_file = Path(input_dir) / (img_name + TSNE_FILE_SUFFIX)
    df_tsne = read_csv_to_pandas(tsne_file)
    return np.array(df_tsne.values.tolist())


def sprm_to_anndata(img_name, input_dir, output_dir):
    """Main method for processing sprm results for an image and writing AnnData

    :param str img_name: Name of the image, like R001_X001_Y001
    :param str input_dir: Path to the image
    :param str output_dir: Path for the output store
    """
    type_x_antigen_dict = get_type_x_antigen_dict(img_name, input_dir)

    adata = anndata.AnnData(
        X=type_x_antigen_dict["cell_x_antigen_mean"],
        layers=type_x_antigen_dict,
        obsm={"xy": get_xy(img_name, input_dir), "tsne": get_tsne(img_name, input_dir)},
        obs=get_cluster_df(img_name, input_dir),
        var=get_antigen_labels(img_name, input_dir),
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
