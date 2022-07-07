from pathlib import Path
from typing import Callable, Dict, Optional

from shapely.geometry import Polygon
import pandas as pd
import numpy as np

SEGMENTATION_X_ANTIGEN_FILE_SUFFIX = ".ome.tiff-SEGMENTATION_TYPE_channel_AGG_TYPE.csv"


def get_centroid(shape):
    """Function for generating a centroid from polygon coordinates
    :param pandas.core.series.Series shape: A pandas series with one entry, Shape, like [[0, 1], [0, 2]]
    :rtype: list The centroid, like [0, 1]
    """
    poly = Polygon(np.asarray(shape[0]))
    # poly.centroid.coords is an object that should be converted to a list, but when we get the 0th
    # coordinate i.e the centroid, it is a tuple which needs to be a list for further processing.
    return list(list(poly.centroid.coords)[0])


def read_csv_to_pandas(input_file: Path, converters: Optional[Dict[str, Callable]] = None) -> pd.DataFrame:
    """Function for reading in a csv file to a dataframe
    :param Path input_file: File path
    :param dict converters: Per columen converters
    :rtype: pandas.core.frame.DataFrame
    """
    if converters is None:
        converters = {}
    df = pd.read_csv(input_file, converters=converters).set_index("ID")
    return df


def get_type_x_antigen_df(img_name: str, input_dir: Path, segmentation_type: str, agg_type: str) -> pd.DataFrame:
    """Converts an input image to a dataframe
    for a given aggregation type (mean/total) and segmentation type (i.e cell_boundaries, nuclei etc.)

    :param str img_name: Name of the image, like R001_X001_Y001
    :param Path input_dir: Path to the image
    :param str segmentation_type: Like cell or nuclei
    :param str agg_type: Like total or mean
    :rtype: pandas.core.frame.DataFrame
    """
    img_file = img_name + SEGMENTATION_X_ANTIGEN_FILE_SUFFIX.replace(
        "SEGMENTATION_TYPE", segmentation_type
    ).replace("AGG_TYPE", agg_type)
    segmentation_quantification_path = input_dir / img_file
    segmentation_quantification_df = read_csv_to_pandas(
        segmentation_quantification_path
    )
    return segmentation_quantification_df
