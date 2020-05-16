import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

import pandas as pd
from enum import Enum


class COLUMNS(Enum):
    # Cell Type Annotation column names
    CELL_ID = "cell_id"
    ANNOTATION = "annotation"
    PREDICTION_SCORE = "prediction_score"


def cell_type_annotations_to_factors(input_path, output_path):
    df = pd.read_csv(input_path)
    
    # Get a list containing all unique annotation values in the input file.
    annotation_values = df[COLUMNS.ANNOTATION.value].unique().tolist()

    # Append a temporary column to the dataframe containing the index of each annotation,
    # relative to the `annotation_values` list.
    df["annotation_index"] = df[COLUMNS.ANNOTATION.value].apply(lambda val: annotation_values.index(val))
    
    # Construct the factors dict, which will be saved as a JSON object.
    factors = {
        "Cell Type Annotations": {
            "map": annotation_values,
            "cells": dict(zip(df[COLUMNS.CELL_ID.value].tolist(), df["annotation_index"].tolist()))
        }
    }
    pretty_json_factors = json.dumps(factors).replace('}},', '}},\n')
    with open(output_path, 'w') as f:
        f.write(pretty_json_factors)


def main(input_dir, output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    for input in glob(input_dir + '/*.csv'):
        input_path = Path(input)
        factors_name = input_path.with_suffix('.factors.json').name
        factors_path = Path(output_dir) / factors_name
        cell_type_annotations_to_factors(input_path, factors_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Transform Cell Type Annotation CSV files into Vitessce factors files.
        ''')
    parser.add_argument(
        '--input_dir', required=True,
        help='directory containing .csv files to read')
    parser.add_argument(
        '--output_dir', required=True,
        help='directory where arrow files should be written')
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
