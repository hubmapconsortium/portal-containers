#!/usr/bin/env python

import argparse
from glob import glob
from pathlib import Path

from aicsimageio import AICSImage
from aicsimageio.writers import ome_tiff_writer
from tifffile import TiffFile



def main(expected_dir, actual_dir):
    # potentially we could have .ome.tiff
    for actual_file in glob(actual_dir + '/*.ome.tif*'):
        # get image structured annotations
        with AICSImage(actual_file) as actual_image:
            structured_annotations_actual = actual_image.metadata.structured_annotations
            # create the new output path for the ome tiff
        expected_file = Path(expected_dir) / Path(actual_file).name
        with AICSImage(expected_file) as expected_image:
            structured_annotations_expected = expected_image.metadata.structured_annotations
        offsets_expected = structured_annotations_expected.get_original_metadata_value(key='IFD_Offsets')
        offsets_actual = structured_annotations_actual.get_original_metadata_value(key='IFD_Offsets')
        assert offsets_expected == offsets_actual
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Add offsets to the OMEXML Metadata
        ''')
    parser.add_argument(
        '--expected_dir', required=True,
        help='directory containing expected ome-tiff outputs')
    parser.add_argument(
        '--actual_dir', required=True,
        help='directory containing actual ome-tiff outputs')
    args = parser.parse_args()
    main(args.expected_dir, args.actual_dir)
