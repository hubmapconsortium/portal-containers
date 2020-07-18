import argparse
from glob import glob
from pathlib import Path
from os import makedirs
from itertools import chain
import json
from os.path import dirname, exists, splitext

from skimage.external import tifffile

def get_offsets(tiff_filepath):
    with tifffile.TiffFile(tiff_filepath) as tif:
        offsets = [page._offset for page in tif.pages]
    return offsets

def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)

    # Find all OME.TIFFs in the input directory.
    for input_path in chain(input_dir.glob('**/*.ome.tif'), input_dir.glob('**/*.ome.tiff')):
        offsets = get_offsets(str(input_path))

        # Set output filename for JSON file and dump to disk:
        output_path = input_path.stem.replace('ome','offsets.json')
        with open(output_path, 'w') as f:
            f.write(json.dumps(offsets))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create a json file for OME-TIFF offsets')
    parser.add_argument(
        '--input_dir', required=True, type=Path,
        help='Directory containing ome-tiff files to read')
    parser.add_argument(
        '--output_dir', required=True, type=Path,
        help='Directory where ome-tiff offsets should be written')
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
