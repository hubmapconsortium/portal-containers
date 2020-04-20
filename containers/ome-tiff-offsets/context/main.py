import argparse
from glob import glob
from pathlib import Path
from os import makedirs
import json

from skimage.external import tifffile


def get_offsets(tiff_filepath):
    with tifffile.TiffFile(tiff_filepath) as tif:
        offsets = [page._offset for page in tif.pages]
    return offsets

def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input_path in glob(input_dir + '/*.ome.tiff') + glob(input_dir + '/*.ome.tif'):
        offsets = get_offsets(input_path)
        # Truncate `.ome.tiff`.
        input_name = Path(input_path).with_suffix('').with_suffix('').name
        output_path = Path(output_dir) / str(input_name + '.offsets.json')
        with open(output_path, 'w') as f:
            f.write(json.dumps(offsets))
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Create a json file for OME-TIFF offsets')
    parser.add_argument(
        '--input_dir', required=True,
        help='Directory containing ome-tiff files to read')
    parser.add_argument(
        '--output_dir', required=True,
        help='Directory where ome-tiff offsets should be written')
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
