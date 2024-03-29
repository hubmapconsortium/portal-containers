import argparse
from glob import glob
from pathlib import Path
from os import makedirs
from itertools import chain
import json
from os.path import dirname, exists, splitext

from tifffile import TiffFile

def get_offsets(tiff_filepath):
    with TiffFile(tiff_filepath) as tif:
        offsets = [page.offset for page in tif.pages]
    return offsets

def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)

    # Find all OME.TIFFs in the input directory.
    tiffs = list(chain(input_dir.glob('**/*.ome.tif'), input_dir.glob('**/*.ome.tiff')))
    if not tiffs:
        raise Exception(f'No OME TIFFs found in {input_dir}')
    for input_path in tiffs:
        offsets = get_offsets(str(input_path))

        # Create output path for each OME.TIFF:
        new_output_dir = (output_dir / input_path.relative_to(input_dir)).parent
        new_output_dir.mkdir(parents=True, exist_ok=True)

        # Set output filename for JSON file and dump to disk:
        output_path = str(output_dir / input_path.relative_to(input_dir).with_suffix('').with_suffix(''))+'.offsets.json'
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
