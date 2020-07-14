import argparse
from glob import glob
from pathlib import Path
from os import makedirs
from os.path import splitext
import json
from os.path import dirname
from os.path import exists
from os import makedirs

from skimage.external import tifffile

def get_offsets(tiff_filepath):
    with tifffile.TiffFile(tiff_filepath) as tif:
        offsets = [page._offset for page in tif.pages]
    return offsets

def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)

    # Find all OME.TIFFs in the input directory.
    for input_path in glob(input_dir + '/**/*ome.tif', recursive=True) + glob(input_dir + '/**/*ome.tiff', recursive=True):
        offsets = get_offsets(input_path)

        # Create output path for each OME.TIFF:
        output_path = dirname(str(Path(output_dir) / Path(input_path).relative_to(input_dir)))
        if not exists(output_path):
            makedirs(output_path)

        # Set output filename for JSON file and dump to disk:
        output_path = str( Path(output_dir) / Path( input_path ).relative_to( input_dir ).with_suffix('').with_suffix('') ) + '.offsets.json'
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
