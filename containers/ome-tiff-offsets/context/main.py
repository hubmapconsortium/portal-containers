import argparse
import os
from glob import glob
from pathlib import Path
from os import makedirs
import json

from skimage.external import tifffile
import xml.dom.minidom

def get_offsets(tiff_filepath):
    with tifffile.TiffFile(tiff_filepath) as tif:
        offsets = [page._offset for page in tif.pages]
    return offsets



def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input_image in glob(input_dir + '/*.ome.tiff') + glob(input_dir + '/*.ome.tif'):
        offset_values = get_offsets(input_image)
        offsets = { "offsetValues": offset_values } 
        with open(Path(output_dir) / str(Path(input_image).with_suffix('').with_suffix('').name + '.json') , 'w') as f:
            f.write(json.dumps(offsets))
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Add offsets to the OMEXML Metadata')
    parser.add_argument(
        '--input_dir', required=True,
        help='Directory containing ome-tiff files to read')
    parser.add_argument(
        '--output_dir', required=True,
        help='Directory where ome-tiff files should be written')
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
