import argparse
from glob import glob
from pathlib import Path
from os import makedirs
from itertools import chain
from ome_types import from_tiff

from create_segments_ome_tiff import create_segments_ome_tiff

def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)

    # Find all OME.TIFFs in the input directory.
    tiffs = list(chain(input_dir.glob('**/*.ome.tif'), input_dir.glob('**/*.ome.tiff')))
    if not tiffs:
        raise Exception(f'No OME TIFFs found in {input_dir}')
    for input_path in tiffs:
        ome = from_tiff(input_path)
        # Create output path for each OME.TIFF:
        new_output_dir = (output_dir / input_path.relative_to(input_dir)).parent
        new_output_dir.mkdir(parents=True, exist_ok=True)
        output_path = (output_dir / input_path.relative_to(input_dir).with_suffix('').with_suffix(''))
        output_path.mkdir(parents=True, exist_ok=True)
        # print(new_output_dir, str(output_path))
        create_segments_ome_tiff(ome, str(output_path))
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Transform GeoMx Ome-tiff to segments")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing GeoMx ome-tiff to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where transformed files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)