import argparse
from glob import glob
from pathlib import Path
from os import makedirs
from itertools import chain
import json
from ome_types import from_tiff


def get_metadata(tiff_file_path):
    extracted_metadata = {}
    try:
        ome_data = from_tiff(tiff_file_path)
        image = ome_data.images[0]
        pixels = image.pixels

        physical_size_x = pixels.physical_size_x
        physical_size_y = pixels.physical_size_y

        physical_size_unit_x = pixels.physical_size_x_unit  
        physical_size_unit_y = pixels.physical_size_y_unit  
        

        physical_size_unit_x = physical_size_unit_x.value if hasattr(physical_size_unit_x, 'value') else physical_size_unit_x
        physical_size_unit_y = physical_size_unit_y.value if hasattr(physical_size_unit_y, 'value') else physical_size_unit_y
        extracted_metadata = {
            "PhysicalSizeX": physical_size_x,
            "PhysicalSizeY": physical_size_y,
            "PhysicalSizeUnitX": physical_size_unit_x, 
            "PhysicalSizeUnitY": physical_size_unit_y, 
        }
        print(f"Extracted metadata: {extracted_metadata}")
    except FileNotFoundError:
        print(f"Error: The file {tiff_file_path} does not exist.")
    except IndexError:
        print("Error: The TIFF file does not contain any images.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return extracted_metadata


def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)

    # Find all OME.TIFFs in the input directory.
    tiffs = list(chain(input_dir.glob('**/*.ome.tif'), input_dir.glob('**/*.ome.tiff')))
    if not tiffs:
        raise Exception(f'No OME TIFFs found in {input_dir}')
    for input_path in tiffs:
        metadata = get_metadata(str(input_path))

        # Create output path for each OME.TIFF:
        new_output_dir = (output_dir / input_path.relative_to(input_dir)).parent
        new_output_dir.mkdir(parents=True, exist_ok=True)

        # Set output filename for JSON file and dump to disk:
        output_path = str(output_dir / input_path.relative_to(input_dir).with_suffix('').with_suffix(''))+'.metadata.json'
        with open(output_path, 'w') as f:
            f.write(json.dumps(metadata))

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
