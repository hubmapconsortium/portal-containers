import argparse
import os
from glob import glob
from pathlib import Path
from os import makedirs

from aicsimageio import AICSImage
from aicsimageio.writers import ome_tiff_writer
from skimage.external import tifffile
import xml.dom.minidom

def get_offsets(tiff_filepath):
    with tifffile.TiffFile(tiff_filepath) as tif:
        offsets = [page._offset for page in tif.pages]
    return offsets



def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input in glob(input_dir + '/*.ome.tif*') + glob(input_dir + '/*.ome.tiff'):
        # Get image metadata and image data.
        input_image = AICSImage(input)            
        omexml = input_image.metadata
        offsets = get_offsets(str(input))
        structured_annotations = omexml.structured_annotations
        structured_annotations.add_original_metadata(key='IFD_Offsets', value=str(offsets))
        # Write the file out to the bound output directory.
        output_path = Path(output_dir)
        with open(output_path / str(os.path.splitext(Path(input).name)[0] + '.xml'), 'w') as xml_write:
            xml_write.write(xml.dom.minidom.parseString(str(omexml)).toprettyxml())
        


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
