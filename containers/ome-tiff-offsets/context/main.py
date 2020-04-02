import argparse
from glob import glob
from pathlib import Path
from os import makedirs

from aicsimageio import AICSImage
from aicsimageio.writers import ome_tiff_writer
from tifffile import TiffFile
import xml.dom.minidom

def get_offsets(tiff_filepath):
    with TiffFile(tiff_filepath) as tif:
        offsets = [page.offset for page in tif.pages]
    return offsets



def main(input_dir, output_dir):
    makedirs(output_dir, exist_ok=True)
    for input in glob(input_dir + '/*.ome.tif*') + glob(input_dir + '/*.ome.tiff'):
        # get image metadata and image data
        with AICSImage(input) as input_image:
            image_data_from_input = input_image.get_image_data()[0]
            omexml = input_image.metadata
        # get offsets and add them to the omexml as structured annotations
        offsets = get_offsets(input)
        structured_annotations = omexml.structured_annotations
        structured_annotations.add_original_metadata(key='IFD_Offsets', value=str(offsets))
        # create the new output path for the ome tiff
        input_path = Path(input)
        new_ome_tiff_path = Path(output_dir) / input_path.name
        # write the file out
        with open(Path(output_dir) / 'ome.xml', 'w') as xml_write:
            xml_write.write(xml.dom.minidom.parseString(str(omexml)).toprettyxml())
        with ome_tiff_writer.OmeTiffWriter(new_ome_tiff_path) as ome_writer:
            ome_writer.save(
                image_data_from_input,
                ome_xml = omexml
            )
        


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
