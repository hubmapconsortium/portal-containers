import argparse
from pathlib import Path
from os import mkdir
import base64
import xml.etree.ElementTree as ET

import numpy as np
from tifffile import imwrite
from ome_types import from_xml, to_dict
import numpy as np
import base64
import io
from PIL import Image

file_name = 'geomx_ome'

def create_mask_bin_data():
    """
    Creates a binary mask in base64-encoded format to be embedded within an OME-XML file.

    The function generates a binary mask of a specified size (113 bytes) using zeroed byte values (`b'\x00'`). 
    It then base64-encodes this binary data and formats it into a BinData element in OME-XML. The mask data 
    is intended to be used for specifying regions of interest (ROI) in OME-TIFF files.

    Returns:
        str: A string containing the BinData XML element, including the base64-encoded mask data and 
             the appropriate 'Length' attribute for the OME-XML structure.
    """
    # the max dimensions in the roi is 30x30=900 bits, which become ~113 bytes (900/8) 
    binary_mask_data =  b'\x00' * 113  
    encoded_bin_data = base64.b64encode(binary_mask_data).decode("utf-8")
    bin_data_element = f'<BinData BigEndian="false" Length="{len(encoded_bin_data)}">{encoded_bin_data}</BinData>'    
    return bin_data_element

def create_ome_tiff_for_geomx(output_path):
    """
    Creates an OME-TIFF file from a defined OME-XML mock sample.

    Returns:
        Saves the file to the specified output_path
    """

    bin_data_element = create_mask_bin_data()
    ome_xml_string  = f"""<?xml version="1.0" encoding="UTF-8"?>
        <OME xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd" Creator="Nanostring GeoMx 2.4.2.2" xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
            <Image ID="Image:0" Name="Image1">
                <Pixels ID="Pixels:0" BigEndian="false" Type="uint16" SignificantBits="16" Interleaved="false" DimensionOrder="XYZCT" PhysicalSizeX="0.400673121" PhysicalSizeY="0.399840057" SizeX="200" SizeY="200" SizeZ="1" SizeC="1" SizeT="1" PhysicalSizeXUnit="µm" PhysicalSizeYUnit="µm">
                    <Channel ID="Channel:0" Name="Channel0" Color="65279" Fluor="SYTO 13" SamplesPerPixel="1"/>
                    <TiffData FirstC="0" FirstT="0" FirstZ="0" IFD="0" />
                </Pixels>
                <ROIRef ID="ROI:0" />
                <ROIRef ID="ROI:1" />
                <AnnotationRef ID="Annotation:3" />
            </Image>
            <StructuredAnnotations>
                <XMLAnnotation ID="Annotation:1">
                    <Value>
                    <ChannelThresholds>
                        <RoiName>001</RoiName>
                        <Bluethreshold>0</Bluethreshold>
                        <Greenthreshold>50</Greenthreshold>
                        <Yellowthreshold>255</Yellowthreshold>
                        <Redthreshold>20</Redthreshold>
                    </ChannelThresholds>
                    </Value>
                </XMLAnnotation>
                <XMLAnnotation ID="Annotation:2">
                    <Value>
                    <ChannelThresholds>
                        <RoiName>002</RoiName>
                        <Bluethreshold>0</Bluethreshold>
                        <Greenthreshold>51</Greenthreshold>
                        <Yellowthreshold>255</Yellowthreshold>
                        <Redthreshold>12</Redthreshold>
                    </ChannelThresholds>
                    </Value>
                </XMLAnnotation>
                <XMLAnnotation ID="Annotation:3">
                    <Value>
                        <SegmentDefinition>
                        <Name>S1</Name>
                        <CollectionOrder>1</CollectionOrder>
                        <DisplayColor>#f28888</DisplayColor>
                        <BlueSelection>3</BlueSelection>
                        <GreenSelection>3</GreenSelection>
                        <YellowSelection>3</YellowSelection>
                        <RedSelection>3</RedSelection>
                        <Erode>2</Erode>
                        <Dilate>2</Dilate>
                        <HoleSize>160</HoleSize>
                        <ParticleSize>25</ParticleSize>
                        </SegmentDefinition>
                    </Value>
                </XMLAnnotation>
            </StructuredAnnotations>
            <ROI ID="ROI:0">
                <Union>
                    <Label ID="Shape:0" Text="001" X="5" Y="5"/>
                    <Polygon ID="Shape:1" Points="10,10 60,10 60,60 10,60 10,10"/>
                    <Mask ID="Shape:2" X="10" Y="10" Width="30" Height="30" Text="S1">
                    {bin_data_element}
                    </Mask>
                </Union>
                <AnnotationRef ID="Annotation:1"/>
            </ROI>
            <ROI ID="ROI:1">
                <Union>
                    <Label ID="Shape:3" Text="002" X="10" Y="12"/>
                    <Polygon ID="Shape:4" Points="20,20 80,20 80,80 20,80 20,20"/>
                    <Mask ID="Shape:5" X="20" Y="20" Width="25" Height="25" Text="S1"> {bin_data_element}
                    </Mask>
                </Union>
                <AnnotationRef ID="Annotation:2"/>
            </ROI>
        </OME>"""

    image_data = np.random.randint(0, 256, size=(200, 200, 4), dtype=np.uint16)
    ome_xml_bytes = ome_xml_string.encode('utf-8')
    imwrite(output_path, image_data, description=ome_xml_bytes, metadata=None)
    print(f"OME-TIFF written to {output_path}")

def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    ome_path = Path(output_dir) / f"{file_name}.ome.tiff"
    create_ome_tiff_for_geomx(ome_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Creates a minimal ome-tiff fixture,
            with a structure similar to what we'll actually receive.
        """
    )
    parser.add_argument(
        "dest", help="Directory where test files should be written")
    args = parser.parse_args()
    main(args.dest)