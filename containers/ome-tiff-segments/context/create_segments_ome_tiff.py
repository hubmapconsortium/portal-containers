import numpy as np
from ome_types import OME, model
import re
from tifffile import imwrite  
#  python context/main.py --input_dir test-input --output_dir test-output-expected

output_file_name= "segmentations.ome.tif"

def decode_binary_mask(bin_data, width, height):
    """
    Decodes binary mask data and reshapes it into a 2D array matching the given dimensions.

    This function takes binary mask data as input, unpacks it into an array of bits, and reshapes 
    it into a 2D NumPy array that corresponds to the specified width and height of the mask. 
    It ensures that the unpacked bit array is resized to match the expected size based on the mask's dimensions.

    Args:
        bin_data (bytes): The binary mask data, typically in byte format, to be decoded.
        width (int): The width of the mask in pixels.
        height (int): The height of the mask in pixels.

    Returns:
        np.ndarray: A 2D NumPy array representing the decoded binary mask, with shape (height, width).
    """

    try:
        image_array = np.frombuffer(bin_data, dtype=np.uint8)
        binarray = np.unpackbits(image_array)
        expected_size = width * height
        imgarr = np.reshape(binarray[:expected_size], (height, width))
        return imgarr
    
    except Exception as e:
        raise ValueError(f'Error in decoding binary mask {str(e)}')

def create_masks_from_rois(ome: OME, size_x: int, size_y: int):
    """
    Places the masks from the ROIs in an OME object into one full-size channel per segment.

    For each mask in the OME object, the function decodes its binary data and places it into the channel
    for its segment (the 'text' field of the mask), with the mask's pixels set to its Shape ID.

    Args:
        ome (OME): The OME object containing ROI and shape data.
        size_x (int): The width of the original image in which to place the masks.
        size_y (int): The height of the original image in which to place the masks.

    Returns:
        tuple: The sorted segment names, and a (segments, size_y, size_x) uint16 array with
            one channel per segment, in that order.
    """
    masks = [shape for roi in ome.rois for shape in roi.union if isinstance(shape, model.Mask)]
    # Same fallback as create_aoi_table. Sorted because portal-visualization labels the channels
    # with the AOI table's segment categories, which anndata stores sorted.
    segments = sorted({mask.text or 'N/A' for mask in masks})
    # Filled in place, so the channels are not copied again to stack them for writing.
    # ponytail: still holds every full-size channel; write page by page if slides outgrow memory.
    channels = np.zeros((len(segments), size_y, size_x), dtype=np.uint16)

    for mask in masks:
        width = int(mask.width)
        height = int(mask.height)
        x_offset = int(mask.x)
        y_offset = int(mask.y)
        mask_array = decode_binary_mask(mask.bin_data.value, width, height)
        # Extract the Shape ID (e.g., Shape:8 -> 8)
        shape_id = int(re.search(r'\d+', mask.id).group())
        max_id = np.iinfo(channels.dtype).max
        if shape_id > max_id:
            raise ValueError(f'Shape ID {shape_id} ({mask.id}) is over the uint16 limit of {max_id}')

        # Clip the mask to the image on both sides, and skip it if none of it is inside.
        x0, y0 = max(x_offset, 0), max(y_offset, 0)
        x1, y1 = min(x_offset + width, size_x), min(y_offset + height, size_y)
        if x0 >= x1 or y0 >= y1:
            print(f'Skipping {mask.id}: entirely outside the {size_x}x{size_y} image')
            continue

        # Assign Shape ID (e.g., AOI ID) to the mask pixels. The mask is uint8, so widen it
        # first: numpy raises OverflowError for a Shape ID over 255.
        region = channels[segments.index(mask.text or 'N/A'), y0:y1, x0:x1]
        np.maximum(
            region,
            mask_array[y0 - y_offset:y1 - y_offset, x0 - x_offset:x1 - x_offset].astype(np.uint16) * shape_id,
            out=region,
        )

    return segments, channels

def save_masks_to_ome_tiff(segments, channels, output_path, img_name):
    try:
        imwrite(output_path, channels,
                metadata={'axes': 'CYX', 'Name': img_name, 'Channel': [{'Name': name} for name in segments]})
    except Exception as e:
        raise ValueError(f'Error in saving the segmentations ome-tiff {str(e)}')


def create_segments_ome_tiff(ome, output_path):
    size_x, size_y = ome.images[0].pixels.size_x, ome.images[0].pixels.size_y
    segments, channels = create_masks_from_rois(ome, size_x, size_y)
    output_ome_tiff = f'{output_path}.{output_file_name}'

    save_masks_to_ome_tiff(segments, channels, output_ome_tiff, img_name="ROI Masks")
    print(f"Masks saved to {output_path}.{output_file_name}")