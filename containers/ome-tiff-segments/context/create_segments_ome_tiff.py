import numpy as np
from ome_types import OME, model
import re
from tifffile import imwrite  
#  python context/main.py --input_dir test-input --output_dir test-output-expected

output_file_name= "segmentations.ome.tiff"

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
    Generates a dictionary of masks from the ROIs in an OME object, placing them into a original image of the specified size.

    For each ROI in the OME object, the function extracts binary mask data, decodes it, and places the mask
    into a larger image (channel) of size (size_y, size_x). Each segment (identified by the 'text' field of a mask)
    is treated as a separate channel.

    Args:
        ome (OME): The OME object containing ROI and shape data.
        size_x (int): The width of the original image in which to place the masks.
        size_y (int): The height of the original image in which to place the masks.

    Returns:
        dict: A dictionary where keys are segment names and values are 2D NumPy arrays representing 
            the masks placed on a original image.
    """
    channel_images = {}

    for roi in ome.rois:
        try:
            for shape in roi.union:
                if isinstance(shape, model.Mask):
                    try:
                        segment = shape.text  # Segment identifier for the mask
                        mask_data = shape.bin_data.value  # Binary mask data
                        width = int(shape.width)
                        height = int(shape.height)
                        x_offset = int(shape.x)
                        y_offset = int(shape.y)
                    except (AttributeError, ValueError) as e:
                        raise ValueError(f"Error extracting mask attributes from shape: {e}")

                    # Decode the binary mask data
                    try:
                        mask_array = decode_binary_mask(mask_data, width, height)
                    except Exception as e:
                        raise ValueError(f"Failed to decode mask data for segment '{segment}': {e}")
                    
                    # Extract the Shape ID (e.g., Shape:8 -> 8)
                    try:
                        shape_id = int(re.search(r'\d+', shape.id).group())
                    except AttributeError as e:
                        raise ValueError(f"Failed to extract Shape ID from shape: {shape.id}, error: {e}")

                    # Initialize channel if it doesn't exist
                    if segment not in channel_images:
                        channel_images[segment] = np.zeros((size_y, size_x), dtype=np.uint16)

                    # Ensure mask placement within bounds of the larger image
                    x_end = min(x_offset + width, size_x)
                    y_end = min(y_offset + height, size_y)

                    # Assign Shape ID (e.g., AOI ID) to the mask pixels
                    try:
                        channel_images[segment][y_offset:y_end, x_offset:x_end] = np.maximum(
                            channel_images[segment][y_offset:y_end, x_offset:x_end],
                            mask_array[:y_end - y_offset, :x_end - x_offset] * shape_id
                        )
                    except Exception as e:
                        raise ValueError(f"Error placing mask for segment '{segment}' on canvas: {e}")
                    
        except Exception as e:
            print(f"Error processing ROI {roi.id}: {e}")
            continue

    return channel_images

def save_masks_to_ome_tiff(mask_dict, output_path, img_name, size_x, size_y):
    channel_names = list(mask_dict.keys())
    try:
        masks_stack = np.stack([mask_dict[segment] for segment in channel_names], axis=0)
        
        imwrite(output_path, masks_stack, 
                shape=(len(channel_names), size_y, size_x), 
                metadata={'axes': 'CYX', 'Name': img_name, 'Channel': [{'Name': name} for name in channel_names]})
    except Exception as e:
        raise ValueError(f'Error in saving the segmentations ome-tiff {str(e)}')
        

def create_segments_ome_tiff(ome, output_path):
    size_x, size_y = ome.images[0].pixels.size_x, ome.images[0].pixels.size_y
    mask_dict = create_masks_from_rois(ome, size_x, size_y)
    output_ome_tiff = f'{output_path}/{output_file_name}'

    save_masks_to_ome_tiff(mask_dict, output_ome_tiff, img_name="ROI Masks", size_x=size_x, size_y=size_y)
    print(f"Masks saved to {output_ome_tiff}")