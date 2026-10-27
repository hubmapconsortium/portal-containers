"""Plain assert checks for create_masks_from_rois, on a synthetic OME built in memory.

test_docker.sh runs this inside the image: python test_create_segments_ome_tiff.py
"""
import numpy as np
from ome_types import OME, model

from create_segments_ome_tiff import create_masks_from_rois


def ome_with_masks(*masks):
    return OME(rois=[model.ROI(union=list(masks))])


def full_mask(shape_id, segment, x, y, size):
    bits = np.packbits(np.ones(size * size, dtype=np.uint8)).tobytes()
    return model.Mask(id=f'Shape:{shape_id}', text=segment, x=x, y=y, width=size, height=size,
                      bin_data=model.BinData(value=bits, big_endian=False, length=len(bits)))


segments, channels = create_masks_from_rois(ome_with_masks(
    full_mask(300, 'S2', 1, 1, 2),    # over 255: overflowed the uint8 mask before
    full_mask(7, 'S1', -2, -1, 4),    # hangs off the near edges
    full_mask(8, 'S1', 8, 8, 4),      # hangs off the far edges
    full_mask(9, 'S1', 12, 0, 4),     # entirely outside: skipped
), size_x=10, size_y=10)

# Sorted, not in the order the masks appear, so the channels line up with the AOI table.
assert segments == ['S1', 'S2'], segments
expected = np.zeros((2, 10, 10), dtype=np.uint16)
expected[0, 0:3, 0:2] = 7
expected[0, 8:10, 8:10] = 8
expected[1, 1:3, 1:3] = 300
assert channels.dtype == np.uint16, channels.dtype
assert (channels == expected).all(), channels

try:
    create_masks_from_rois(ome_with_masks(full_mask(70000, 'S1', 0, 0, 1)), size_x=10, size_y=10)
except ValueError as e:
    assert 'Shape ID 70000' in str(e), e
else:
    raise AssertionError('Shape ID over 65535 was not rejected')

print('create_masks_from_rois checks passed')
