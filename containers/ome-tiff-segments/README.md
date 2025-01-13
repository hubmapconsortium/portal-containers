# ome-tiff-segments
This container creates different component datasets needed for the visualization of GeoMx Assays Ome-tiff files which includes regions of interests (ROIs) and Areas of Interests (AOI). These datasets are created by converting the ome-tiff to an ome-xml file.

## Input
The input to the container is the ome-tiff file from GeoMx assay.

## Output
The following output files are generated to support the visualization of GeoMx assays.

- ROIs as `obsSegmentations.json` file by extracting the vertices from Polygon tags within each ROI.

- A Segmentation OME-TIFF file extracted from the Bitmask within each ROI’s mask and grouped by the segment (text field within the mask). 

- AOIs as Zarr store with obs representing the segment, roi-id, and aoi-id. The aoi-id has composite values (e.g., Shape:2), so the index is the numeric part extracted from this composite value.

- ROIs as Zarr store with obs having channel thresholds for the ROIs extracted from the annotations in ome-xml.
 
## Normalization
None

## Example
Example of a hubmap dataset using this container for data conversion would be `TODO` 