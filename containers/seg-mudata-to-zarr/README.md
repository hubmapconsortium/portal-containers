# seg-mudata-to-zarr

This container saves [a MuData store](https://mudata.readthedocs.io/en/latest/api/generated/mudata.read_h5mu.html#mudata.read_h5mu) in `zarr` format for the segmentation mask for viewing in the browser.
For Vitessce, the index of the zarr store needs to be a number so that it can be mapped to the segmentation mask image for Visualization purposes (e.g., enabiling hovering) and mudata contains a compound index (e.g., tubules-202) due to multiple masks being detected at the same index.  

## Input
The input to the container is an [MuData store](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html).


## Output
The output includes the converted `zarr` store for each of the mask name present in the `mudata` as well as  `metadata.json` file that includes the mask names needed to retrieve the corresponding zarr files for visualization.
 

## Normalization
The compound index from the mudata is transformed into numeric index by extracting the numeric part from it.

## Example
Example of a hubmap dataset using this container for data conversion would be `TODO` 