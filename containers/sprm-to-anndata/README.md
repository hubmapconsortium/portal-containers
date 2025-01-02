# sprm-to-anndata

This container takes [SPRM](https://docs.google.com/document/d/1c7UR0Pe1newpVhQY2HEFkfV8O7GAj9Vk4XnuSiSnDeY/edit) output CSVs and converts them into an [AnnData zarr store](https://anndata.readthedocs.io/en/latest/anndata.AnnData.write_zarr.html#anndata.AnnData.write_zarr).


## Input
The input to the container SPRM output csvs.


## Output
The output of the container is the zarr store created from all the input csv files.
 

## Normalization
None

## Example
Example of a hubmap dataset using this container for data conversion would be `HBM522.BSZT.385` 