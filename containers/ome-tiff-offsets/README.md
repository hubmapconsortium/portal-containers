# ome-tiff-offsets

This docker container creates a JSON list of byte offsets for each TIFF from an input directory. This is needed for visualizing image datasets as it makes visualization much more efficient by allowing requesting specific IFDs and their tiles more efficiently.

## Input
The input to the container is one or more ome-tiff image files.

## Output
The output is a json file that includes an array/list of byte offsets for every input ome-tiff image.

## Normalization
None

## Example 
Example of a hubmap dataset using this container for offsets generation for Vitessce (visualization) would be
`https://portal.hubmapconsortium.org/browse/dataset/1c7e9d1b6de5263aad5e4b95096a6ec5`.

