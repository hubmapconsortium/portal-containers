# mudata-to-ui

This container saves [a MuData store](https://mudata.readthedocs.io/en/latest/api/generated/mudata.read_h5mu.html#mudata.read_h5mu) in `zarr` format for viewing in the browser.  

It also selects an appropriate subset of genes to be used for visualization and generates genomic profiles for all the present clusters.

## Input
The input to the container is an [AnnData file in h5ad format](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html).

## Output
The output is the converted zarr store.

## Normalization
All data from the input is scaled to [zero-mean unit-variance] (https://github.com/hubmapconsortium/salmon-rnaseq/blob/master/bin/analysis/scanpy_entry_point.py#L31-L33) `TODO: update line number in the link`.
The `X` is replaced with the log-normalized raw counts to be visualized by Vitessce.

## Example 
Example of a hubmap dataset using this container for data conversion for Vitessce (visualization) would be 
`https://portal.hubmapconsortium.org/browse/dataset/845e7b1c35e8f4926e53b4ef862c0ce7`
