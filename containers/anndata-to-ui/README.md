# anndata-to-ui

This container saves [an AnnData store](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html) in `zarr` format for viewing in the browser due to it's scalability, performance, and flexibility features.  It also selects an appropriate subset of genes to be used for visualization.

## Input
The input to the container is an [AnnData file in h5ad format](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html).

## Output
The output is the converted `zarr` store.

## Normalization
All data from the input is scaled to [zero-mean unit-variance] (https://github.com/hubmapconsortium/salmon-rnaseq/blob/main/bin/analysis/scanpy_entry_point.py#L31-L33).
The `X` is replaced with the log-normalized raw counts to be visualized by Vitessce.

## Example
Example of a hubmap dataset using this container for data conversion would be 
`https://portal.hubmapconsortium.org/browse/HBM856.HVWM.567`