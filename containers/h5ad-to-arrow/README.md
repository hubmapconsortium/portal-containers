# h5ad-to-arrow

This container translates [anndata's h5ad](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html) to [Apache Arrow](https://arrow.apache.org/),
as well as CSV, and Vitessce JSON which conforms to our [schemas](https://github.com/hubmapconsortium/vitessce/tree/master/src/schemas).
The arrow format is a columnar format optimized for analytical workloads like querying and aggregations and is faster than AnnData's row-based storage for certain operations.

## Input
The input to the container is an [an AnnData file in h5ad format](https://anndata.readthedocs.io/en/latest/anndata.read_h5ad.html).


## Output
The output includes the converted `arrow` file, a csv file representing the arrow file for readability purposes, and json files representing cells and cell sets crucial for Vitessce visualization. 


## Normalization
None

## Example 
Example of a hubmap dataset using this container for data conversion for Vitessce (visualization) would be `HBM768.NCSB.762`
