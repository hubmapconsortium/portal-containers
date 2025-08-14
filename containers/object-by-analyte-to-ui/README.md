# mudata-to-ui

This container saves [a MuData store](https://mudata.readthedocs.io/en/latest/api/generated/mudata.read_h5mu.html#mudata.read_h5mu) in `zarr` format for viewing in the browser.

It also selects an appropriate subset of genes to be used for visualization and generates genomic profiles for all the present clusters.

## Input

The input to the container is a [MuData file in h5mu format](https://muon.readthedocs.io/en/latest/io/output.html#id2).

## Output

The output is the converted zarr store in zip-zarr format.

## Normalization

All data from the input is scaled to [zero-mean unit-variance] (https://github.com/hubmapconsortium/multiome-rna-atac-pipeline/blob/a52b6bb37f56dcd78d45ceef1868095d59ef1aac/bin/downstream.py#L30-L37).
