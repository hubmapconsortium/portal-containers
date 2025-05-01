import argparse
from pathlib import Path
from os import path
import warnings

import zarr
from scipy import sparse
from numpy import asarray
from anndata import read_h5ad

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
SECONDARY_ANALYSIS = "secondary_analysis.h5ad"
SCVELO_ANNOTATED = "scvelo_annotated.h5ad"


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5ad_file in ["secondary_analysis.h5ad", "scvelo_annotated.h5ad"]:
        print(f"Processing {h5ad_file}")
        # Check if input file exists, skip it if it doesn't exist
        input_path = path.join(input_dir, h5ad_file)
        if not path.exists(input_path):
            print(f"Input file {h5ad_file} does not exist.")
            continue
        adata = read_h5ad(input_dir / h5ad_file)
        if "rank_genes_groups" in adata.uns:
            # Handle marker genes by putting top n per cluster in `obs` for `factors` visualization.
            marker_genes = []
            for i in range(NUM_MARKER_GENES_TO_VISUALIZE):
                adata.obs[f"marker_gene_{str(i)}"] = ["" for v in adata.obs.index]
                for cluster in adata.obs["leiden"]:
                    marker_gene = adata.uns["rank_genes_groups"]["names"][i][cluster]
                    adata.obs.loc[adata.obs["leiden"] == cluster, f"marker_gene_{str(i)}"] = marker_gene
                    marker_genes.append(marker_gene)
            adata.var["marker_genes_for_heatmap"] = [
                gene in marker_genes for gene in adata.var.index
            ]
        if "dispersions_norm" in adata.var:
            top_dispersion = adata.var["dispersions_norm"].iloc[
                sorted(
                    range(len(adata.var["dispersions_norm"])),
                    key=lambda k: adata.var["dispersions_norm"].iloc[k],
                )[-len(adata.obs['leiden'].unique()) * NUM_MARKER_GENES_TO_VISUALIZE:][0]
            ]
            adata.var["top_highly_variable"] = (
                adata.var["dispersions_norm"] > top_dispersion
            )
        for layer in adata.layers:
            if isinstance(adata.layers[layer], sparse.spmatrix):
                adata.layers[layer] = adata.layers[layer].tocsc()
    
        # All data from secondary_analysis is scaled at the moment to zero-mean unit-variance
        # https://github.com/hubmapconsortium/salmon-rnaseq/blob/master/bin/analysis/scanpy_entry_point.py#L31-L33
        # We currently cannot visualize this in Vitessce so we replace `X` with the log-normalized raw counts:
        # https://github.com/hubmapconsortium/salmon-rnaseq/commit/9cf1dd4dbe4538b565a0355f56399d3587827eff
        # Ideally, we should be able to manage the `layers` and `X` simultaneously in `zarr` but currently we cannot:
        # https://github.com/theislab/anndata/issues/524
        if (SECONDARY_ANALYSIS == h5ad_file):
            adata.layers['scaled'] = adata.X.copy()
            adata.X = adata.layers['unscaled'].copy()

        # If the matrix is sparse, it's best for performance to
        # use non-sparse formats to keep the portal responsive.
        # In the future, we should be able to use CSC sparse data natively
        # and get equal performance:
        # https://github.com/theislab/anndata/issues/524 
        if isinstance(adata.X, sparse.spmatrix):
            adata.X = asarray(adata.X.todense())
        
        # It is now possible for adata.X to be empty and have shape (0, 0)
        # so we need to check for that here, otherwise there will
        # be a division by zero error during adata.write_zarr
        # Reference: https://github.com/hubmapconsortium/salmon-rnaseq/blob/dfb0e2a/bin/analysis/scvelo_analysis.py#L69
        chunks = (adata.shape[0], VAR_CHUNK_SIZE) if adata.shape[1] >= VAR_CHUNK_SIZE else None
        zip_path = output_dir / (Path(h5ad_file).stem + ".zarr.zip")

        with zarr.ZipStore(str(zip_path), mode='w') as store:
            with warnings.catch_warnings():
                # To suppress the duplicate warning https://github.com/zarr-developers/zarr-python/issues/129
                warnings.filterwarnings("ignore", category=UserWarning)
                adata.write_zarr(store, chunks=chunks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Transform Anndata into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing AnnData .h5ad files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where AnnData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
