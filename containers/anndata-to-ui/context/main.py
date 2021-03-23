import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

import zarr
from scipy import sparse
from anndata import read_h5ad
import scanpy as sc
from constants import Assay

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
SECONDARY_ANALYSIS = "secondary_analysis.h5ad"
SCVELO_ANNOTATED = "scvelo_annotated.h5ad"


def main(input_dir, output_dir, assay):
    output_dir.mkdir(exist_ok=True)
    for h5ad_file in ["secondary_analysis.h5ad", "scvelo_annotated.h5ad"]:
        adata = read_h5ad(input_dir / h5ad_file)
        if "rank_genes_groups" in adata.uns:
            # Handle marker genes by putting top n per cluster in `obs` for `factors` visualization.
            marker_genes = []
            for i in range(NUM_MARKER_GENES_TO_VISUALIZE):
                adata.obs[f"marker_gene_{str(i)}"] = ["" for v in adata.obs.index]
                for cluster in adata.obs["leiden"]:
                    marker_gene = adata.uns["rank_genes_groups"]["names"][i][cluster]
                    adata.obs[f"marker_gene_{str(i)}"][
                        adata.obs["leiden"] == cluster
                    ] = marker_gene
                    marker_genes.append(marker_gene)
            adata.var["marker_genes_for_heatmap"] = [
                gene in marker_genes for gene in adata.var.index
            ]
        if "dispersions_norm" in adata.var:
            top_dispersion = adata.var["dispersions_norm"][
                sorted(
                    range(len(adata.var["dispersions_norm"])),
                    key=lambda k: adata.var["dispersions_norm"][k],
                )[-len(adata.obs['leiden'].unique()) * NUM_MARKER_GENES_TO_VISUALIZE:][0]
            ]
            adata.var["top_highly_variable"] = (
                adata.var["dispersions_norm"] > top_dispersion
            )
        for layer in adata.layers:
            if isinstance(adata.X, sparse.spmatrix):
                adata.X = adata.X.tocsc()
    
        # All data from secondary_analysis is scaled at the moment to zero-mean unit-variance
        # https://github.com/hubmapconsortium/salmon-rnaseq/blob/master/bin/analysis/scanpy_entry_point.py#L47
        # We currently cannot visaulize this in Vitessce so we replace `X` with the raw counts.
        # Ideally, we should be able to manage the `layers` and `X` simultaneously in `zarr` but currently we cannot:
        # https://github.com/theislab/anndata/issues/524
        if (SECONDARY_ANALYSIS == h5ad_file):
            # Need .copy() to prevent scanpy from altering the original layer.
            adata.layers['scaled'] = adata.X.copy()
            adata.X = adata.layers[assay.secondary_analysis_layer].copy()
            # From:
            # https://github.com/hubmapconsortium/salmon-rnaseq/blob/4db77ce1de9019ac997dd5da553ea0f11a15d9f9/bin/analysis/scanpy_entry_point.py#L39-L40
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        zarr_path = output_dir / (Path(h5ad_file).stem + ".zarr")

        # If the matrix is sparse, it's best for performance to
        # use non-sparse formats to keep the portal responsive.
        # In the future, we should be able to use CSC sparse data natively
        # and get equal performance:
        # https://github.com/theislab/anndata/issues/524 
        if isinstance(adata.X, sparse.spmatrix):
            adata.X = adata.X.todense()
        adata.write_zarr(zarr_path, [adata.shape[0], VAR_CHUNK_SIZE])


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
    parser.add_argument(
        "--assay",
        required=True,
        # From: https://github.com/hubmapconsortium/salmon-rnaseq/blob/54a5a126ce7d4aefe8677623f4552f34b98ab008/bin/common/common.py
        choices=list(Assay),
        type=Assay,
        help="Assay name",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir, args.assay)
