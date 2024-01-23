import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ, path
import json

import zarr
from scipy import sparse
from anndata import read_h5ad
from mudata import read_h5mu

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
SECONDARY_ANALYSIS = "secondary_analysis.h5ad"
SCVELO_ANNOTATED = "scvelo_annotated.h5ad"
MUDATA = "multiome_downstream_tfidf.h5mu"


# Select marker genes by putting top n per cluster in `obs` for `factors` visualization.
def rank_genes_groups(adata):
    marker_genes = []
    for i in range(NUM_MARKER_GENES_TO_VISUALIZE):
        adata.obs[f"marker_gene_{str(i)}"] = ["" for v in adata.obs.index]
        for cluster in adata.obs["leiden"]:
            marker_gene = adata.uns["rank_genes_groups"]["names"][i][cluster]
            adata.obs.loc[adata.obs["leiden"] == cluster,
                          f"marker_gene_{str(i)}"] = marker_gene
            marker_genes.append(marker_gene)
    adata.var["marker_genes_for_heatmap"] = [
        gene in marker_genes for gene in adata.var.index
    ]


# Identify the top highly variable genes by normalized dispersion values and store them in `var`
# for use in the heatmap visualization.
def calculate_top_dispersion(adata):
    top_dispersion = adata.var["dispersions_norm"][
        sorted(
            range(len(adata.var["dispersions_norm"])),
            key=lambda k: adata.var["dispersions_norm"][k],
        )[-len(adata.obs['leiden'].unique()) * NUM_MARKER_GENES_TO_VISUALIZE:][0]
    ]
    adata.var["top_highly_variable"] = (
        adata.var["dispersions_norm"] > top_dispersion
    )


# All data from secondary_analysis is scaled at the moment to zero-mean unit-variance
# https://github.com/hubmapconsortium/salmon-rnaseq/blob/master/bin/analysis/scanpy_entry_point.py#L47
# We currently cannot visaulize this in Vitessce so we replace `X` with the log-normalized raw counts:
# https://github.com/hubmapconsortium/salmon-rnaseq/commit/9cf1dd4dbe4538b565a0355f56399d3587827eff
# Ideally, we should be able to manage the `layers` and `X` simultaneously in `zarr` but currently we cannot:
# https://github.com/scverse/anndata/issues/524
def rescale(adata):
    adata.layers['scaled'] = adata.X.copy()
    adata.X = adata.layers['unscaled'].copy()


# If the matrix is sparse, it's best for performance to
# use non-sparse formats to keep the portal responsive.
# In the future, we should be able to use CSC sparse data natively
# and get equal performance:
# https://github.com/scverse/anndata/issues/524
def densify_sparse_matrix(adata):
    if isinstance(adata.X, sparse.spmatrix):
        adata.X = adata.X.todense()


# It is now possible for adata.X to be empty and have shape (0, 0)
# so we need to check for that here, otherwise there will
# be a division by zero error during adata.write_zarr
# Reference: https://github.com/hubmapconsortium/salmon-rnaseq/blob/dfb0e2a/bin/analysis/scvelo_analysis.py#L69
def write_zarr(adata, zarr_path):
    chunks = (
        adata.shape[0], VAR_CHUNK_SIZE) if adata.shape[1] >= VAR_CHUNK_SIZE else None

    adata.write_zarr(zarr_path, chunks=chunks)


def file_path_is_valid(input_file_path):
    return path.exists(input_file_path) and path.isfile(input_file_path)


def handle_mudata(input_file_path):
    mdata = read_h5mu(input_file_path)
    rna, atac = mdata.mod['rna'], mdata.mod['atac_cbg']

    # Copy clusterings from other modalities to RNA assuming the indices are exactly aligned
    rna.obs['leiden_from_atac_cbg'] = atac.obs['leiden']
    rna.obs['leiden_wnn_from_root'] = mdata.obs['leiden_wnn']

    # Marker gene logic copied from anndata-to-ui
    if "rank_genes_groups" in rna.uns:
        rank_genes_groups(rna)

    # Dispersion normalization logic copied from anndata-to-ui
    if "dispersions_norm" in rna.var:
        calculate_top_dispersion(rna)

    # Convert sparse layer matrices to CSC format for performance
    for modality in [rna, atac]:
        for layer in modality.layers:
            if isinstance(modality.layers[layer], sparse.spmatrix):
                modality.layers[layer] = modality.layers[layer].tocsc()

    for data_layer in [mdata, rna, atac]:
        densify_sparse_matrix(data_layer)

    return mdata


def handle_anndata(input_file_path, should_rescale=False):
    adata = read_h5ad(input_file_path)
    if "rank_genes_groups" in adata.uns:
        rank_genes_groups(adata)
    if "dispersions_norm" in adata.var:
        calculate_top_dispersion(adata)
    for layer in adata.layers:
        if isinstance(adata.layers[layer], sparse.spmatrix):
            adata.layers[layer] = adata.layers[layer].tocsc()
    densify_sparse_matrix(adata)
    if (should_rescale):
        adata = rescale(adata)
    return adata


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5ad_file in [SECONDARY_ANALYSIS, SCVELO_ANNOTATED, MUDATA]:
        # Check if input file exists, skip it if it doesn't exist
        input_path = path.join(input_dir, h5ad_file)
        if not file_path_is_valid(input_path):
            print(f"Input file {h5ad_file} does not exist.")
            continue

        if (h5ad_file.endswith('.h5mu')):
            result = handle_mudata(input_path)
        else:
            result = handle_anndata(input_path,
                                    should_rescale=(SECONDARY_ANALYSIS == h5ad_file))

        zarr_path = output_dir / (Path(h5ad_file).stem + ".zarr")
        write_zarr(result, zarr_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"Transform Anndata/Mudata into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing AnnData (.h5ad)/MuData (.h5mu) files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
