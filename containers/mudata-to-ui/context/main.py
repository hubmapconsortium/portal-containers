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
INPUT_FILE_NAME = "multiome_downstream_tfidf.h5mu"


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5mu_file in [INPUT_FILE_NAME]:
        # Check if input file exists, skip it if it doesn't exist
        input_path = path.join(input_dir, h5mu_file)
        if not path.exists(input_path):
            print(f"Input file {h5mu_file} does not exist in input directory.")
            continue
        mdata = read_h5mu(input_dir / h5mu_file)

        rna, atac = mdata.mod['rna'], mdata.mod['atac_cbg']

        # Copy clusterings from other modalities to RNA assuming the indices are exactly aligned
        rna.obs['leiden_from_atac_cbg'] = atac.obs['leiden']
        rna.obs['leiden_wnn_from_root'] = mdata.obs['leiden_wnn']

        # Marker gene logic copied from anndata-to-ui
        if "rank_genes_groups" in rna.uns:
            # Handle marker genes by putting top n per cluster in `obs` for `factors` visualization.
            marker_genes = []
            for i in range(NUM_MARKER_GENES_TO_VISUALIZE):
                rna.obs[f"marker_gene_{str(i)}"] = ["" for v in rna.obs.index]
                for cluster in rna.obs["leiden"]:
                    marker_gene = rna.uns["rank_genes_groups"]["names"][i][cluster]
                    rna.obs.loc[rna.obs["leiden"] == cluster, f"marker_gene_{str(i)}"] = marker_gene
                    marker_genes.append(marker_gene)
            rna.var["marker_genes_for_heatmap"] = [
                gene in marker_genes for gene in mdata.var.index
            ]
        
        # Dispersion normalization logic copied from anndata-to-ui
        if "dispersions_norm" in rna.var:
            top_dispersion = rna.var["dispersions_norm"][
                sorted(
                    range(len(rna.var["dispersions_norm"])),
                    key=lambda k: rna.var["dispersions_norm"][k],
                )[-len(rna.obs['leiden'].unique()) * NUM_MARKER_GENES_TO_VISUALIZE:][0]
            ]
            rna.var["top_highly_variable"] = (
                rna.var["dispersions_norm"] > top_dispersion
            )

        # Convert sparse layer matrices to CSC format for performance
        for modality in [rna, atac]:
            for layer in modality.layers:
                if isinstance(modality.layers[layer], sparse.spmatrix):
                    modality.layers[layer] = modality.layers[layer].tocsc()

        # If the main matrix is sparse, it's best for performance to
        # use non-sparse formats to keep the portal responsive.
        # In the future, we should be able to use CSC sparse data natively
        # and get equal performance:
        # https://github.com/theislab/anndata/issues/524 
        for data_layer in [mdata, rna, atac]:
            if isinstance(data_layer.X, sparse.spmatrix):
                data_layer.X = data_layer.X.todense()

        # It is now possible for adata.X to be empty and have shape (0, 0)
        # so we need to check for that here, otherwise there will
        # be a division by zero error during adata.write_zarr
        # Reference: https://github.com/hubmapconsortium/salmon-rnaseq/blob/dfb0e2a/bin/analysis/scvelo_analysis.py#L69
        chunks = (mdata.shape[0], VAR_CHUNK_SIZE) if mdata.shape[1] >= VAR_CHUNK_SIZE else None
        
        zarr_path = output_dir / (Path(h5mu_file).stem + ".zarr")
        mdata.write_zarr(zarr_path, chunks=chunks)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=f"Transform Mudata into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing MuData .h5mu files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where MuData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
