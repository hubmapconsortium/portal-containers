import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

import zarr
from anndata import read_h5ad


def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5ad_file in ["secondary_analysis.h5ad", "scvelo_annotated.h5ad"]:
        adata = read_h5ad(input_dir / h5ad_file)
        if "rank_genes_groups" in adata.uns:
            # Handle marker genes by putting top 5 per cluster in `obs` for `factors` visualization.
            marker_genes = []
            for i in range(5):
                for cluster in adata.obs["leiden"]:
                    marker_gene = adata.uns["rank_genes_groups"]["names"][i][cluster]
                    adata.obs[f"marker_gene_{str(i)}"] = ["" for v in adata.obs.index]
                    adata.obs[f"marker_gene_{str(i)}"][
                        adata.obs["leiden"] == cluster
                    ] = marker_gene
                    marker_genes.append(marker_gene)
            adata.var["marker_genes_for_heatmap"] = [
                gene in marker_genes for gene in adata.var.index
            ]
        # Chunk the anndata object to be available efficiently for use.
        adata.write_zarr(
            output_dir / (Path(h5ad_file).stem + ".zarr"), chunks=[adata.shape[0], 10]
        )


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
