import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

from anndata import read_h5ad

def main(input_dir, output_dir):
    output_dir.mkdir(exist_ok=True)
    for h5ad_file in ['secondary_analysis.h5ad', 'scvelo_annotated.h5ad']
        adata = read_h5ad(h5ad_file)
        # We load 4 channels of Uint16 (i.e 2 bytes per pixel) at around 1200 x 1200 pixels in Viv, so this is a heuristic for what we can show smoothly in Vitessce.
        if len(adata.var) * len(adata.obs) > ((1500 * 1500) / 4) * 4:
            hvg_200 = sc.pp.highly_variable_genes(adata, n_top_genes=200, inplace=False).rename(columns={'means': 'means_hvg_200', 'highly_variable': 'highly_variable_top_200' , 'dispersions': 'dispersions_top_200', 'dispersions_norm': 'dispersions_norm_top_200', 'mean_bin': 'mean_bin_top_200'})
            for col in hvg_200:
                if(col != 'mean_bin_top_200'):
                    adata.var[col] = hvg_200[col].values
            adata.obsm['X_top_200_genes'] = adata[:, adata.var['highly_variable_top_200']].X.copy()
        if 'rank_genes_groups' in adata.uns:
            # Handle marker genes
            for cluster in adata.obs['louvain']:
                adata.obs['marker_genes'][adata.obs['louvain'] == cluster] = adata.uns['rank_genes_groups']['names'][0][cluster]
        adata.write_zarr(output_dir / (Path(h5ad_file).stem + '.zarr'))
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Transform Anndata into zarr.
        ''')
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

