import argparse
from pathlib import Path
from os import mkdir

from anndata import AnnData
from numpy import array
from pandas import DataFrame

def create_h5ad(h5ad_path):
    h5ad = AnnData(
        X=array([[], [], []]),
        obs=DataFrame(
            # data=[1, 2, 3],
            index=['CAT', 'TAG', 'ATG']
            # columns=['leiden']
        ),
        obsm={
            'X_umap': array([
                [-1, -1],
                [0, 0],
                [1, 1]
            ])
        }
    )
    print(h5ad)
    # The real data looks like this:
    #
    # AnnData object with n_obs × n_vars = 5078 × 9308
    # obs: 'n_genes', 'n_counts', 'leiden'
    # var: 'n_cells'
    # uns: 'leiden', 'leiden_colors', 'neighbors', 'rank_genes_groups'
    # obsm: 'X_pca', 'X_umap'
    h5ad.write(h5ad_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5ad_path = Path(output_dir) / 'fake.h5ad'
    create_h5ad(h5ad_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Creates a minimal anndata input fixture,
            with a structure similar to what we'll actually receive.
        ''')
    parser.add_argument(
        'dest',
        help='Directory where arrow files should be written')
    args = parser.parse_args()
    main(args.dest)
