import argparse
from pathlib import Path
from os import mkdir

from anndata import AnnData
from numpy import array
from pandas import DataFrame

def create_h5ad(h5ad_path):
    X = DataFrame(data=[
        {"ENSG00000282431.1": 0.0, "ENSG00000228985.1": 1.0, "ENSG00000237235.2": 10.0, "ENSG00000223997.1": 2.5},
        {"ENSG00000282431.1": 4.0, "ENSG00000228985.1": 3.0, "ENSG00000237235.2": 3.0, "ENSG00000223997.1": 3.5},
        {"ENSG00000282431.1": 5.5, "ENSG00000228985.1": 8.0, "ENSG00000237235.2": 9.0, "ENSG00000223997.1": 4.0},
    ], index=["CAT", "TAG", "ATG"], columns=["ENSG00000282431.1", "ENSG00000228985.1", "ENSG00000237235.2", "ENSG00000223997.1"])
    var = DataFrame(index=["ENSG00000282431.1", "ENSG00000228985.1", "ENSG00000237235.2", "ENSG00000223997.1"])
    var.index = var.index.rename("index")
    obs = DataFrame(index=['CAT', 'TAG', 'ATG'], columns=['leiden'], data=array([[0], [1], [2]]))
    obs.index = obs.index.rename("index")
    h5ad = AnnData(
        X=X,
        obs=obs,
        obsm={
            'X_umap': array([
                [-1, -1],
                [0, 0],
                [1, 1]
            ])
        },
        var=var,
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
