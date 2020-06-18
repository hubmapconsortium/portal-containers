import argparse
from pathlib import Path

import numpy as np
import pandas as pd

def create_csv(csv_path: Path):
    data = np.array(
        [
            [-1, -1, 0],
            [0, 0, 1],
            [1, 1, 2],
        ],
        dtype=float,
    )

    d = pd.DataFrame(
        data,
        index=['CAT', 'TAG', 'ATG'],
        columns=["umap.1", "umap.2", "cluster"],
    )
    d.cluster = d.cluster.astype(int)

    h5ad = AnnData(
        X=array([[], [], []]),
        obs=DataFrame(
            index=['CAT', 'TAG', 'ATG']
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


def main(output_dir: Path):
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_path = Path(output_dir) / 'fake.csv'
    create_csv(csv_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Creates a minimal CSV input fixture,
            with a structure similar to what we'll actually receive.
        ''')
    parser.add_argument(
        'dest',
        help='Directory where arrow files should be written',
        type=Path
    )
    args = parser.parse_args()

    main(args.dest)
