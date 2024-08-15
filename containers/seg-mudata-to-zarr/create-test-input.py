import argparse
from pathlib import Path
from os import mkdir

from mudata import MuData
from anndata import AnnData
from numpy import array, float32
from pandas import DataFrame


def create_h5mu(h5mu_path):
    obs_dim = 15
    var_dim = 3
    data = array(
        [[i for i in range(var_dim)] for _ in range(obs_dim)],
        dtype=float32
    )

    obs_metadata = {
        'Source file': ['file_{}'.format(i) for i in range(obs_dim)],
        'Mask name': ['mask_{}'.format(i) for i in range(obs_dim)],
        'Mask ID': ['mask_id_{}'.format(i) for i in range(obs_dim)],
        'Protocol for mask creation (DOI)': ['doi_{}'.format(i) for i in range(obs_dim)],
        'Annotation tool': ['tool_{}'.format(i) for i in range(obs_dim)],
    }

    obsm_data = {
        'X_spatial': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
        'morphology': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
        'ontology': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
    }

    adata = AnnData(
        X=data,
        obs=DataFrame(obs_metadata, index=range(obs_dim)),
        obsm=obsm_data
    )

    mdata = MuData({'default': adata})

    # mdata.var_names_make_unique()
    mdata.write(h5mu_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5mu_path = Path(output_dir) / "secondary_analysis.h5mu"
    create_h5mu(h5mu_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Creates a minimal mudata input fixture,
            with a structure similar to what we'll actually receive.
        """
    )
    parser.add_argument(
        "dest", help="Directory where test files should be written")
    args = parser.parse_args()
    main(args.dest)
