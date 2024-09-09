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
    num_repeats = 5

    data = array(
        [[i for i in range(var_dim)] for _ in range(obs_dim)],
        dtype=float32
    )

    obs_metadata = {
        'source file': ['file_{}'.format(i % num_repeats) for i in range(obs_dim)],
        'mask name': ['mask_{}'.format(i % num_repeats) for i in range(obs_dim)],
        'mask id': ['mask_id_{}'.format(i % num_repeats) for i in range(obs_dim)],
        'protocol for mask creation (DOI)': ['doi_{}'.format(i % num_repeats) for i in range(obs_dim)],
        'annotation tool': ['tool_{}'.format(i % num_repeats) for i in range(obs_dim)],
    }

    obsm_data = {
        'X_spatial': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
        'morphology': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
        'ontology': array([[0, 1, 1] for _ in range(obs_dim)], dtype=float32),
    }
    # To avoid anndata throwing warning on numeric index
    obs_df = DataFrame(obs_metadata, index=[str(i) for i in range(obs_dim)])

    adata = AnnData(
        X=data,
        obs=obs_df,
        obsm=obsm_data
    )

    mdata = MuData({'default': adata})
    # to remove default: in the obs column-names
    mdata.obs.columns = mdata.obs.columns.str.replace('default:', '')
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
