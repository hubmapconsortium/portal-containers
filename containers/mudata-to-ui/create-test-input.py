import argparse
from pathlib import Path
from os import mkdir
import math

from mudata import MuData
from anndata import AnnData
from numpy import array, float32
from pandas import DataFrame
from scipy.sparse import csr_matrix


def create_h5mu(h5mu_path):
    data = array(
        [[i for i in range(15)],
         [i for i in range(15)],
         [i for i in range(15)]], dtype=float32)
    log_data = array(
        [[math.log(1 + i) for i in range(15)],
         [math.log(1 + i) for i in range(15)],
         [math.log(1 + i) for i in range(15)]]
    )

    index = ["TCG", "TAC", "GTC"]
    obs_data = [0, 1, 1]

    rna_var = DataFrame(
        index=[f"gene_{i}" for i in range(15)],
        data={"dispersions_norm": [i for i in range(15)]}
    )
    cbg_var = DataFrame(
        index=[f"gene_atac_{i}" for i in range(15)],
        data={"highly_variable": [i % 2 == 0 for i in range(15)]},
    )
    uns = {
        "rank_genes_groups": {
            "names": [
                ["gene_1", "gene_2"],
                ["gene_3", "gene_4"],
                ["gene_4", "gene_5"],
                ["gene_6", "gene_7"],
                ["gene_7", "gene_8"],
            ]
        }
    }
    adata_rna = AnnData(
        X=data,
        obs=DataFrame(
            index=index,
            data={"leiden": obs_data}
        ),
        var=rna_var,
        uns=uns,
        layers={
            "spliced": csr_matrix(log_data),
            "spliced_unspliced_sum": csr_matrix(log_data),
            "unspliced": csr_matrix(log_data),
        },
    )
    adata_atac = AnnData(
        X=data,
        obs=DataFrame(
            index=index,
            data={"Clusters": obs_data}
        ),
        var=cbg_var,
        uns=uns,
        layers={
            "smoothed": csr_matrix(log_data),
        },
    )
    h5mu = MuData(
        {
            'rna': adata_rna,
            'atac_cbg': adata_atac,
        },
    )
    h5mu.obs['leiden_wnn'] = obs_data
    h5mu.update()
    h5mu.var_names_make_unique()
    h5mu.write(h5mu_path)


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
