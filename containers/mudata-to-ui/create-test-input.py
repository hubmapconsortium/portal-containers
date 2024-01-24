import argparse
from pathlib import Path
from os import mkdir
import math

from mudata import MuData
from anndata import AnnData
from numpy import array
from pandas import DataFrame
from scipy.sparse import csr_matrix


def create_h5mu(h5mu_path):
    data = array(
        [[i for i in range(15)], [i for i in range(15)],
         [i for i in range(15)]]
    )
    log_data = array(
        [[math.log(1 + i) for i in range(15)], [math.log(1 + i)
                                                for i in range(15)], [math.log(1 + i) for i in range(15)]]
    )
    layers = {'unscaled': csr_matrix(log_data)}

    obs = DataFrame(index=["TCG", "TAC", "GTC"], data={"leiden": [0, 1, 1]})
    var = DataFrame(
        index=[f"gene_{i}" for i in range(15)],
        data={"dispersions_norm": [i for i in range(15)]},
    )
    obsm = {"X_umap": array([[-1, -1], [0, 0], [1, 1]])}
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
        obs=obs,
        var=var,
        obsm=obsm,
        uns=uns,
        layers=layers,
    )
    adata_atac = AnnData(
        X=data,
        obs=obs,
        var=var,
        obsm=obsm,
        uns=uns,
        layers=layers,
    )
    h5mu = MuData(
        {
            'rna': adata_rna,
            'atac_cbg': adata_atac
        }
    )
    h5mu.obs['leiden_wnn'] = [0, 1, 1]
    h5mu.write(h5mu_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5mu_path = Path(output_dir) / "multiome_downstream_tfidf.h5mu"
    create_h5mu(h5mu_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Creates a minimal mudata input fixture,
            with a structure similar to what we'll actually receive.
        """
    )
    parser.add_argument(
        "dest", help="Directory where zarr files should be written")
    args = parser.parse_args()
    main(args.dest)
