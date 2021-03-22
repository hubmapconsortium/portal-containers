import argparse
from pathlib import Path
from os import mkdir

from anndata import AnnData
from numpy import array
from pandas import DataFrame
from scipy.sparse import csr_matrix

def create_h5ad_secondary_analysis(h5ad_path):
    data = array(
        [[i for i in range(15)], [i for i in range(15)], [i for i in range(15)]]
    )
    layers={ 'spliced_unspliced_sum': csr_matrix(data) }
    h5ad = AnnData(
        X=data,
        obs=DataFrame(index=["TCG", "TAC", "GTC"], data={"leiden": [0, 1, 1]}),
        var=DataFrame(
            index=[f"gene_{i}" for i in range(15)],
            data={"dispersions_norm": [i for i in range(15)]},
        ),
        obsm={"X_umap": array([[-1, -1], [0, 0], [1, 1]])},
        uns={
            "rank_genes_groups": {
                "names": [
                    ["gene_1", "gene_2"],
                    ["gene_3", "gene_4"],
                    ["gene_4", "gene_5"],
                    ["gene_6", "gene_7"],
                    ["gene_7", "gene_8"],
                ]
            }
        },
        layers=layers,
    )
    print(h5ad)
    h5ad.write(h5ad_path)


def create_h5ad_scvelo(h5ad_path):
    data = array(
        [[i for i in range(15)], [i for i in range(15)], [i for i in range(15)]]
    )
    layers={ 'spliced_unspliced_sum': csr_matrix(data) }
    h5ad = AnnData(
        X=data,
        obs=DataFrame(index=["CTG", "GCA", "CTG"], data={"leiden": [1, 1, 2]}),
        var=DataFrame(index=[f"gene_{i}" for i in range(15)]),
        obsm={"X_umap": array([[-1, -1], [0, 0], [1, 1]])},
        layers=layers
    )
    print(h5ad)
    h5ad.write(h5ad_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5ad_path = Path(output_dir) / "secondary_analysis.h5ad"
    create_h5ad_secondary_analysis(h5ad_path)
    h5ad_path = Path(output_dir) / "scvelo_annotated.h5ad"
    create_h5ad_scvelo(h5ad_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"""
            Creates a minimal anndata input fixture,
            with a structure similar to what we'll actually receive.
        """
    )
    parser.add_argument("dest", help="Directory where arrow files should be written")
    args = parser.parse_args()
    main(args.dest)
