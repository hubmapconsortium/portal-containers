import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json
import re

from anndata import read_h5ad
import pyarrow as pa
from pandas import DataFrame
import zarr
from numcodecs import Zlib
import scipy.cluster
from pyensembl import Genome as EnsemblGenome


def h5ad_to_arrow(h5ad_file, arrow_file):
    ann_data = read_h5ad(h5ad_file)
    umap = ann_data.obsm['X_umap'].transpose()
    leiden = ann_data.obs['leiden'].to_numpy().astype('uint8')
    index = ann_data.obs.index

    df = DataFrame(
        data={'umap_x': umap[0], 'umap_y': umap[1], 'leiden': leiden},
        index=index
    )
    table = pa.Table.from_pandas(df)

    writer = pa.RecordBatchFileWriter(arrow_file, table.schema)
    writer.write(table)
    writer.close()


def arrow_to_csv(arrow_file, csv_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df.to_csv(csv_file)


def arrow_to_json(arrow_file, **kwargs):
    umap_json  = kwargs['umap_json']
    leiden_json  = kwargs['leiden_json']
    cell_sets_json = kwargs['cell_sets_json']
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df_items = df.T.to_dict().items()

    id_to_umap = {
        k: {
            "mappings": {"UMAP": [v['umap_x'], v['umap_y']]},
            "factors": {"Leiden Cluster": str(int(v['leiden']))}
        }
        for (k,v) in df_items
    }
    pretty_json_umap = json.dumps(id_to_umap).replace('}},', '}},\n')
    with open(umap_json, 'w') as f:
        f.write(pretty_json_umap)

    leiden_clusters = sorted(df['leiden'].unique().astype('uint8'))
    id_to_factors = {
        'Leiden Cluster': {
            'map': [str(cluster) for cluster in leiden_clusters],
            'cells': { k: v['leiden'] for (k,v) in df_items }
        }
    }
    pretty_json_factors = json.dumps(id_to_factors).replace('}},', '}},\n')
    with open(leiden_json, 'w') as f:
        f.write(pretty_json_factors)
    
    # Construct the tree, according to the following schema:
    # https://github.com/hubmapconsortium/vitessce/blob/d5f63aa1d08aa61f6b20f6ad6bbfba5fceb6b5ef/src/schemas/cell_sets.schema.json
    cell_sets = {
        "datatype": "cell",
        "version": "0.1.2",
        "tree": [
            {
                "name": "Leiden Cluster",
                "children": [
                    {
                        "name": f"Cluster {cluster}",
                        "set": df.loc[df['leiden'] == cluster].index.values.tolist()
                    }
                    for cluster in leiden_clusters
                ]
            }
        ]
    }
    with open(cell_sets_json, 'w') as f:
        json.dump(cell_sets, f, indent=1)

def ensembl_gene_ids_to_gene_names(gene_ids):
    genome = EnsemblGenome(
        reference_name='GRCh38',
        annotation_name='ensembl',
        gtf_path_or_url=environ['ENSEMBL_URL']
    )

    ensembl_id_regex = re.compile(r"(?P<gene_id>ENSG\d+)(\.\d+)?")

    def gene_id_to_name(gene_id_with_version):
        try:
            match = re.match(ensembl_id_regex, gene_id_with_version)
            if match:
                gene_id = match.group('gene_id')
                return genome.gene_by_id(gene_id).gene_name
        except (ValueError, IndexError):
            pass
        return None

    return [ gene_id_to_name(gene_id) for gene_id in gene_ids ]


def h5ad_to_zarr(h5ad_file, **kwargs):
    expression_matrix_zarr = kwargs['expression_matrix_zarr']

    gexp = read_h5ad(h5ad_file)
    gexp_arr = gexp.X
    gexp_df = gexp.to_df()

    # Re-scale the gene expression values between 0 and 255 (one byte ints).
    gexp_arr_min = gexp_arr.min()
    gexp_arr_max = gexp_arr.max()
    gexp_arr_range = gexp_arr_max - gexp_arr_min
    gexp_arr_ratio = 255 / gexp_arr_range
    gexp_norm_arr = (gexp_arr - gexp_arr_min) * gexp_arr_ratio

    # Perform hierarchical clustering along the genes axis.
    Z = scipy.cluster.hierarchy.linkage(gexp_norm_arr.T, method="ward")
    labels = gexp.var.index.values

    # Get the hierarchy-based ordering of genes.
    leaf_index_list = scipy.cluster.hierarchy.leaves_list(Z)
    leaf_list = labels[leaf_index_list].tolist()

    # Create a new *ordered* gene expression dataframe.
    gexp_norm_df = DataFrame(
        index=gexp_df.index.values.tolist(),
        columns=gexp_df.columns.values.tolist(),
        data=gexp_norm_arr
    )
    sorted_gexp_norm_df = gexp_norm_df[leaf_list]
    sorted_gene_ids = sorted_gexp_norm_df.columns.values.tolist()
    sorted_cell_ids = sorted_gexp_norm_df.index.values.tolist()

    # Convert gene ENSEMBL IDs to gene names/symbols.
    sorted_gene_names = ensembl_gene_ids_to_gene_names(sorted_gene_ids)

    # Save the data to the output file.
    z = zarr.open(
        str(expression_matrix_zarr),
        mode='w',
        shape=sorted_gexp_norm_df.shape,
        dtype='uint8',
        compressor=Zlib(level=1)
    )
    # Store the matrix.
    z[:] = sorted_gexp_norm_df.values
    # Store the rows/observations (cell IDs).
    z.attrs["rows"] = sorted_cell_ids
    # Store the columns/variables (gene IDs).
    z.attrs["cols"] = sorted_gene_ids
    z.attrs["colsAlt"] = sorted_gene_names


def main(input_dir, output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    for input in glob(input_dir + '/*.h5ad'):
        input_path = Path(input)
        arrow_name = input_path.with_suffix('.arrow').name
        arrow_path = Path(output_dir) / arrow_name
        h5ad_to_arrow(input_path, arrow_path)
        arrow_to_csv(arrow_path, arrow_path.with_suffix('.csv'))
        arrow_to_json(
            arrow_file=arrow_path,
            umap_json=arrow_path.with_suffix('.cells.json'),
            leiden_json=arrow_path.with_suffix('.factors.json'),
            cell_sets_json=arrow_path.with_suffix('.cell-sets.json'),
        )
        h5ad_to_zarr(
            h5ad_file=input_path,
            expression_matrix_zarr=arrow_path.with_suffix('.expression-matrix.zarr'),
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Transform H5AD into Arrow, JSON, and CSV.
        ''')
    parser.add_argument(
        '--input_dir', required=True,
        help='directory containing h5ad files to read')
    parser.add_argument(
        '--output_dir', required=True,
        help='directory where arrow files should be written')
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
