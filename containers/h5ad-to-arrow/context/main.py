import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

from anndata import read_h5ad
import pyarrow as pa
from pandas import DataFrame


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
            cell_sets_json=arrow_path.with_suffix('.cell-sets.json')
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
