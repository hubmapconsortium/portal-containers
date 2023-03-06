import argparse
from pathlib import Path
import json

import pandas as pd
import pyarrow as pa


def csv_to_arrow(csv_file: Path, arrow_file: Path):
    df = pd.read_csv(csv_file, index_col=0)
    df.index.name = 'index'
    df.columns = ['umap_x', 'umap_y', 'leiden']

    table = pa.Table.from_pandas(df)

    writer = pa.RecordBatchFileWriter(arrow_file, table.schema)
    writer.write(table)
    writer.close()


def arrow_to_csv(arrow_file, csv_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df.to_csv(csv_file)

def try_str_to_int(val):
    try:
        # Keep only numeric characters.
        val_numeric = "".join(filter(str.isdigit, val))
        if len(val_numeric) > 0:
            return int(val_numeric)
    except: 
        return val

# Big TODO: deduplicate this with h5ad-to-arrow
def arrow_to_json(arrow_file, **kwargs):
    umap_json  = kwargs['umap_json']
    leiden_json  = kwargs['leiden_json']
    cell_sets_json = kwargs['cell_sets_json']
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df_items = df.T.to_dict().items()

    # It is possible for the cluster names to not be integers.
    df['leiden'] = df['leiden'].astype(str)
    leiden_clusters = sorted(df['leiden'].unique(), key=try_str_to_int)

    id_to_umap = {
        k: {
            "mappings": {"UMAP": [v['umap_x'], v['umap_y']]},
            "factors": {"Leiden Cluster": v['leiden']}
        }
        for (k,v) in df_items
    }
    pretty_json_umap = json.dumps(id_to_umap).replace('}},', '}},\n')
    with open(umap_json, 'w') as f:
        f.write(pretty_json_umap)

    id_to_factors = {
        'Leiden Cluster': {
            'map': leiden_clusters,
            'cells': { k: leiden_clusters.index(v['leiden']) for (k,v) in df_items }
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
                        "set": df.loc[df['leiden'] == cluster].index.values.tolist(),
                    }
                    for cluster in leiden_clusters
                ]
            }
        ]
    }
    with open(cell_sets_json, 'w') as f:
        json.dump(cell_sets, f, indent=1)


def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(exist_ok=True, parents=True)
    csv_glob = '**/umap_coords_clusters.csv'
    csvs = list(input_dir.glob(csv_glob))
    if not csvs:
        raise Exception(f'No match for {csv_glob} in {input_dir}')
    for input_path in csvs:
        arrow_name = input_path.with_suffix('.arrow').name
        arrow_path = output_dir / arrow_name
        csv_to_arrow(input_path, arrow_path)
        arrow_to_csv(arrow_path, arrow_path.with_suffix('.csv'))
        arrow_to_json(
            arrow_file=arrow_path,
            umap_json=arrow_path.with_suffix('.cells.json'),
            leiden_json=arrow_path.with_suffix('.factors.json'),
            cell_sets_json=arrow_path.with_suffix('.cell-sets.json'),
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Transform raw pipeline CSV into Arrow, JSON, and standardized CSV.',
    )
    parser.add_argument(
        '--input_dir',
        required=True,
        help='directory containing csv files to read',
        type=Path,
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='directory where arrow files should be written',
        type=Path,
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
