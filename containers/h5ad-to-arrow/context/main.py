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
    umap = ann_data.obsm['X_umap']
    index = ann_data.obs.index

    df = DataFrame(data=umap, index=index)
    table = pa.Table.from_pandas(df)

    writer = pa.RecordBatchFileWriter(arrow_file, table.schema)
    writer.write(table)
    writer.close()


def arrow_to_csv(arrow_file, csv_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df.to_csv(csv_file)


def arrow_to_json(arrow_file, json_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    id_to_pair = {
        k: { "mappings": {"UMAP": [v[0], v[1]]} }
        for (k,v) in df.T.to_dict().items()
    }
    pretty_json = json.dumps(id_to_pair).replace(']}},', ']}},\n')
    with open(json_file, 'w') as f:
        f.write(pretty_json)


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
        arrow_to_json(arrow_path, arrow_path.with_suffix('.json'))


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