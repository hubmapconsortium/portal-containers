import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ

from anndata import read_h5ad
import pyarrow as pa


def h5ad_to_arrow(h5ad_file, arrow_file):
    ann_data = read_h5ad(h5ad_file)
    dataframe = ann_data.to_df()
    table = pa.Table.from_pandas(dataframe)

    writer = pa.RecordBatchFileWriter(arrow_file, table.schema)
    writer.write(table)
    writer.close()


def arrow_to_text(arrow_file, text_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df.to_csv(text_file)


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
        arrow_to_text(arrow_path, arrow_path.with_suffix('.csv'))


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
