import argparse
from pathlib import Path
from os import mkdir

from anndata import read_h5ad


def create_h5ad(h5ad_path):
    print(h5ad_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5ad_path = Path(output_dir) / 'fake.h5ad'
    create_h5ad(h5ad_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=f'''
            Creates a minimal anndata input fixture,
            with a structure similar to what we'll actually receive.
        ''')
    parser.add_argument(
        'dest',
        help='Directory where arrow files should be written')
    args = parser.parse_args()
    main(args.dest)
