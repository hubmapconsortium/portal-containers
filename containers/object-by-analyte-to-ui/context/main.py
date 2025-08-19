import argparse
from pathlib import Path
from os import path, walk
import json
import shutil
from numpy import asarray
from scipy import sparse
from mudata import read_h5mu, MuData
from mudata._core.mudata import ModDict
from repro_zipfile import ReproducibleZipFile

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
INPUT_FILE_NAMES = ["secondary_analysis.h5mu"]


# Retrieves the metadata for individual modalities
def get_modality_metadata(modality: ModDict, key: str) -> dict:
    return {
        'name': key,
        'n_obs': modality.n_obs,
        'n_vars': modality.n_vars,
        'obs_keys': list(modality.obs.keys()),
        'var_keys': list(modality.var.keys())
    }


# Retrieves the metadata for the mudata file and its modalities
def get_metadata(mdata: MuData) -> dict:
    return {
        'shape': mdata.shape,
        'n_obs': mdata.n_obs,
        'n_vars': mdata.n_vars,
        'obs_keys': list(mdata.obs.keys()),
        'var_keys': list(mdata.var.keys()),
        'modalities': [get_modality_metadata(mod, key)
                       for key, mod in mdata.mod.items()]
    }


# Basic conversion from mudata to zarr zip
def main(input_dir: str, output_dir: str):
    output_dir.mkdir(exist_ok=True)
    for h5mu_file in INPUT_FILE_NAMES:
        # Check if input file exists, skip it if it doesn't exist
        input_path = path.join(input_dir, h5mu_file)
        if not path.exists(input_path):
            print(f"Input file {h5mu_file} does not exist in input directory.")
            continue
        mdata = read_h5mu(input_dir / h5mu_file)
        if mdata is None:
            print(f"Failed to read {h5mu_file}.")
            continue

        print("################################################")
        print(f"Processing {h5mu_file}")
        metadata = get_metadata(mdata)

        # Convert sparse layer matrices to CSC format for performance
        for key, modality in mdata.mod.items():
            for layer in modality.layers:
                if isinstance(modality.layers[layer], sparse.spmatrix):
                    modality.layers[layer] = modality.layers[layer].tocsc()

        # If the main matrix is sparse, it's best for performance to
        # use non-sparse formats to keep the portal responsive.
        # In the future, we should be able to use CSC sparse data natively
        # and get equal performance:
        # https://github.com/theislab/anndata/issues/524
        # for data_layer in mdata.mod:
        #     print('data_layer: {data_layer}')
        #     if isinstance(data_layer.X, sparse.spmatrix):
        #         data_layer.X = asarray(data_layer.X.todense())

        # # It is now possible for adata.X to be empty and have shape (0, 0)
        # # so we need to check for that here, otherwise there will
        # # be a division by zero error during adata.write_zarr
        # # Reference: https://github.com/hubmapconsortium/salmon-rnaseq/blob/dfb0e2a/bin/analysis/scvelo_analysis.py#L69
        chunks = (
            mdata.shape[0], VAR_CHUNK_SIZE
        ) if mdata.shape[1] >= VAR_CHUNK_SIZE else None

        zarr_path = output_dir / (Path(h5mu_file).stem + ".zarr")
        mdata.write_zarr(zarr_path, chunks=chunks)
        print('Zarr store created')

        zip_path = output_dir / (Path(h5mu_file).stem + ".zarr.zip")

        with ReproducibleZipFile(zip_path, "w") as zf:
            for root, dirs, files in walk(zarr_path):
                for file in sorted(files):
                    full_path = path.join(root, file)
                    arcname = path.relpath(full_path, start=zarr_path)
                    zf.write(full_path, arcname=arcname)
            print("Zip zarr created")

        # # Save potentially useful information to a json for easy access
        metadata_dir = output_dir / (Path(h5mu_file).stem + "_metadata.json")
        with open(metadata_dir, "w") as f:
            json.dump(metadata, f)
            print("JSON Metadata created.")

        # Clean up non-zipped zarr store
        shutil.rmtree(zarr_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"Transform Object by Analyte MuData into zarr.")
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing MuData .h5mu files to read",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where MuData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
