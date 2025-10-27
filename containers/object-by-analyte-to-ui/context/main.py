import argparse
from pathlib import Path
from os import path, walk
import json
import shutil
from scipy import sparse
from mudata import read_h5mu, MuData
from mudata._core.mudata import ModDict
from repro_zipfile import ReproducibleZipFile
import scanpy as sc

NUM_MARKER_GENES_TO_VISUALIZE = 5
VAR_CHUNK_SIZE = 10
INPUT_FILE_NAMES = ["secondary_analysis.h5mu"]

# Keeping in case it's necessary
# uberon_mapping_for_object_types = {
#     'cell': 'CL:0000000',
#     'tissue': 'UBERON:0000479',
#     'organ': 'UBERON:0000062',
#     'organism': 'UBERON:0000468'
# }


def get_annotations(modality: ModDict) -> dict:
    """
    Retrieves the list of annotation methods for the given modality.
    """
    if "annotation" in modality.obsm:
        return list(modality.obsm["annotation"].keys())
    else:
        return []


def get_object_types(modality: ModDict) -> list:
    """
    Retrieves the list of object types for the given modality.
    """
    return list(set(modality.obs["object_type"]))


def get_modality_metadata(modality: ModDict, key: str) -> dict:
    """
    Retrieves the metadata for a specific modality.
    """
    return {
        "name": key,
        "n_obs": modality.n_obs,
        "n_vars": modality.n_vars,
        "obs_keys": sorted(list(modality.obs.keys())),
        "obsm_keys": sorted(list(modality.obsm.keys())),
        "var_keys": sorted(list(modality.var.keys())),
        "annotations": sorted(get_annotations(modality)),
    }


def get_metadata(mdata: MuData) -> dict:
    """
    Retrieves detailed metadata for the mudata file and its modalities.
    Provides additional information for downstream visualization.
    """
    return {
        "shape": mdata.shape,
        "n_obs": mdata.n_obs,
        "n_vars": mdata.n_vars,
        "obs_keys": sorted(list(mdata.obs.keys())),
        "obsm_keys": sorted(list(mdata.obsm.keys())),
        "var_keys": sorted(list(mdata.var.keys())),
        "epic_type": sorted(list(set(mdata.uns["epic_type"]))),
        "modalities": [
            get_modality_metadata(mod, key) for key, mod in mdata.mod.items()
        ],
    }


def get_calculated_metadata(mdata: MuData) -> dict:
    """
    Retrieves the calculated metadata for the mudata file.
    This is provided to the search index's calculated_metadata for the dataset.
    """

    modalities = mdata.mod.values()
    return {
        # Flatten and deduplicate the returned list
        "object_types": sorted(
            list(
                set(
                    object_type
                    for modality in modalities
                    for object_type in modality.obs["object_type"]
                )
            )
        ),
        "annotation_tools": sorted(
            list(
                set(
                    tool
                    for modality in modalities
                    for tool in get_annotations(modality)
                )
            )
        ),
        "epic_type": sorted(list(set(mdata.uns["epic_type"]))),
    }


# Basic conversion from mudata to zarr zip
def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(exist_ok=True)

    # Find the h5mu file to process (only one per directory)
    h5mu_file = None

    # First, check for the specific file
    for candidate_file in INPUT_FILE_NAMES:
        input_path = input_dir / candidate_file
        if input_path.exists():
            h5mu_file = candidate_file
            break

    # If not found, look for any file ending with 'processed.h5mu'
    if h5mu_file is None:
        processed_files = list(input_dir.glob("*_processed.h5mu"))
        if processed_files:
            h5mu_file = processed_files[0].name  # Take the first one found

    if h5mu_file is None:
        print("No suitable h5mu file found in input directory.")
        return

    mdata = read_h5mu(input_dir / h5mu_file)
    if mdata is None:
        print(f"Failed to read {h5mu_file}.")
        return

    print("################################################")
    print(f"Processing {h5mu_file} Metadata")
    metadata = get_metadata(mdata)
    print(f"Processing {h5mu_file} Calculated Metadata")
    calculated_metadata = get_calculated_metadata(mdata)

    # Convert sparse layer matrices to CSC format for performance
    for key, modality in mdata.mod.items():
        for layer in modality.layers:
            if isinstance(modality.layers[layer], sparse.spmatrix):
                modality.layers[layer] = modality.layers[layer].tocsc()

        if modality.n_vars > 300:
            # Select 100 highest dispersion variables using scanpy
            vars_subset_size = 100

            # Use scanpy's highly_variable_genes function
            sc.pp.highly_variable_genes(
                modality, n_top_genes=vars_subset_size, subset=False, inplace=True
            )

            # Mark the top highly variable genes
            modality.var["top_highly_variable"] = modality.var["highly_variable"]

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
    chunks = (
        (mdata.shape[0], VAR_CHUNK_SIZE) if mdata.shape[1] >= VAR_CHUNK_SIZE else None
    )

    zarr_path = output_dir / (Path(h5mu_file).stem + ".zarr")
    mdata.write_zarr(zarr_path, chunks=chunks)
    print("Zarr store created")

    zip_path = output_dir / (Path(h5mu_file).stem + ".zarr.zip")

    with ReproducibleZipFile(zip_path, "w") as zf:
        for root, dirs, files in walk(zarr_path):
            for file in sorted(files):
                full_path = path.join(root, file)
                arcname = path.relpath(full_path, start=zarr_path)
                zf.write(full_path, arcname=arcname)
        print("Zip zarr created")

    # # Save potentially useful information to a json for easy access
    metadata_dir = output_dir / "secondary_analysis_metadata.json"
    with open(metadata_dir, "w") as f:
        json.dump(metadata, f, sort_keys=True)
        print("Additional Metadata JSON created.")

    calculated_metadata_dir = output_dir / "calculated_metadata.json"
    with open(calculated_metadata_dir, "w") as f:
        json.dump(calculated_metadata, f, sort_keys=True)
        print("Calculated Metadata JSON created.")

    # Clean up non-zipped zarr store
    shutil.rmtree(zarr_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Transform Object by Analyte MuData into zarr."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        type=Path,
        help="directory containing MuData .h5mu file to read. Will process 'secondary_analysis.h5mu' if it exists, otherwise the first file ending with 'processed.h5mu' found.",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="directory where MuData zarr files should be written",
    )
    args = parser.parse_args()
    main(args.input_dir, args.output_dir)
