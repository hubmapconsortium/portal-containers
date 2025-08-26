import argparse
from pathlib import Path
from os import mkdir
import json

from mudata import MuData
from anndata import AnnData
from numpy import float32, random, array
from pandas import DataFrame
from scipy.sparse import csr_matrix


def load_parameters(parameters_path):
    """Load parameters from JSON file."""
    with open(parameters_path, 'r') as f:
        return json.load(f)


def create_sample_data(n_obs, n_vars, seed=42):
    """Create sample data arrays with specified dimensions."""
    random.seed(seed)
    # Create sparse-like data with mostly zeros and some positive values
    data = random.poisson(lam=1.0, size=(n_obs, n_vars)).astype(float32)
    return data


def create_obs_dataframe(obs_keys, n_obs, modality_prefix=""):
    """Create observation DataFrame with specified keys."""
    obs_data = {}

    for key in obs_keys:
        if modality_prefix:
            clean_key = key.replace(f"{modality_prefix}:", "")
        else:
            clean_key = key

        if any(x in clean_key.lower() for x in ['leiden', 'cluster']):
            # Cluster assignments (categorical)
            obs_data[clean_key] = [f"cluster_{i % 3}" for i in range(n_obs)]
        elif 'object_type' in clean_key.lower():
            # Always set object_type to 'cell'
            obs_data[clean_key] = ['cell'] * n_obs
        elif any(x in clean_key.lower() for x in
                 ['age', 'bmi', 'weight', 'height']):
            # Numeric data
            obs_data[clean_key] = [20 + (i * 5) % 60 for i in range(n_obs)]
        elif any(x in clean_key.lower() for x in ['sex', 'race', 'tissue']):
            # Categorical data
            categories = ['A', 'B', 'C']
            obs_data[clean_key] = [categories[i % len(categories)]
                                   for i in range(n_obs)]
        elif 'id' in clean_key.lower():
            # ID fields
            obs_data[clean_key] = [f"{clean_key}_{i:03d}"
                                   for i in range(n_obs)]
        elif any(x in clean_key.lower() for x in
                 ['score', 'density', 'mean']):
            # Float scores
            obs_data[clean_key] = [0.1 + (i * 0.1) % 1.0 for i in range(n_obs)]
        else:
            # Default string values
            obs_data[clean_key] = [f"value_{i}" for i in range(n_obs)]

    return DataFrame(obs_data,
                     index=[f"cell_{i:03d}" for i in range(n_obs)])


def create_obsm_data(obsm_keys, n_obs, seed=42):
    """Create obsm multi-dimensional arrays for specified keys."""
    random.seed(seed)
    obsm_data = {}

    for key in obsm_keys:
        if key in ['X_pca', 'X_umap']:
            # Dimensional reduction data (2D for UMAP, could be higher for PCA)
            n_dims = 2 if 'umap' in key.lower() else 10
            obsm_data[key] = random.normal(
                0, 1, size=(n_obs, n_dims)).astype(float32)
        elif key == 'annotation':
            # For annotation, we need to create a structured array or matrix
            # that can be accessed with .keys() method
            # Create a simple 2D array for now
            obsm_data[key] = random.normal(
                0, 1, size=(n_obs, 2)).astype(float32)
        elif 'azimuth_label' in key:
            # Cell type labels
            cell_types = [f"cell_type_{i % 5}" for i in range(n_obs)]
            obsm_data[key] = array(cell_types, dtype=object)
        elif 'leiden' in key:
            # Leiden cluster assignments
            clusters = [f"cluster_{i % 3}" for i in range(n_obs)]
            obsm_data[key] = array(clusters, dtype=object)
        else:
            # Default: create some random data
            obsm_data[key] = random.normal(
                0, 1, size=(n_obs, 2)).astype(float32)

    return obsm_data


def create_var_dataframe(var_keys, n_vars, modality_name=""):
    """Create variable DataFrame with specified keys."""
    var_data = {}

    for key in var_keys:
        if 'symbol' in key.lower() or 'gene' in key.lower():
            # Gene symbols
            var_data[key] = [f"GENE_{i:05d}" for i in range(n_vars)]
        elif 'uniprot' in key.lower():
            # UniProt IDs
            var_data[key] = [f"P{i:05d}" for i in range(n_vars)]
        elif any(x in key.lower() for x in ['mean', 'std']):
            # Numeric statistics
            var_data[key] = [0.1 + (i * 0.01) % 10.0 for i in range(n_vars)]
        elif 'n_cells' in key.lower():
            # Cell counts
            var_data[key] = [i % 100 + 1 for i in range(n_vars)]
        else:
            # Default values
            var_data[key] = [f"var_{i}" for i in range(n_vars)]

    prefix = f"{modality_name}_" if modality_name else ""
    return DataFrame(var_data,
                     index=[f"{prefix}feature_{i:05d}"
                            for i in range(n_vars)])


def create_h5mu(h5mu_path, parameters_path):
    """Create h5mu file based on parameters."""
    params = load_parameters(parameters_path)

    # Create modalities first
    modalities = {}
    for modality in params['modalities']:
        mod_name = modality['name']
        mod_n_obs = modality['n_obs']
        mod_n_vars = modality['n_vars']

        # Create data matrix for this modality
        mod_data = create_sample_data(mod_n_obs, mod_n_vars)

        # Create obs and var DataFrames for this modality
        mod_obs = create_obs_dataframe(modality['obs_keys'], mod_n_obs,
                                       mod_name)
        mod_var = create_var_dataframe(modality['var_keys'], mod_n_vars,
                                       mod_name)

        # Create AnnData for this modality
        adata = AnnData(
            X=csr_matrix(mod_data),
            obs=mod_obs,
            var=mod_var
        )

        # Add obsm data if specified in parameters
        if 'obsm_keys' in modality:
            obsm_data = create_obsm_data(modality['obsm_keys'], mod_n_obs)
            for key, value in obsm_data.items():
                if key == 'annotation':
                    # Handle annotation specially - create a DataFrame
                    # that can be accessed with .keys()
                    annotation_data = {}
                    if 'annotations' in modality:
                        for annotation_key in modality['annotations']:
                            if annotation_key == 'azimuth_label':
                                annotation_data[annotation_key] = [
                                    f"cell_type_{i % 5}" for i in range(
                                        mod_n_obs)]
                            elif annotation_key == 'leiden':
                                annotation_data[annotation_key] = [
                                    f"cluster_{i % 3}" for i in range(
                                        mod_n_obs)]
                            else:
                                # Default annotation data
                                annotation_data[annotation_key] = [
                                    f"annotation_{i}" for i in range(
                                        mod_n_obs)]

                    # Create DataFrame with proper index
                    annotation_df = DataFrame(
                        annotation_data,
                        index=[f"cell_{i:03d}" for i in range(mod_n_obs)]
                    )
                    adata.obsm[key] = annotation_df
                else:
                    adata.obsm[key] = value

        # Add some layers for variety
        adata.layers['raw'] = csr_matrix(mod_data * 0.8)
        adata.layers['normalized'] = csr_matrix(mod_data / mod_data.max())

        modalities[mod_name] = adata

    # Create MuData object
    h5mu = MuData(modalities)

    # Update to get the proper dimensions
    h5mu.update()

    # Now create global obs based on actual dimensions
    global_obs = create_obs_dataframe(params['obs_keys'], h5mu.n_obs)

    # For global var, we need to match the actual concatenated variables
    # Use the existing var index from MuData
    global_var_data = {}
    for key in params['var_keys']:
        if 'symbol' in key.lower() or 'gene' in key.lower():
            # Gene symbols - use existing index names
            global_var_data[key] = [f"GENE_{i:05d}"
                                    for i in range(h5mu.n_vars)]
        elif 'uniprot' in key.lower():
            # UniProt IDs
            global_var_data[key] = [f"P{i:05d}" for i in range(h5mu.n_vars)]
        elif any(x in key.lower() for x in ['mean', 'std']):
            # Numeric statistics
            global_var_data[key] = [0.1 + (i * 0.01) % 10.0
                                    for i in range(h5mu.n_vars)]
        elif 'n_cells' in key.lower():
            # Cell counts
            global_var_data[key] = [i % 100 + 1 for i in range(h5mu.n_vars)]
        else:
            # Default values
            global_var_data[key] = [f"var_{i}" for i in range(h5mu.n_vars)]

    # Create global var DataFrame with the correct index
    global_var = DataFrame(global_var_data, index=h5mu.var.index)

    # Set global obs and var
    h5mu.obs = global_obs
    h5mu.var = global_var

    # Add epic_type to uns if specified in parameters
    if 'epic_type' in params:
        h5mu.uns['epic_type'] = params['epic_type']

    # Final update and write
    h5mu.var_names_make_unique()
    h5mu.write(h5mu_path)


def main(output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    h5mu_path = Path(output_dir) / "secondary_analysis.h5mu"
    parameters_path = Path(__file__).parent / "test-input-parameters.json"
    create_h5mu(h5mu_path, parameters_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
            Creates a minimal mudata input fixture based on parameters.json,
            with a structure similar to what we'll actually receive.
        """
    )
    parser.add_argument(
        "dest", help="Directory where test files should be written")
    args = parser.parse_args()
    main(args.dest)
