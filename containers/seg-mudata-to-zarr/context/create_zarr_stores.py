import mudata as md
import pandas as pd
import numpy as np

VAR_CHUNK_SIZE = 10
index_mapping = None
mask_name_col = "mask name"

# Function to extract numeric part from index
def extract_numeric_index(index):
    numeric_part = index.str.extract(r'(\d+)', expand=False).astype(int)
    return numeric_part


def convert_obs (mask_data):
    global index_mapping 
    mask_data.obs = mask_data.obs.loc[~mask_data.obs.index.duplicated(keep='first')]
    numeric_index = extract_numeric_index(mask_data.obs.index)
    index_mapping = pd.Series(numeric_index.values, index=mask_data.obs.index)
    mask_data.obs.index = numeric_index
    return mask_data

def convert_obsm(obsm_matrix, obs_data):
    """
        Reindexes and deduplicates the obsm data to match obs_data.

    Parameters
    ----------
    obsm_matrix : DataFrame or nparray
        The data in the .obsm component of mudata 

    obs_data : DataFrame
        The .obs component of mudata
    Returns
    -------
        The reindexed and deduplicated obsm_matrix.
    """
    if isinstance(obsm_matrix, pd.DataFrame):
            obsm_matrix = obsm_matrix.reindex(index_mapping.index)
            obsm_matrix.index = index_mapping.values 
            # Getting error without this deduplication step, although did not find any duplicates
            obsm_matrix = obsm_matrix.loc[~obsm_matrix.index.duplicated(keep='first')]
            if not obs_data.index.equals(obsm_matrix.index):
                raise ValueError("Indices do not match between mask_data.obs and obsm_matrix after re-indexing.")
            return obsm_matrix
        # Some keys such as default are only 1-dimensional arrays
    elif isinstance(obsm_matrix, np.ndarray):
        return obsm_matrix
    else:
        raise ValueError("The shape of obsm_matrix does not match the number of observations in mask_data.obs.")
            
    

def create_zarr_for_masks(mdata, output_path):
    """
        Creates Zarr stores for each unique mask in the given MuData object.

    Parameters
    ----------
    mdata : MuData
        A MuData object containing the mask data and the `obs` attribute to contain a column with mask names

    output_path : str or Path
        The directory where the Zarr stores for each mask will be written.
    
    Returns
    -------
        None
    """
    try:
        mask_names = mdata.obs[mask_name_col].unique()
        # Create Zarr stores for each mask name
        for mask_name in mask_names:
            mask_data = mdata[mdata.obs[mask_name_col] == mask_name.lower()]
            mask_data = convert_obs(mask_data)
            for key in mask_data.obsm.keys():
                obsm_matrix = mask_data.obsm[key]
                mask_data.obsm[key] = convert_obsm(obsm_matrix, mask_data.obs)
            for mod_key in mask_data.mod.keys():
                mod_data = mask_data.mod[mod_key]
                for obsm_key in mod_data.obsm.keys():
                    obsm_matrix = mod_data.obsm[obsm_key]
                    mod_data.obsm[obsm_key] = convert_obsm(obsm_matrix, mask_data.obs)

            chunks = (mask_data.shape[0], VAR_CHUNK_SIZE) if mask_data.shape[1] >= VAR_CHUNK_SIZE else None
            zarr_store_path = f'{output_path}/{mask_name}.zarr'
            mask_data.write_zarr(zarr_store_path, chunks=chunks)
            print(f'Created Zarr store for the mask: {mask_name}')
    except Exception as e:
        print(f'Error in getting mask names from mudata {str(e)}')
