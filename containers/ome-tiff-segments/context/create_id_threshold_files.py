import re
import pandas as pd
import numpy as np
import json
import zarr
from anndata import read_zarr
from numcodecs import JSON 
from ome_types import model, OME

def create_aoi_table(ome: OME, output_path:str):
    '''
    Creates a csv file of ids, that is, AOI-id, ROI-id, and segment name.

    Args:
        ome (OME): OME object containing ROI data.
        output_path(str): location where the csv file should be created

    Returns:
        creates a csv file with the ids
    '''                
    rows = []  
    try:
        for roi in ome.rois:
            roi_id = roi.id
            if roi.union:
                for shape in roi.union:
                    if isinstance(shape, model.Mask):
                        mask_id = shape.id
                        text = shape.text if shape.text else 'N/A' 
                        rows.append({
                            'roi_id': roi_id,
                            'aoi_id': mask_id,
                            'segment': text,
                            'aoi_id_numeric':int(re.search(r'\d+', shape.id).group())
                        })
    except Exception as e:
        raise ValueError(f'Error in extracting ids from ROIs {str(e)}')
    
    zarr_file_path=f'{output_path}/aoi.zarr'
    convert_to_zarr(rows, zarr_file_path)
    # The index need to be converted to the aoi_numeric_id for Vitessce to read the segmentation ome.tiff
    convert_index_in_zarr(zarr_file_path)



def create_roi_table(ome:OME, output_path:str):
    '''
    Creates a CSV file containing ROI names and their associated threshold values from an OME object.

    The function processes the ROIs (Regions of Interest) and structured annotations in the OME object.
    For each ROI, it extracts its associated annotation reference and then retrieves the corresponding
    annotation's elements to identify the ROI name and threshold values.

    Args:
        ome (OME): The OME object containing ROIs and structured annotations.
        output_path(str): location where the csv file should be created

    Returns:
        None: The function writes the processed ROI data to a CSV file named 'roi.csv'.

    '''
    roi_data = {}
    try:
        for roi in ome.rois:
            annotation_id = roi.annotation_refs[0].id
            roi_id = roi.id
            roi_data[annotation_id] = roi_id
    except Exception as e:
        raise ValueError(f'Error in extracting annotation refs from ROIs {str(e)}')

    rows = []
    try:
        for annotation in ome.structured_annotations:
            roi_name = None
            thresholds = {}

            if annotation.id in roi_data:
                any_elements = annotation.value.any_elements
                for element in any_elements:
                        for child in element.children:
                            qname = child.qname
                            text = child.text
                            if 'RoiName' in qname:
                                roi_name = roi_data[annotation.id]
                            elif 'threshold' in qname.lower():
                                threshold_name = qname.lower().replace('threshold', '').strip()
                                threshold_name_formatted = re.sub(r'\{.*\}', '', threshold_name)
                                if threshold_name_formatted:
                                    thresholds[f'{threshold_name_formatted}_threshold'] = int(text)
                                else: 
                                    thresholds['threshold'] = int(text)
            if roi_name:
                row = {'roi_id': roi_name}
                row.update(thresholds)
                rows.append(row)

    except Exception as e:
        raise ValueError(f'Error in extracting channel thresholds from annotations {str(e)}')
    zarr_file_path=f'{output_path}/roi.zarr'
    convert_to_zarr(rows, zarr_file_path)

def create_mask_vertices_from_rois(ome: OME, output_path:str):
    '''
    Create a dictionary of segment names with corresponding polygon vertices from the ROIs in the OME object.

    Args:
        ome (OME): OME object containing ROI data.
        output_path(str): location where the json file should be created

    Returns:
        Creates a json file from the dictionary of segment names and polygon vertices
    '''
    roi_ids = {}
    try:
        for roi in ome.rois:
            for shape in roi.union:
                if isinstance(shape, model.Polygon):
                    roi_id = roi.id
                    if roi_id not in roi_ids:
                        roi_ids[roi_id] = []
                    
                    points = np.array([list(map(float, p.split(','))) for p in shape.points.split()])
                    
                    roi_ids[roi_id].append(points.tolist())
    except Exception as e:
        raise ValueError(f'Error in extracting polygon vertices for ROIs {str(e)}')

    with open(f'{output_path}/obsSegmentations.json', 'w') as json_file:
            json.dump(roi_ids, json_file, indent=4)
            print(f'obsSegmentations.json saved')


def convert_to_zarr(rows, zarr_file_path):
    df = pd.DataFrame(rows) 
    try:
        zarr_store = zarr.open(zarr_file_path, mode='w')
        obs_group = zarr_store.create_group('obs')

        for column in df.columns:
            data = df[column].values
            # Use JSON codec for string/mixed types, otherwise store without a codec
            if data.dtype.kind in {'O', 'U'}:
                obs_group.array(column, data, object_codec=JSON())
            else:
                obs_group.array(column, data)

        print(f'obs saved at {zarr_file_path}')

    except Exception as e:
        raise ValueError(f'Error storing rows as Zarr stores: {str(e)}')
    
def convert_index_in_zarr(zarr_file_path):
    try:
        adata = read_zarr(zarr_file_path)
        adata.obs['aoi_id_numeric'] = adata.obs['aoi_id_numeric'].astype(str)
        adata.obs = adata.obs.set_index("aoi_id_numeric")
        adata.write_zarr(zarr_file_path)
        print('Index converted in aoi.zarr')
    except Exception as e:
        raise ValueError(f'Error converting aoi.zarr\'s index: {str(e)}')