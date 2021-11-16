import argparse
from glob import glob
from pathlib import Path
from os import mkdir, environ
import json

from anndata import read_h5ad
import pyarrow as pa
from pandas import DataFrame


PREDICTED_ASCT_CELLTYPE = 'predicted.ASCT.celltype'
PREDICTED_ASCT_CELLTYPE_SCORE = 'predicted.ASCT.celltype.score'

def has_cell_type_annotations(adata):
    return ('annotation_metadata' in adata.uns
        and 'is_annotated' in adata.uns['annotation_metadata']
        and adata.uns['annotation_metadata']['is_annotated']
        and PREDICTED_ASCT_CELLTYPE in adata.obs.columns.values.tolist()
        and PREDICTED_ASCT_CELLTYPE_SCORE in adata.obs.columns.values.tolist()
    )

def h5ad_to_arrow(h5ad_file, arrow_file):
    adata = read_h5ad(h5ad_file)
    umap = adata.obsm['X_umap'].transpose()
    leiden = adata.obs['leiden'].to_numpy().astype('uint8')
    index = adata.obs.index

    # TODO: condition on is_annotated
    adata_is_annotated = has_cell_type_annotations(adata)
    if adata_is_annotated:
        predicted_cell_type = adata.obs[PREDICTED_ASCT_CELLTYPE].astype(str)
    else:
        predicted_cell_type = None

    df = DataFrame(
        data={
            'umap_x': umap[0],
            'umap_y': umap[1],
            'leiden': leiden,
            **({
                PREDICTED_ASCT_CELLTYPE: predicted_cell_type,
            } if adata_is_annotated else {})
        },
        index=index
    )
    table = pa.Table.from_pandas(df)

    writer = pa.RecordBatchFileWriter(arrow_file, table.schema)
    writer.write(table)
    writer.close()


def arrow_to_csv(arrow_file, csv_file):
    df = pa.ipc.open_file(arrow_file).read_pandas()
    df.to_csv(csv_file)

def h5ad_to_json(h5ad_file, **kwargs):
    umap_json  = kwargs['umap_json']
    leiden_json  = kwargs['leiden_json']
    cell_sets_json = kwargs['cell_sets_json']
    adata = read_h5ad(h5ad_file)
    df = adata.obs
    df["umap_x"] = adata.obsm["X_umap"].T[0]
    df["umap_y"] = adata.obsm["X_umap"].T[1]
    df_items = df.T.to_dict().items()

    adata_is_annotated = has_cell_type_annotations(adata)

    id_to_umap = {
        k: {
            "mappings": {"UMAP": [v['umap_x'], v['umap_y']]},
            "factors": {
                "Leiden Cluster": str(int(v['leiden'])),
                **({
                    "Cell Type Prediction": str(v[PREDICTED_ASCT_CELLTYPE]),
                    "Cell Type Prediction Score": v[PREDICTED_ASCT_CELLTYPE_SCORE]
                } if adata_is_annotated else {})
            }
        }
        for (k,v) in df_items
    }
    pretty_json_umap = json.dumps(id_to_umap).replace('}},', '}},\n')
    with open(umap_json, 'w') as f:
        f.write(pretty_json_umap)

    leiden_clusters = sorted(df['leiden'].unique().astype('uint8'))

    if adata_is_annotated:
        predicted_cell_types = sorted(df[PREDICTED_ASCT_CELLTYPE].unique().astype(str))
    else:
        predicted_cell_types = None
    id_to_factors = {
        'Leiden Cluster': {
            'map': [str(cluster) for cluster in leiden_clusters],
            'cells': { k: v['leiden'] for (k,v) in df_items }
        },
        **({
            'Predicted ASCT Cell Type': {
                'map': [str(predicted_cell_type) for predicted_cell_type in predicted_cell_types],
                'cells': { k: v[PREDICTED_ASCT_CELLTYPE] for (k,v) in df_items }
            }
        } if adata_is_annotated else {})
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
                        "set": df.loc[df['leiden'] == cluster].index.values.tolist()
                    }
                    for cluster in leiden_clusters
                ]
            },
            *([
                {
                    "name": "Predicted ASCT Cell Type",
                    "children": [
                        {
                            "name": predicted_cell_type,
                            "set": df.loc[df[PREDICTED_ASCT_CELLTYPE] == predicted_cell_type].index.values.tolist()
                        }
                        for predicted_cell_type in predicted_cell_types
                    ]
                }
            ] if adata_is_annotated else [])
        ]
    }
    with open(cell_sets_json, 'w') as f:
        json.dump(cell_sets, f, indent=1)


def main(input_dir, output_dir):
    try:
        mkdir(output_dir)
    except FileExistsError:
        pass
    input_path = Path(input_dir + '/secondary_analysis.h5ad')
    arrow_name = input_path.with_suffix('.arrow').name
    arrow_path = Path(output_dir) / arrow_name
    h5ad_to_arrow(input_path, arrow_path)
    arrow_to_csv(arrow_path, arrow_path.with_suffix('.csv'))
    h5ad_to_json(
        h5ad_file=input_path,
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
