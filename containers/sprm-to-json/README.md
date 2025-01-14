# sprm-to-json

This container takes [SPRM](https://docs.google.com/document/d/1c7UR0Pe1newpVhQY2HEFkfV8O7GAj9Vk4XnuSiSnDeY/edit) output CSVs and converts them into Vitessce JSON which conforms to our [schemas](https://github.com/hubmapconsortium/vitessce/tree/master/src/schemas).

## Input
The input to the container SPRM output csvs.

## Output
The output includes json files for cells and cell-sets that conform to Vitessce schemas.
 
## Normalization
None

## Example
Example of a hubmap dataset using this container for data conversion would be 
`https://portal.hubmapconsortium.org/browse/dataset/776bd778f20409d8ce2dd8f5bf5789ae` 