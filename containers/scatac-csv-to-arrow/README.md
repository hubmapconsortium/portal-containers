# scatac-csv-to-arrow

This container takes the CSV from the [HuBMAP scATAC-seq pipeline](https://github.com/hubmapconsortium/sc-atac-seq-pipeline) and converts it to [Apache Arrow](https://arrow.apache.org/), as well as a normalized CSV, and Vitessce JSON which conforms to our [schemas](https://github.com/hubmapconsortium/vitessce/tree/master/src/schemas).


## Input
The input to the container is a csv file.


## Output
The output is an arrow file, along with csv and json files representing cells and cell sets for Vitessce.


## Normalization
To be continue...

## Example 
Example of a hubmap dataset using this container for arrow file generation for Vitessce (visualization) would be `HBM443.PWNF.335`.
