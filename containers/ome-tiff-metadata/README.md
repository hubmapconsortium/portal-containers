# ome-tiff-metadata

This docker container creates a JSON object for metadata (physical sizes and units) extracted from each TIFF within an input directory. This is needed particularly for laying segmentation masks over base-images when there is a misalignment between the physical sizes fo both images. Note that, even with same pixel sizes, if physical sizes are different, misalignment can occur. This metadata is useful to add scaling to the segmentation mask when visualized using Vitessce.

## Input

The input to the container is one or more ome-tiff image files.

## Output

The output is a json file that includes an object with following structure for every input ome-tiff image.

```
{"PhysicalSizeX": 0.5, "PhysicalSizeY": 0.5, "PhysicalSizeUnitX": "\u00b5m", "PhysicalSizeUnitY": "\u00b5m"}

```

## Normalization

None

## Example

Example of a hubmap dataset using this container for metadata generation for Vitessce (visualization) would be
`TODO`.
