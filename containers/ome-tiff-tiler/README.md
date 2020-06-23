# ome-tiff-tiler

This docker container contains a tiler based on [bioformats2raw](https://github.com/glencoesoftware/bioformats2raw) and [raw2ometiff](https://github.com/glencoesoftware/raw2ometiff).
The CWL allows for an input directory `-i` and optional bindings for workers with a `-w` and then `rgb` with `-r`. The input is a directory of `OME-TIFF` files and the output is a `OME-TIFF` pyramid and a `n5` directory for quicker local analysis.
