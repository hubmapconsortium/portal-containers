#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', './anndata-zarr', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-sprm-to-anndata:0.0.13
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: anndata-zarr
