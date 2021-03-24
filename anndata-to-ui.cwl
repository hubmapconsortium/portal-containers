#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', './anndata-zarr']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-anndata-to-ui:0.0.2
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 5
        prefix: --input_dir
outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: anndata-zarr
