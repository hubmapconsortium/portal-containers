#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', './mudata-zarr', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-mudata-to-ui:0.0.6
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: mudata-zarr
