#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', 'output_offsets', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-ome-tiff-offsets:0.0.5
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  json:
    type: Directory
    outputBinding:
      glob: output_offsets
