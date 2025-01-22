#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', 'image_metadata', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-ome-tiff-metadata:0.0.1
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  json:
    type: Directory
    outputBinding:
      glob: image_metadata
