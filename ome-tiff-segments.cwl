#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', 'output_ome_segments', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-ome-tiff-segments:0.0.2
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  json:
    type: Directory
    outputBinding:
      glob: output_ome_segments