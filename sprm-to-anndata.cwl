#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
# TODO: Make main.py executable?
baseCommand: ['python', '/main.py', '--output_dir', './hubmap_ui', '--input_dir']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-sprm-to-anndata:0.0.1
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 6
outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: hubmap_ui
