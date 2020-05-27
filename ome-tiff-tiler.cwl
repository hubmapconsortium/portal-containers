#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/opt/process.sh', '-o', 'output_tiff_pyramids']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-ome-tiff-tiler:0.0.2
inputs:
  input_directory:
    type: Directory
    inputBinding:
      prefix: -i
      position: 3
  workers:
    type: int
    inputBinding:
      prefix: -w
      position: 4
  rgb:
    type: boolean
    inputBinding:
      prefix: -r
      position: 5
outputs:
  output_directory:
    type: Directory
    outputBinding:
      glob: output_tiff_pyramids
