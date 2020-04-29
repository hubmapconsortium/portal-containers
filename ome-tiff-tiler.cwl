#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/opt/process.sh', '-o', '.', '-i']
hints:
  DockerRequirement:
    dockerPull: hubmap/portal-container-ome-tiff-tiler:0.0.1
inputs:
  input_directory:
    type: Directory
    inputBinding:
        position: 5
outputs:
  tiff:
    type: File
    outputBinding:
      glob: '*.ome.tif*'
  n5:
    type: Directory
    outputBinding:
      glob: '*.n5'
