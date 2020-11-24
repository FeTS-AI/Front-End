cwlVersion: v1.0
class: CommandLineTool
baseCommand: OpenFLCLI
inputs:
  dataDir:
    type: Directory
    label: Dir with Read/Write access
    inputBinding:
      position: 1
      prefix: -d
    doc: Input data directory.
  modelName:
    type: File
    label: Model file
    inputBinding:
      position: 1
      prefix: -m
    doc: Input model weights file.
  training:
    type: boolean
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -t
    doc: Whether performing training or inference.1==Train and 0==Inference.
  LoggingDir:
    type: Directory
    label: Dir with write access
    inputBinding:
      position: 1
      prefix: -L
    doc: Location of logging directory.
  runtest:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -rt
    doc: Runs the tests.
  cwl:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -cwl
    doc: Generates a .cwl file for the software.
  gpu:
    type: boolean?
    label: 0-1
    inputBinding:
      position: 1
      prefix: -g
    doc: Whether to run the process on GPU or not.Defaults to '0'.
  colName:
    type: string?
    label: none
    inputBinding:
      position: 1
      prefix: -c
    doc: Common name of collaborator.Required for training.
hints:
  SoftwareRequirement:
    packages:
      OpenFLCLI:
        version:
          - 0.0.1