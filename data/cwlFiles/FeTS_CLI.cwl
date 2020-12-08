cwlVersion: v1.0
class: CommandLineTool
baseCommand: FeTS_CLI
inputs:
  dataDir:
    type: Directory
    label: Dir with Read/Write access
    inputBinding:
      position: 1
      prefix: -d
    doc: Input data directory.
  training:
    type: boolean
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -t
    doc: Whether performing training or inference.1==Train and 0==Inference.
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
  LoggingDir:
    type: Directory?
    label: Dir with write access
    inputBinding:
      position: 1
      prefix: -L
    doc: Location of logging directory.
  archs:
    type: string?
    label: 3DResUNet,3DUNet,deepMedic,deepscan
    inputBinding:
      position: 1
      prefix: -a
    doc: "The architecture(s) to infer/train on.Only a single architecture is supported for training.Comma-separated values for multiple options.Defaults to: 3dresunet."
  labelFuse:
    type: string?
    label: STAPLE,ITKVoting,SIMPLE,MajorityVoting
    inputBinding:
      position: 1
      prefix: -lF
    doc: "The label fusion strategy to follow for multi-arch inference.Comma-separated values for multiple options.Defaults to: STAPLE."
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
      FeTS_CLI:
        version:
          - 0.0.1