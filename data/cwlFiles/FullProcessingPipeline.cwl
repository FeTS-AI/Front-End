cwlVersion: v1.0
class: CommandLineTool
baseCommand: FullProcessingPipeline
inputs:
  inputCSV:
    type: File
    label: Input CSV file
    inputBinding:
      position: 1
      prefix: -i
    doc: Input CSV file which contains paths to structural images.Headers should be 'PatientID,T1,T1GD,T2,T2FLAIR'.
  outputDir:
    type: Directory
    label: Directory
    inputBinding:
      position: 1
      prefix: -o
    doc: "Output directory for final output.This will write 2 folders: 'DataForFeTS' and 'DataForQC'.Former contains only the files needed for FeTS inference/training and .latter contains all intermediate files from this processing."
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
  archs:
    type: string?
    label: 3DResUNet,deepMedic,deepscan
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
  LoggingDir:
    type: Directory?
    label: Dir with write access
    inputBinding:
      position: 1
      prefix: -L
    doc: Location of logging directory.
hints:
  SoftwareRequirement:
    packages:
      FullProcessingPipeline:
        version:
          - 0.0.2