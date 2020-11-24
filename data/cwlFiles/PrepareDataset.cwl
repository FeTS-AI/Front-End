cwlVersion: v1.0
class: CommandLineTool
baseCommand: PrepareDataset
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
hints:
  SoftwareRequirement:
    packages:
      PrepareDataset:
        version:
          - 0.0.1