cwlVersion: v1.0
class: CommandLineTool
baseCommand: DeepMedic
inputs:
  inputImages:
    type: string
    label: NIfTI files
    inputBinding:
      position: 1
      prefix: -i
    doc: "Input images provided in a comma-separated manner.Should be in the same order as trained model.Example: '-i /path/t1.nii.gz,/path/t1c.nii.gz'."
  modelDir:
    type: Directory
    label: none
    inputBinding:
      position: 1
      prefix: -md
    doc: The trained model to use.See examples in '${CaPTk_installDir}/data/deepMedic/saved_models'.
  output:
    type: Directory
    label: none
    inputBinding:
      position: 1
      prefix: -o
    doc: The output directory.
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
  mask:
    type: File?
    label: none
    inputBinding:
      position: 1
      prefix: -m
    doc: The Optional input mask file..This is needed for normalization only.
  zScoreNorm:
    type: boolean?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -zn
    doc: Z-Score normalization.Set to '0' if you are passing normalized images.
  zNormQuant:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -zq
    doc: "The Lower-Upper Quantile range to remove.Default: 5.000000,95.000000."
  zNormCut:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -zc
    doc: "The Lower-Upper Cut-off (multiple of stdDev) to remove.Default: 3.000000,3.000000."
  resizeResolution:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -rr
    doc: "[Resample] Isotropic resampling resolution to change to.Defaults to 1.000000.If '0' is passed, no resampling is done.."
  debugModel:
    type: boolean?
    label: 0-1
    inputBinding:
      position: 1
      prefix: -d
    doc: Enable/disable debug mode for extra information on console.
  Logger:
    type: string?
    label: log file which user has write access to
    inputBinding:
      position: 1
      prefix: -L
    doc: Full path to log file to store console outputs.By default, only console output is generated.
hints:
  SoftwareRequirement:
    packages:
      DeepMedic:
        version:
          - 0.0.1.nonRelease