cwlVersion: v1.0
class: CommandLineTool
baseCommand: Preprocessing
inputs:
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
  inputImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image for processing.
  maskImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -m
    doc: Input Mask for processing.
  outputImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Output Image for processing.
  histoMatch:
    type: File?
    label: NIfTI Target
    inputBinding:
      position: 1
      prefix: -hi
    doc: Match inputImage with the file provided in with this parameter.Pass the target image after '-hi'.
  hMatchBins:
    type: int?
    label: 1-1000
    inputBinding:
      position: 1
      prefix: -hb
    doc: Number of histogram bins for histogram matching.Only used for histoMatching.Defaults to 100.
  hMatchQnts:
    type: int?
    label: 1-1000
    inputBinding:
      position: 1
      prefix: -hq
    doc: Number of quantile values to match for histogram matching.Only used for histoMatching.Defaults to 40.
  zScoreNorm:
    type: boolean?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -zn
    doc: Z-Score normalization.
  zNormQuant:
    type: string?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -zq
    doc: "The Lower-Upper Quantile range to remove.Default: 5,95."
  zNormCut:
    type: string?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -zc
    doc: "The Lower-Upper Cut-off (multiple of stdDev) to remove.Default: 3,3."
  n3BiasCorr:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n3
    doc: "Runs the N3 bias correction.Optional parameters: mask or bins, spline order, filter noise level, fitting levels, max iterations, full-width-at-half-maximum."
  n4BiasCorr:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -n4
    doc: "Runs the N4 bias correction.Optional parameters: mask or bins, spline order, filter noise level, fitting levels, full-width-at-half-maximum."
  nSplineOrder:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nS
    doc: The spline order for the bias correction.Defaults to 3.
  nFilterNoise:
    type: float?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nF
    doc: The filter noise level for the bias correction.Defaults to 0.010000.
  nBiasBins:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nB
    doc: If no mask is specified, N3/N4 bias correction makes one using Otsu.This parameter specifies the number of histogram bins for Otsu.Defaults to 10.
  nFittingLevels:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nFL
    doc: The number of fitting levels to use for bias correction.Defaults to 4.
  nMaxIterations:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nMI
    doc: The maximum number of iterations for bias correction (only works for N3).Defaults to 100.
  nFullWidthHalfMaximum:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -nFWHM
    doc: Set the full-width-at-half-maximum value for bias correction.Defaults to 0.150000.
  susanSmooth:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ss
    doc: Susan smoothing of an image.
  susanSigma:
    type: float?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssS
    doc: Susan smoothing Sigma.Defaults to 0.500000.
  susanRadius:
    type: int?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssR
    doc: Susan smoothing Radius.Defaults to 1.
  susanThresh:
    type: float?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -ssT
    doc: Susan smoothing Intensity Variation Threshold.Defaults to 80.000000.
  p1p2norm:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -p12
    doc: P1-P2 normalization required for skull stripping.
  registration:
    type: string?
    label: Affine-DOF | Deformable | Rigid
    inputBinding:
      position: 1
      prefix: -reg
    doc: "The kind of registration to perform.Defaults to Affine.Can use Mask File with '-m' and multiple moving images with '-i'.For Affine, the second number defines the degrees of freedom, eg: '-ref Affine-12'."
  regFixedImg:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -rFI
    doc: The Fixed Image for the registration.Needed for registration.
  regMetrics:
    type: string?
    label: SSD | MI | NMI | NCC-AxBxC
    inputBinding:
      position: 1
      prefix: -rME
    doc: "The kind of metrics to use: SSD (Sum of Squared Differences) or MI (Mutual Information) or.NMI (Normalized Mutual Information) or NCC-AxBxC (Normalized Cross correlation with integer radius for 3D image).Defaults to NMI."
  regNoIters:
    type: string?
    label: N1,N2,N3
    inputBinding:
      position: 1
      prefix: -rNI
    doc: The number of iterations per level of multi-res.Defaults to 100,50,5.
  regInterSave:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -rIS
    doc: Whether the intermediate files are to be saved or not.Defaults to 0.
  regSegMoving:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -rSg
    doc: Whether the Moving Image(s) is a segmentation file.If 1, the 'Nearest Label' Interpolation is applied.Defaults to 0.
  regInterAffn:
    type: File?
    label: mat
    inputBinding:
      position: 1
      prefix: -rIA
    doc: The path to the affine transformation to apply to moving image.If this is present, the Affine registration step will be skipped.Also used for rigid transformation.
  regInterDefm:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -rID
    doc: The path to the deformable transformation to apply to moving image.If this is present, the Deformable registration step will be skipped.
  rescaleImage:
    type: string?
    label: Output Intensity range
    inputBinding:
      position: 1
      prefix: -rsc
    doc: The output intensity range after image rescaling.Defaults to 0.000000:1000.000000.If multiple inputs are passed (comma-separated), the rescaling is done in a cumulative manner,.i.e., stats from all images are considered for the scaling.
  debugMode:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -d
    doc: "Enabled debug mode.Default: 0."
hints:
  SoftwareRequirement:
    packages:
      Preprocessing:
        version:
          - 0.0.1