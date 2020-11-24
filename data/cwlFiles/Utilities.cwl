cwlVersion: v1.0
class: CommandLineTool
baseCommand: Utilities
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
    type: string?
    label: File or Dir
    inputBinding:
      position: 1
      prefix: -i
    doc: Input Image (all CaPTk supported images) for processing.Directory to a single series DICOM only.
  maskImage:
    type: string?
    label: File or Dir
    inputBinding:
      position: 1
      prefix: -m
    doc: Input Mask (all CaPTk supported images) for processing.Directory to a single series DICOM only.
  outputImage:
    type: File?
    label: NIfTI
    inputBinding:
      position: 1
      prefix: -o
    doc: Output Image for processing.
  dicomDirectory:
    type: Directory?
    label: none
    inputBinding:
      position: 1
      prefix: -df
    doc: Absolute path of directory containing single dicom series.
  resize:
    type: int?
    label: 10-500
    inputBinding:
      position: 1
      prefix: -r
    doc: "Resize an image based on the resizing factor given.Example: -r 150 resizes inputImage by 150%.Defaults to 100, i.e., no resizing.Resampling can be done on image with 100."
  resampleResolution:
    type: string?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -rr
    doc: "[Resample] Resolution of the voxels/pixels to change to.Defaults to 1.0,1.0,1.0.Use '-rf' for a reference file."
  resampleReference:
    type: File?
    label: NIfTI image
    inputBinding:
      position: 1
      prefix: -rf
    doc: "[Resample] Reference image on which resampling is to be done.Resize value needs to be 100.Use '-ri' for resize resolution."
  resampleInterp:
    type: string?
    label: NEAREST:NEARESTLABEL:LINEAR:BSPLINE:BICUBIC
    inputBinding:
      position: 1
      prefix: -ri
    doc: "[Resample] The interpolation type to use for resampling or resizing.Defaults to LINEAR.Use NEARESTLABEL for multi-label masks."
  resampleMask:
    type: boolean?
    label: 0 or 1
    inputBinding:
      position: 1
      prefix: -rm
    doc: "[Resample] Rounds the output of the resample, useful for resampling masks.Defaults to '0'."
  sanityCheck:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -s
    doc: Do sanity check of inputImage with the file provided in with this parameter.Performs checks on size, origin & spacing.Pass the target image after '-s'.
  information:
    type: boolean?
    label: true or false
    inputBinding:
      position: 1
      prefix: -inf
    doc: Output the information in inputImage.If DICOM file is detected, the tags are written out.
  cast:
    type: string?
    label: (u)char, (u)int, (u)long, (u)longlong, float, double
    inputBinding:
      position: 1
      prefix: -c
    doc: "Change the input image type.Examples: '-c uchar', '-c float', '-c longlong'."
  uniqueVals:
    type: boolean?
    label: true or false
    inputBinding:
      position: 1
      prefix: -un
    doc: Output the unique values in the inputImage.Pass value '1' for ascending sort or '0' for no sort.Defaults to '1'.
  boundingBox:
    type: File?
    label: NIfTI Mask
    inputBinding:
      position: 1
      prefix: -b
    doc: Extracts the smallest bounding box around the mask file.With respect to inputImage.Writes to outputImage.
  boundingIso:
    type: boolean?
    label: Isotropic Box or not
    inputBinding:
      position: 1
      prefix: -bi
    doc: Whether the bounding box is Isotropic or not.Defaults to true.
  testBase:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -tb
    doc: Baseline image to compare inputImage with.
  testRadius:
    type: int?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -tr
    doc: Maximum distance away to look for a matching pixel.Defaults to 0.
  testThresh:
    type: float?
    label: 0-5
    inputBinding:
      position: 1
      prefix: -tt
    doc: Minimum threshold for pixels to be different.Defaults to 0.0.
  createMask:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -cm
    doc: "Create a binary mask out of a provided (float) thresholds.Format: -cm lower,upper.Output is 1 if value >= lower or <= upper.Defaults to 1,Max."
  changeValue:
    type: string?
    label: N.A.
    inputBinding:
      position: 1
      prefix: -cv
    doc: "Change the specified pixel/voxel value.Format: -cv oldValue1xoldValue2,newValue1xnewValue2.Can be used for multiple number of value changes.Defaults to 3,4."
  dicom2Nifti:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -d2n
    doc: If path to reference is present, then image comparison is done.Use '-i' to pass input DICOM image.Use '-o' to pass output image file.Pass a directory to '-o' if you want the JSON information.
  nifi2dicom:
    type: Directory?
    label: DICOM Reference
    inputBinding:
      position: 1
      prefix: -n2d
    doc: A reference DICOM is passed after this parameter.The header information from the DICOM reference is taken to write output.Use '-i' to pass input NIfTI image.Use '-o' to pass output DICOM directory.
  nifi2dicomDirc:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -ndD
    doc: "The direction tolerance for DICOM writing.Because NIfTI images have issues converting directions,.Ref: https://github.com/InsightSoftwareConsortium/ITK/issues/1042.this parameter can be used to override checks.Defaults to '0.000000'."
  nifi2dicomOrign:
    type: float?
    label: 0-100
    inputBinding:
      position: 1
      prefix: -ndS
    doc: "The spacing tolerance for DICOM writing.Because NIfTI images have issues converting spacings,.Ref: https://github.com/InsightSoftwareConsortium/ITK/issues/1042.this parameter can be used to override checks.Defaults to '0.000000'."
  dcmSeg:
    type: Directory?
    label: DICOM Reference
    inputBinding:
      position: 1
      prefix: -ds
    doc: A reference DICOM is passed after this parameter.The header information from the DICOM reference is taken to write output.Use '-i' to pass input NIfTI image.Use '-o' to pass output DICOM file.
  dcmSegJSON:
    type: File?
    label: JSON file for Metadata
    inputBinding:
      position: 1
      prefix: -dsJ
    doc: The extra metadata needed to generate the DICOM-Seg object.Use http://qiicr.org/dcmqi/#/seg to create it.Use '-i' to pass input NIfTI segmentation image.Use '-o' to pass output DICOM file.
  orient:
    type: string?
    label: Desired 3 letter orientation
    inputBinding:
      position: 1
      prefix: -or
    doc: The desired orientation of the image.See the following for supported orientations (use last 3 letters only):.https://itk.org/Doxygen/html/namespaceitk_1_1SpatialOrientation.html#a8240a59ae2e7cae9e3bad5a52ea3496e.Use the -bv or --bvec option to reorient an accompanying bvec file..
  bvec:
    type: File?
    label: bvec file to reorient
    inputBinding:
      position: 1
      prefix: -bv
    doc: The bvec file to reorient alongside the corresponding image.For correct output, the given file should be in the same orientation as the input image.This option can only be used alongside the -or or --orient options..
  threshAbove:
    type: float?
    label: Desired_Threshold
    inputBinding:
      position: 1
      prefix: -thA
    doc: The intensity ABOVE which pixels of the input image will be.made to OUTSIDE_VALUE (use '-tOI').Generates a grayscale image.
  threshBelow:
    type: float?
    label: Desired_Threshold
    inputBinding:
      position: 1
      prefix: -thB
    doc: The intensity BELOW which pixels of the input image will be.made to OUTSIDE_VALUE (use '-tOI').Generates a grayscale image.
  threshAnB:
    type: string?
    label: Lower_Threshold,Upper_Threshold
    inputBinding:
      position: 1
      prefix: -tAB
    doc: The intensities outside Lower and Upper thresholds will be.made to OUTSIDE_VALUE (use '-tOI').Generates a grayscale image.
  threshOtsu:
    type: boolean?
    label: 0-1
    inputBinding:
      position: 1
      prefix: -thO
    doc: Whether to do Otsu threshold.Generates a binary image which has been thresholded using Otsu.Optional mask to localize Otsu search area.
  thrshBinary:
    type: string?
    label: Lower_Threshold,Upper_Threshold
    inputBinding:
      position: 1
      prefix: -tBn
    doc: The intensity BELOW and ABOVE which pixels of the input image will be.made to OUTSIDE_VALUE (use '-tOI').Default for OUTSIDE_VALUE=0.
  threshOutIn:
    type: string?
    label: Outside_Value,Inside_Value
    inputBinding:
      position: 1
      prefix: -tOI
    doc: The values that will go inside and outside the thresholded region.Defaults to '0,1', i.e., a binary output.
  image2world:
    type: string?
    label: x,y,z
    inputBinding:
      position: 1
      prefix: -i2w
    doc: "The image coordinates that will be converted to world coordinates for the input image.Example: '-i2w 10,20,30'."
  world2image:
    type: string?
    label: i,j,k
    inputBinding:
      position: 1
      prefix: -w2i
    doc: "The world coordinates that will be converted to image coordinates for the input image.Example: '-w2i 10.5,20.6,30.2'."
  joined2extracted:
    type: boolean?
    label: 0-1
    inputBinding:
      position: 1
      prefix: -j2e
    doc: "Axis to extract is always the final axis (axis '3' for a 4D image).The '-o' parameter can be used for output: '-o /path/to/extracted_'."
  extracted2joined:
    type: float?
    label: 0-10
    inputBinding:
      position: 1
      prefix: -e2j
    doc: The spacing in the new direction.Pass the folder containing all images in '-i'.
  labelSimilarity:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -ls
    doc: Calculate similarity measures for 2 label maps.Pass the reference map after '-ls' and the comparison will be done with '-i'.For images with more than 2 labels, individual label stats are also presented.
  lSimilarityBrats:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -lsb
    doc: Calculate BraTS similarity measures for 2 brain labels.Pass the reference map after '-lsb' and the comparison will be done with '-i'.Assumed labels in image are '1,2,4' and missing labels will be populate with '0'.
  hausdorffDist:
    type: File?
    label: NIfTI Reference
    inputBinding:
      position: 1
      prefix: -hd
    doc: Calculate the Hausdorff Distance for the input image and.the one passed after '-hd'.
  collectInfo:
    type: boolean?
    label: Dir with read
    inputBinding:
      position: 1
      prefix: -co
    doc: Collects information about all images in input directory.Input directory passed using '-i'.Recursion defined using '-co 1'.Output CSV passed using '-o'.
  collectFileName:
    type: string?
    label: File pattern
    inputBinding:
      position: 1
      prefix: -cF
    doc: The file pattern to check for in every file when collecting information.Defaults to check all.
  collectFileExt:
    type: string?
    label: File extension
    inputBinding:
      position: 1
      prefix: -cFE
    doc: The file extension to check for in every file when collecting information.Defaults to check all.
  collectProperties:
    type: string?
    label: 0-2
    inputBinding:
      position: 1
      prefix: -cP
    doc: "Requested image property.0: spacings, 1: size, 2: origin.Defaults to 0.Defaults to '-cP 0,1'."
hints:
  SoftwareRequirement:
    packages:
      Utilities:
        version:
          - 0.0.1