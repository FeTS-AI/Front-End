import os, argparse, sys, platform
from datetime import date
import SimpleITK as sitk
import numpy as np

def imageSanityCheck(targetImageFile, inputImageFile) -> bool:
  '''
  This function does sanity checking of 2 images
  '''
  targetImage = sitk.ReadImage(targetImageFile)
  inputImage = sitk.ReadImage(inputImageFile)

  if targetImage.GetDimension() != 3:
    print('Dimension for target image, \'' + targetImageFile + '\' is not 3', file = sys.stderr)
    return False

  if inputImage.GetDimension() != 3:
    print('Dimension for input image, \'' + inputImageFile + '\' is not 3', file = sys.stderr)
    return False

  commonMessage = ' mismatch for target image, \'' + targetImageFile + '\' and input image, \'' + inputImageFile + '\''
  problemsIn = ''
  returnTrue = True

  if targetImage.GetSize() != inputImage.GetSize():
    if not problemsIn:
      problemsIn += 'Size'
    else:
      problemsIn += ', Size'
    returnTrue = False
    
  if targetImage.GetOrigin() != inputImage.GetOrigin():
    if not problemsIn:
      problemsIn += 'Origin'
    else:
      problemsIn += ', Origin'
    returnTrue = False

  if targetImage.GetSpacing() != inputImage.GetSpacing():
    if not problemsIn:
      problemsIn += 'Spacing'
    else:
      problemsIn += ', Spacing'
    returnTrue = False

  if returnTrue:
    return True
  else:
    print(problemsIn + commonMessage, file = sys.stderr)
    return False

def main():
  # copyrightMessage = 'Contact: software@cbica.upenn.edu/n/n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use./nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
  # parser = argparse.ArgumentParser(prog='SanityCheck', formatter_class=argparse.RawTextHelpFormatter, description = 'This application performs rudimentary sanity checks the input data folder for FeTS training./n/n' + copyrightMessage)
  # parser.add_argument('-inputDir', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)

  # args = parser.parse_args()
  # inputDir = args.inputDir

  inputDir = "C:/Users/sarth/Downloads/fets_sanity_check"

  if not os.path.isdir(inputDir):
    sys.exit('The specified inputDir is not present, please try again')

  errorMessage = 'Subject_ID,Label_File,Label_Value,Number_of_Voxels\n'
  numberOfProblematicCases = 0

  # initialize modality dict
  files_to_check = {
    'T1': '_t1.nii.gz',
    'T1CE': '_t1ce.nii.gz',
    'T2': '_t2.nii.gz',
    'FL': '_flair.nii.gz',
    'MASK': '_final_seg.nii.gz'
  }

  label_values_expected = np.array([0,1,2,4]) # initialize label array
  
  for dirs in os.listdir(inputDir):
    if dirs != 'logs': # don't perform sanity check for the 'logs' folder
      currentSubjectDir = os.path.join(inputDir, dirs)
      if os.path.isdir(currentSubjectDir): # for detected subject dir
        filesInDir = os.listdir(currentSubjectDir) # get all files in each directory
        files_for_subject = {}
        for i in range(len(filesInDir)):
          for modality in files_to_check: # check all modalities
            if filesInDir[i].endswith(files_to_check[modality]): # if modality detected, populate subject dict
              files_for_subject[modality] = os.path.abspath(os.path.join(currentSubjectDir, filesInDir[i]))
        
        if len(files_for_subject) != 5: # if all modalities are not present, add exit statement
          numberOfProblematicCases += 1
          errorMessage += dirs + ',Not_all_modalities_present,N.A.,N.A.\n'

        first, *rest = files_for_subject.items() # split the dict
        for i in range(0, len(rest)):
          if not(imageSanityCheck(first[1], rest[i][1])): # image sanity check
            numberOfProblematicCases += 1
            errorMessage += dirs + ',Image_dimension_mismatch_between_' + first[0] + '_and_' + rest[i][0] + ',N.A.,N.A.\n'
        
        if 'MASK' in files_for_subject:
          currentLabelFile = files_for_subject['MASK']
          mask_array = sitk.GetArrayFromImage(sitk.ReadImage(currentLabelFile))
          unique, counts = np.unique(mask_array, return_counts=True) # get unique elements and their counts
          if not(np.array_equal(unique,label_values_expected)): # this is for the case where the label contains numbers other than 0,1,2,4
            numberOfProblematicCases += 1
            for j in range(0,len(unique)): # iterate over a range to get counts easier
              if not(unique[j] in label_values_expected):
                errorMessage += dirs + ',' + currentLabelFile + ',' + str(unique[j]) + ',' + str(counts[j]) + '\n'
        else:
          numberOfProblematicCases += 1
          errorMessage += dirs + ',Label_absent,N.A.,N.A.\n'

  if numberOfProblematicCases > 0:
    print('There were problematic cases found in the dataset. Please see the following:')
    sys.exit(errorMessage)
  else:
    print('Congratulations, all subjects are fine and ready to train!')

if __name__ == '__main__':
  if platform.system() == 'Darwin':
    sys.exit('macOS is not supported')
  else:
    main()
