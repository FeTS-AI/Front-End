import os, argparse, sys, csv, platform, subprocess, shutil
from pathlib import Path
from datetime import date

def GetCSVContents(filename):
  '''
  Read filename and return a list of dictionaries that have the csv contents
  '''
  with open(filename, 'r') as csvfile:
    datareader = csv.reader(csvfile)
    
    parserHeader = True
    headers = [] # save headers
    csvContents = [] # csv contents
    for row in datareader:

      if parserHeader: # parser headers first

        for col in row:
          temp = col.lower() # convert to lower case
          if ((temp == 'patientid') or (temp == 'subjectid') or (temp == 'subject') or (temp == 'subid')):
            headers.append('ID')
          elif ((temp == 't1gd') or (temp == 't1ce') or (temp == 't1post')):
            headers.append('T1GD')
          elif ((temp == 't1') or (temp == 't1pre')):
            headers.append('T1')
          elif ((temp == 't2')):
            headers.append('T2')
          elif ((temp == 't2flair') or (temp == 'flair') or (temp == 'fl') or ('fl' in temp) or ('t2fl' in temp)):
            headers.append('FLAIR')

        parserHeader = False

      else:
        if len(headers) != 5:
          sys.exit('All required headers were not found in CSV. Please ensure the following are present: \'PatientID,T1,T1GD,T2,T2FLAIR\'')

        col_counter = 0
        currentRow = {}
        for col in row: # iterate through columns
          if ' ' in col:
            sys.exit('Please ensure that there are no spaces in the file paths.')
          else:
            currentRow[headers[col_counter]] = col # populate header with specific identifiers
          col_counter += 1

        csvContents.append(currentRow) # populate csv rows
  
  return csvContents

def copyFilesToCorrectLocation(interimOutputDir, finalSubjectOutputDir, subjectID):
  '''
  This function copies the intermediate files and final outputs to correct location and if these are absent, returns a bool flag stating that brats pipeline needs to run again
  '''

  # copy files to correct location for inference and training
  runBratsPipeline = False
  output_t1c_brain_file_inter = os.path.join(interimOutputDir, "brain_T1CE.nii.gz")
  output_t1c_brain_file_final = os.path.join(finalSubjectOutputDir, subjectID + "_brain_t1ce.nii.gz")
  if os.path.exists(output_t1c_brain_file_inter):
    shutil.copyfile(output_t1c_brain_file_inter, output_t1c_brain_file_final)
  else:      
    output_t1c_brain_file_inter = os.path.join(interimOutputDir, "brain_T1GD.nii.gz")
    if os.path.exists(output_t1c_brain_file_inter):
      shutil.copyfile(output_t1c_brain_file_inter, output_t1c_brain_file_final)
    else:
      runBratsPipeline = True

    
  output_t1_brain_file_inter = os.path.join(interimOutputDir, "brain_T1.nii.gz")
  output_t1_brain_file_final = os.path.join(finalSubjectOutputDir, subjectID + "_brain_t1.nii.gz")
  if os.path.exists(output_t1_brain_file_inter):
    shutil.copyfile(output_t1_brain_file_inter, output_t1_brain_file_final)
  else:
    runBratsPipeline = True
    
  output_t2_brain_file_inter = os.path.join(interimOutputDir, "brain_T2.nii.gz")
  output_t2_brain_file_final = os.path.join(finalSubjectOutputDir, subjectID + "_brain_t2.nii.gz")
  if os.path.exists(output_t2_brain_file_inter):
    shutil.copyfile(output_t2_brain_file_inter, output_t2_brain_file_final)
  else:
    runBratsPipeline = True
    
  output_fl_brain_file_inter = os.path.join(interimOutputDir, "brain_FL.nii.gz")
  output_fl_brain_file_final = os.path.join(finalSubjectOutputDir, subjectID + "_brain_flair.nii.gz")
  if os.path.exists(output_fl_brain_file_inter):
    shutil.copyfile(output_fl_brain_file_inter, output_fl_brain_file_final)
  else:
    runBratsPipeline = True

  return runBratsPipeline

def main():
  copyrightMessage = 'Contact: software@cbica.upenn.edu/n/n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use./nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
  parser = argparse.ArgumentParser(prog='PrepareDataset', formatter_class=argparse.RawTextHelpFormatter, description = 'This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately./n/n' + copyrightMessage)
  parser.add_argument('-inputCSV', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)
  parser.add_argument('-outputDir', type=str, help = 'The output file to write the results', required=True)

  args = parser.parse_args()
  outputDir_qc = os.path.normpath(args.outputDir + '/DataForQC')
  outputDir_final = os.path.normpath(args.outputDir + '/DataForFeTS')

  Path(args.outputDir).mkdir(parents=True, exist_ok=True)
  Path(outputDir_qc).mkdir(parents=True, exist_ok=True)
  Path(outputDir_final).mkdir(parents=True, exist_ok=True)

  bratsPipeline_exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'BraTSPipeline')
  if platform.system() == 'Windows':
    bratsPipeline_exe += '.exe'

  csvContents = GetCSVContents(args.inputCSV)

  for row in csvContents:
    interimOutputDir_actual = os.path.join(outputDir_qc, row['ID'])
    finalSubjectOutputDir_actual = os.path.join(outputDir_final, row['ID'])
    runBratsPipeline = copyFilesToCorrectLocation(interimOutputDir_actual, finalSubjectOutputDir_actual, row['ID'])

    if runBratsPipeline:
      command = bratsPipeline_exe + ' -t1 ' + row['T1'] + ' -t1c ' + row['T1GD'] + ' -t2 ' + row['T2'] + ' -fl ' + row['FLAIR'] + ' -o ' + interimOutputDir + ' -s 1'
      print('Command: ', command)
      subprocess.Popen(command, shell=True).wait()
    
    if copyFilesToCorrectLocation(interimOutputDir, finalSubjectOutputDir, row['ID']):
      print('BraTSPipeline failed for subject \'', row['ID'], file = sys.stderr)

if __name__ == '__main__':
  if platform.system() == 'Darwin':
    sys.exit('macOS is not supported')
  else:
    main()
