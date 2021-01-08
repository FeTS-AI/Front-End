import os, argparse, sys, pkg_resources
from pathlib import Path
from datetime import date
import csv

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
          if ((temp == "patientid") or (temp == "subjectid") or (temp == "subject") or (temp == "subid")):
            headers.append('ID')
          elif ((temp == "t1gd") or (temp == "t1ce") or (temp == "t1post")):
            headers.append("T1GD")
          elif ((temp == "t1") or (temp == "t1pre")):
            headers.append("T1")
          elif ((temp == "t2")):
            headers.append("T2")
          elif ((temp == "t2flair") or (temp == "flair") or (temp == "fl") or ('fl' in temp) or ('t2fl' in temp)):
            headers.append("FLAIR")

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

def main():
  copyrightMessage = 'Contact: software@cbica.upenn.edu\n\n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
  parser = argparse.ArgumentParser(prog='PrepareDataset', formatter_class=argparse.RawTextHelpFormatter, description = "This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately.\n\n" + copyrightMessage)
  parser.add_argument('-inputCSV', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)
  parser.add_argument('-outputDir', type=str, help = 'The output file to write the results', required=True)

  args = parser.parse_args()
  
  GetCSVContents(args.inputCSV)

if __name__ == '__main__':
  main()
