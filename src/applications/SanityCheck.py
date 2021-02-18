import os, argparse, sys, csv, platform, subprocess, shutil
from pathlib import Path
from datetime import date

def main():
  copyrightMessage = 'Contact: software@cbica.upenn.edu/n/n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use./nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
  parser = argparse.ArgumentParser(prog='SanityCheck', formatter_class=argparse.RawTextHelpFormatter, description = 'This application performs rudimentary sanity checks the input data folder for FeTS training./n/n' + copyrightMessage)
  parser.add_argument('-inputDir', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)

  args = parser.parse_args()
  inputDir = args.inputDir

  
if __name__ == '__main__':
  if platform.system() == 'Darwin':
    sys.exit('macOS is not supported')
  else:
    main()
