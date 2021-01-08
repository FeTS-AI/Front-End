import os, argparse, sys, pkg_resources
from pathlib import Path
from datetime import date


def main():
    copyrightMessage = 'Contact: software@cbica.upenn.edu\n\n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
    parser = argparse.ArgumentParser(prog='PrepareDataset', formatter_class=argparse.RawTextHelpFormatter, description = "This application calls the BraTSPipeline for all input images and stores the final and intermediate files separately.\n\n" + copyrightMessage)
    parser.add_argument('-inputCSV', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)
    parser.add_argument('-outputDir', type=str, help = 'The output file to write the results', required=True)

if __name__ == '__main__':
    main()
