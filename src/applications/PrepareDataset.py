import os, argparse, sys, pkg_resources
from pathlib import Path
from datetime import date


def main():
    copyrightMessage = 'Contact: software@cbica.upenn.edu\n\n' + 'This program is NOT FDA/CE approved and NOT intended for clinical use.\nCopyright (c) ' + str(date.today().year) + ' University of Pennsylvania. All rights reserved.' 
    parser = argparse.ArgumentParser(prog='LabelFusion', formatter_class=argparse.RawTextHelpFormatter, description = "Fusion of different labels together.\n\n" + copyrightMessage)
    parser.add_argument('-inputs', type=str, help = 'The absolute, comma-separated paths of labels that need to be fused', required=True)
    parser.add_argument('-classes', type=str, help = 'The expected labels; for example, for BraTS, this should be \'0,1,2,4\'', required=True)
    parser.add_argument('-method', type=str, help = 'The method to apply; currently available: STAPLE | ITKVoting | MajorityVoting | SIMPLE', required=True)
    parser.add_argument('-output', type=str, help = 'The output file to write the results', required=True)

if __name__ == '__main__':
    main()
