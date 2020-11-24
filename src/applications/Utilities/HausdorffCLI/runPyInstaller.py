import os
import PyInstaller
from shutil import copyfile
os.chdir('./')
print("Creating exe")

os.system("pyinstaller --onefile Hausdorff95.spec")

print("All Done")