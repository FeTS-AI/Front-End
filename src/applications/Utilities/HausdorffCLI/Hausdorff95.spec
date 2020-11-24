# -*- mode: python ; coding: utf-8 -*-

import sys
sys.setrecursionlimit(5000)

import os

block_cipher = None


a = Analysis(['Hausdorff95.py'],
             pathex=[os.getcwd()],
             hiddenimports=['numpy.core._dtype_ctypes', 'pkg_resources.py2_warn'],
             binaries=[],
             datas=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='Hausdorff95',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
