# -*- mode: python -*-
a = Analysis([os.path.join(HOMEPATH,'support/_mountzlib.py'), os.path.join(HOMEPATH,'support/useUnicode.py'), '/Users/farrell/python/DataPipeline/PipelineApp.py'],
             pathex=['/Users/farrell/pyinstaller-1.5.1'])
a.datas += [('/Users/farrell/python/PEATDB/Ekin/models.dict')]
pyz = PYZ(a.pure)
exe = EXE( pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name=os.path.join('dist', 'PipelineApp'),
          debug=False,
          strip=False,
          upx=True,
          console=True )
