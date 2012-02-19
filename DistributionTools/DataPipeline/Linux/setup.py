from distutils.core import setup
import sys,os
home=os.path.expanduser('~')
pythonpath=os.path.join(home,'python')

setup(name='DataPipeline', version='1.0',
	description='DataPipeline Tool',
	author='Nielsen Group, UCD',
	author_email='damien.farrell[at]ucd.ie',
	url='http://enzyme.ucd.ie/main/index.php/DataPipeline',
	license='For academic and non-profit use only. No redistribution, no alterations.',
	packages=['DataPipeline','PEATDB','PEATDB.Ekin'],
	package_dir={'':pythonpath},
	package_data={'DataPipeline':['app.ico']},
	requires=['Pmw','numpy','matplotlib'],
)
