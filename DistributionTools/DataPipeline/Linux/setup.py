from setuptools import setup
import sys,os
home=os.path.expanduser('~')

#requires we copy the source dir to the current directory first
setup(name='DataPipeline', version='1.2',
	description='DataPipeline Tool',
	author='Damien Farrell',
	author_email='farrell.damien[at]gmail.com',
	url='http://pypi.python.org/packages/source/D/DataPipeline/',
	license='GPL v3',
	packages=['DataPipeline'],
    include_package_data = True,
	package_data={'DataPipeline':['app.ico','testfiles/*.txt']},
	install_requires=['numpy>=1.1',
                      'matplotlib>=0.98.5.2',
                      'PEATDB','xlrd'],
    dependency_links = [
          "http://download.sourceforge.net/pmw/Pmw.1.3.tar.gz"],
    entry_points = { 'gui_scripts': [
                     'pipelineapp = DataPipeline.PipelineApp:main'],
                     'console_scripts': [
                     'pipelinecommand = DataPipeline.PipelineCommand:main']}
)
