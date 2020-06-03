from setuptools import setup
from pytmt import __version__

setup(
   name='pytmt',
   version=__version__,
   description='pytmt returns ms2 tmt quantification values from Crux Percolator output',
   author='Edward Lau',
   url='https://www.maggielab.org',
   author_email='edward.lau@cuanschutz.edu',
   packages=['pytmt'],  #same as name
   install_requires=['pandas', 'pymzml', 'tqdm'], #external packages as dependencies
   entry_points={
      'console_scripts': ['pytmt=pytmt.__main__:main']
   },
)
