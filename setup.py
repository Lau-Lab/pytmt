""" setup module """

from setuptools import setup, find_packages
import os.path
import codecs

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


# Get version number from __init__
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name='pytmt',
    version=get_version("pytmt/__init__.py"),
    description='pytmt returns ms2 tmt quantification values from Crux Percolator output',

    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/Molecular-Proteomics/pytmt',

    author='Edward Lau',
    author_email='edward.lau@cuanschutz.edu',

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3 :: Only',
    ],

    keywords='scientific proteomics mass-spectrometry',  # Optional

    packages=find_packages(),

    python_requires='>=3.6, <4',

    install_requires=['pandas>=1,<2',
                      'pymzml>=2,<3',
                      'tqdm>=4,<5'],  # external packages as dependencies
    entry_points={
        'console_scripts': [
            'pytmt=pytmt.main:main',
        ],
    },

    project_urls={
        'Source': 'https://github.com/Molecular-Proteomics/pytmt',
        'Maggie Lam Lab': 'http://www.maggielab.irg',
    },

)
