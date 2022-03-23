""" setup module """

from setuptools import setup, find_packages
import os.path

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='pytmt',
    version='0.3.0',
    description='pytmt returns ms2 tmt quantification values from Crux Percolator output',

    long_description=long_description,
    long_description_content_type='text/markdown',

    url='https://github.com/Molecular-Proteomics/pytmt',

    author='Edward Lau, Maggie Lam',
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

    install_requires=['pandas>=1.0,<2',
                      'pymzml>=2,<3',
                      'tqdm>=4,<5'],  # external packages as dependencies
    entry_points={
        'console_scripts': [
            'pytmt=pytmt.__main__:main',
        ],
    },

    project_urls={
        'Source': 'https://github.com/Molecular-Proteomics/pytmt',
        'Maggie Lam Lab': 'http://www.maggielab.org',
    },

    data_files=[('tests',
                 [os.path.join('tests', 'data', 'percolator', 'percolator.target.mzid'),
                  os.path.join('tests', 'data', 'percolator', 'percolator.target.psms.txt'), ]),
                ],

)
