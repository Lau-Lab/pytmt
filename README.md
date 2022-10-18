# PyTMT v0.3.1

PyTMT returns ms2 tandem mass tag quantification values from Crux/Percolator output.

## Getting Started


#### Requirements

Install Python 3.7+ and pip. See instructions on Python website for specific instructions for your operating system.

Install pytmt from PyPI. We recommend using a virtual environment.
		
	$ pip install pytmt


#### Running
	
Launch pytmt as a Python module (Usage/Help)
	
	$ python -m pytmt

Console entry point:

    $ pytmt
    
Example command: 
	
	$ python -m pytmt /path/to/mzml/ /path/to/percolator/psms.txt -o /path/to/output/
	
To test installation of test data files and download two test mzml files from ProteomeXchange:

    $ pip install tox
    $ tox

To run pytmt on the test data files and print the file to Desktop:
    
    $ pytmt tests/data/mzml tests/data/percolator/percolator.target.psms.txt -o ~/Desktop/pytmt2
    

### Prerequisites

pytmt is tested in Python 3.7, and 3.8 and requires the following packages:

```
pandas>=1.0.4
pymzml==2.4.6
tqdm>=4.46.0
```

## TMT tags 

Below are the accurate masses used to find the TMT tags signals, and the channels corresponding
to each option in the `-m` argument:

| Tag mass      | Tag name | 0   | 2 | 6 | 10 | 11 | 16 |  18 |
| ----------- | ----------- | --- | --- | --- | --- | --- | --- | --- |
| 126.127726  | 126 |x | x|x |x |x | x| x|
| 127.124761  | 127N| | x| | x| x| x| x|
| 127.131081  | 127C| | | x| x|x | x| x|
| 128.128116  | 128N| | | | x|x | x| x|
| 128.134436  | 128C| | |x | x|x | x| x|
| 129.131471  | 129N| | | x| x|x | x| x|
| 129.137790  |  129C| | | | x|x | x| x|
| 130.134825  | 130N| | | | x| x| x| x|
| 130.141145  | 130C| | | x| x|x | x| x|
| 131.138180  | 131N| | | x| x|x | x| x|
| 131.144500  | 131C (11-plex) | | | | |x | x| x|
| 132.141535  | 132N (Pro)| | | | | | x| x|
| 132.147855  |  132C (Pro)| | | | | | x| x|
| 133.144890  |  133N (Pro)| | | | | | x| x|
| 133.141210  | 133C (Pro)| | | | | | x| x|
| 134.148245  | 134N (Pro)| | | | | | x| x|
| 134,154565  | 134C (Pro-18)| | | | | | | x|
| 135.151600   | 135N (Pro-18)| | | | | | | x|

## All Options

```angular2html
python -m pytmt -h
usage: __main__.py [-h] [-u] [-q QVALUE] [-m MULTIPLEX] [-p PRECISION] [-o OUT] [-v] [-c CONTAM] [-n] mzml id

pytmt returns ms2 tmt quantification valuesfrom Percolator output and perform contaminationcorrection

positional arguments:
  mzml                  path to folder containing mzml files
  id                    path to percolator target psms output file

optional arguments:
  -h, --help            show this help message and exit
  -u, --unique          quantify unique peptides only
  -q QVALUE, --qvalue QVALUE
                        quantify peptides with q value below this threshold [default: 1.0]
  -m MULTIPLEX, --multiplex MULTIPLEX
                        TMT-plex (0, 2, 6, 10, 11, 16, 18) [default:10]
  -p PRECISION, --precision PRECISION
                        ms2 spectrum mass shift tolerance in ppm [default: 10]
  -o OUT, --out OUT     name of the output directory [default: tmt_out]
  -v, --version         show program's version number and exit
  -c CONTAM, --contam CONTAM
                        Path to contaminant matrix csv file. Leave blank to get tmt output without correction
  -n, --nnls            uses non-negative least square for contamination correction

```

## Contributing
Please contact us if you wish to contribute, and submit pull requests to us.


## Authors
* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam](https://github.com/Maggie-Lam)


## License
This project is licensed under the MIT License - see the LICENSE file for details
