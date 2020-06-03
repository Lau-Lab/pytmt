PyTMT returns ms2 tandem mass tag quantification values from Crux/Percolator output.

## Getting Started


#### Requirements

Install Python 3.6+ and pip. See instructions on Python website for specific instructions for your operating system.

Install pytmt from PyPI. We recommend using a virtual environment.
		
	$ pip install pytmt


#### Running
	
Launch pytmt as a Python module (Usage/Help)
	
	$ python -m pytmt

Console entry point:

    $ pytmt
    
Example command: 
	
	$ python -m pytmt /path/to/mzml/ /path/to/percolator/ -o /path/to/output/
	

### Prerequisites

pytmt is tested in Python 3.6, 3.7, and 3.8 and requires the following packages:

```
pandas>=1.0.4
pymzml==2.4.6
tqdm>=4.46.0
```

## Contributing
Please contact us if you wish to contribute, and submit pull requests to us.


## Authors
* **Edward Lau, PhD** - *Code/design* - [ed-lau](https://github.com/ed-lau)
* **Maggie Lam, PhD** - *Code/design* - [Maggie-Lam](https://github.com/Maggie-Lam)


## License
This project is licensed under the MIT License - see the LICENSE file for details
