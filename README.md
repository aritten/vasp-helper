# Description

Vasp Helper is a program that handles file creation for a variety of differential analyses:

* **Bader Charge Analysis** 
* **Charge Differential Analysis**
* **PDOS Electron Differential Analysis**

## Installation

Vasp-Helper can be installed on a computer with Python 3.10 and beyond. 

### Dependencies
* `pandas`

To install, use
```
git clone https://github.com/aritten/vasp-helper
```

Then, use pip to install vasp-helper and its dependencies.

## Future Plans

This program is a work in progress and the following functions will be added
* Add new modes of file creation
    * Create files to calculate the core level binding energy shifts
    * Freeze atoms based on position within the super cell
* Add functionality to allow for post-processing steps
    * Convolute core level binding shift data to produce theoretical XPS spectra
    * Color atoms by core level binding energy, number of electrons and position within super cell
    * Plot output of DOS and PDOS

Besides functionality, the code will be continuously improved upon to create the best user experience possible.