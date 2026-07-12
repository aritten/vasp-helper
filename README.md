## Description

`vasp-helper` is a program that handles file creation for VASP workflows. It is used to handle large numbers of files to greatly reduce introduction of human error and speed up the process. 

It handles files creation for the following processes:

* **Bader Charge Analysis** 
* **Charge Differential Analysis**
* **PDOS Electron Differential Analysis**
* **Core Level Binding Energy Shifts Calculations**

Additionally, it can be used to freeze atoms in CONTCAR by layer and z-position.

## Contents

- [Why?](#why)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Future Plans](#future-plans)

## Why?

I needed a tool to process large numbers of files quickly and reduce human error when creating files. There are tools that do some of what is accomplished in this tool but often they are not built to handle large numbers of files. Additionally, it is hard to find tools that have command line functionality so this tool can be called in a script to automate many processes.

## Installation

Vasp-Helper can be installed on a computer with Python 3.10 and beyond.

To install, use
```
git clone https://github.com/aritten/vasp-helper
```

Then, use pip to install vasp-helper and its dependencies.

## Dependencies
* `pandas`

## Usage

This tool can be used by calling `vasphelper` and navigating menu options to get all necessary information to run the program or indivdual functions can be called by entering their name into command line with all the required tags.

**Differential Input File Maker**

`diffinputmaker {}

**ICORE Input File Maker**

`icoreinputmaker

### Atom Freezer

For a run that freezes by layer
```
atomfreezer {filename} layer {Number of Layers in Surface} {Number of relaxed layers}
```
**Optional Tags**

`-n` or `--num_ads` - Specifies the number of adsorbate species present in the surface

`-t` or `--tolerance` - Specifies tolerance to determine the presence of a layer. Default value is 0.01.

## Future Plans

This program is a work in progress and the following functions will be added:
* Add file creation functionality
    * Create large throughput geometry optimziation  
* Add functionality to allow for post-processing steps
    * Convolute core level binding shift data to produce theoretical XPS spectra
    * Color atoms by core level binding energy, number of electrons and position within super cell
    * Plot output of DOS and PDOS

Besides functionality, the code will be continuously improved upon to create the best user experience possible.