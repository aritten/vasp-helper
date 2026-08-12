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
    - [Differential Input File Maker](#differential-input-file-maker)
        - [Bader Charge](#bader-charge)
        - [Charge Density](#charge-density)
        - [Electron Distribution](#electron-distribution)
    - [ICORE Input File Maker](#icore-input-file-maker)
        - [Bulk Unit Cell](#bulk-unit-cell)
        - [Surface](#surface)
        - [Surface with Adsorbates](#surface-with-adsorbates)
    - [Atom Freezer](#atom-freezer)
        - [Freeze by Layer](#freeze-by-layer)
        - [Freeze by z-position](#freeze-by-z-position)
    - [Data Plotter](#data-plotter)
        - [Plot XPS](#plot-xps)
        - [Plot DOS](#plot-dos)
    - [Atomic Data Visualizer](#atomic-data-visualizer)
        - [Visualize CLBE](#visualize-clbes)
        - [Visualize PDOS](#visualize-pdos)
        - [Visualize Bader Charge](#visualize-bader-charge)
- [Future Plans](#future-plans)

## Why?

I needed a tool to process large numbers of files quickly and reduce human error when creating files. There are tools that do some of what is accomplished in this tool but often they are not built to handle large numbers of files. Additionally, it is hard to find tools that have command line functionality so this tool can be called in a script to automate many processes.

## Installation

`vasp-helper` can be installed on a computer with Python 3.10 and beyond.

To install, use
```
git clone https://github.com/aritten/vasp-helper
```

Then, use pip to install vasp-helper and its dependencies.

## Dependencies
* `pandas`

## Usage

This tool can be used by calling `vasphelper` and navigating menu options to get all necessary information to run the program or indivdual functions can be called by entering their name into command line with all the required tags.

### Differential Input File Maker

Run the program in the desired file location where CONTCARs are present. This can be run for just surface `surf`, for just adsorbates `ads`, for both surface and ads `both` or for surface, ads and overall `all`.

**Required Files**
- CONTCARs for each case
- POTCARs for each element with element name matching name in CONTCAR (ex: for La -> POTCAR_La)
- INCAR from geometry calculation
- KPOINTS with correct accuracy for calculation type

#### <ins>Bader Charge</ins>

Makes files for Bader Charge Analysis.
```
diffinputmaker bader {ads, surf, both, all} {Number of adsorbate species in the unit cell}
```

#### <ins>Charge Density</ins>

Makes files for charge density analysis.
```
diffinputmaker chg {ads, surf, both, all} {Number of adsorbate species in the unit cell}
```

#### <ins>Electron Distribution</ins>

Makes files for electron distribution calculations.
```
diffinputmaker pdos {ads, surf, both, all} {Number of adsorbate species in the unit cell}
```

### ICORE Input File Maker

Run in folder with desired position file. It makes folders for all relaxed atoms within the unit cell.

**Required Files**
- CONTCAR from geomentry calculation
- POTCARs for each element with element name matching name in CONTCAR (ex: for La -> POTCAR_La)
- INCAR from geometry calculation
- KPOINTS with correct accuracy for calculation type

#### <ins>Bulk Unit Cell</ins>

Makes files for all atoms within a bulk unit cell.

```
icoreinputmaker {filename} bulk
```

#### <ins>Surface</ins>

Makes files for all relaxed atoms within a surface.
```
icoreinputmaker {filename} surf
```

#### <ins>Surface with Adsorbates</ins>

There are two modes that can be used if you have a suface with adsorbates. The first mode makes input files for all relaxed atoms. The second mode makes input files for only atoms surrounding a choosen atom and its nearest neighbors of each species.

##### <ins>All</ins>

Makes files for all relaxed atoms in the surface. 
```
icoreinputmaker {filename} ads all {Number of adsorbate species}
```
##### <ins>Partial</ins>

Makes files for only atoms around a specified atom.
```
icoreinputmaker {filename} ads partial {Number of adsorbate species} {Number of Choosen Atom}
```
**Optional Tags**

`-s` `--num-surr-atoms` - Specifies the number of atoms from each specie surrounding the choosen atom that files should be made for. Default is 7

### Atom Freezer

Run in folder where desired position file is located. 

**Required Files**
- Position file in VASP format

#### <ins>Freeze by Layer</ins>
```
atomfreezer {filename} layer {Number of Layers in Surface} {Number of relaxed layers}
```
**Optional Tags**

`-n` or `--num_ads` - Specifies the number of adsorbate species present in the surface

`-t` or `--tolerance` - Specifies tolerance to determine the presence of a layer. Default value is 0.01.

#### <ins>Freeze by z-position</ins>
```
atomfreezer {filename} zpos {z-position to freeze surface at}
```
### Data Plotter

Run in folder where desired position file is located. 

#### <ins>Plot XPS</ins>

Solves the following equation:
$$
CLBE = (E_{n-1}- E_{surface}) - (E_{n-1}_{ref} - E_{surf}) + E_{experimental}
$$
where $E_{n-1}$ is the final approximation energy, $E_{surf}$ is the surface energy after geometry optimization, $E_{n-1}_{ref} is the reference final approximation energy and $E_{experimental}$ is the experimental energy shift. 

Convolutes spectra using Gaussians broadened by user inputted value of `{Gaussian Broadening Width}`

**Required Files**
- Core Level Binding Shift Output File 'cbles_*.dat' (Created from clbes_job.sh in repositiory)

```
dataplotter xps {Surface Energy} {Gaussian Broadening Width}
```
**Optional Tags**

`-e` or `--exp_energy` - Specifies the experimental energy shift. Default is 0.

`-r` or `--ref_energy` - Specifies the reference energy for core level binding energy equation.

#### <ins>Plot DOS</ins>

Plots the Density of States and cubic spline of Density of States.

**Required Files**
- DOSCAR files for each case

```
dataplotter dos
```

### Atomic Data Visualizer

Run in folder where desired position file is located. 

#### <ins>Visualize CLBES</ins>

Solves the following equation:
$$
CLBE = (E_{n-1}- E_{surface}) - (E_{n-1}_{ref} - E_{surf}) + E_{experimental}
$$
where $E_{n-1}$ is the final approximation energy, $E_{surf}$ is the surface energy after geometry optimization, $E_{n-1}_{ref} is the reference final approximation energy and $E_{experimental}$ is the experimental energy shift. 

Colors atoms in .vesta file by binding energy values

**Required Files**
- Core Level Binding Shift Output File 'cbles_*.dat' for each Case (Created from clbes_job.sh in repositiory)
- CONTCAR from VASP for each Case
- CONTCAR.vesta for each Case
- corelevel.csv with orbital information for each atom (example in tests)

```
atomicdatavisualizer all {Increment Width} clbes {Surface Energy}
```
**Optional Tags**

`-e` or `--exp_energy` - Specifies the experimental energy shift. Default is 0.

`-r` or `--ref_energy` - Specifies the reference energy for core level binding energy equation.

#### <ins>Visualize PDOS</ins>

Visualizes the number of electrons per choosen orbital for atoms in unit cell.

**Required Files**
- DOSCAR for each case
- CONTCAR from VASP for each case
- CONTCAR.vesta for each case

```
atomicdatavisualizer all {Increment Width} pdos
```
**Optional Tags**

`-o` or `--orbitals` - Specifies orbitals to display in .vesta file. Default is 't' for total. Options 's', 'p', 'd', 'f' and 't'.

#### <ins>Visualize Bader Charge</ins>

Visualizes the charge for atoms in unit cell.

**Required Files**
- 'ACF.dat' for each case
- CONTCAR from VASP for each case
- CONTCAR.vesta for each case

```
atomicdatavisualizer all {Increment Width} bader
```

## Future Plans

This program is a work in progress and the following functions will be added:
* Add file creation functionality
    * Create large throughput geometry optimziation  
* Add functionality to allow for post-processing steps
    * Plot PDOS
    * Allow for atoms be colored for Atomic Data Visualizer by differences between case and reference case and by atom type

Besides functionality, the code will be continuously improved upon to create the best user experience possible.