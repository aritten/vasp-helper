#!/usr/bin/python3

import argparse
from vasphelper import diff_input_maker
from vasphelper import icore_input_maker
from vasphelper import atom_freezer
from vasphelper import atomic_data_visualizer
from typing import Any
from pathlib import Path

def get_choice(prompt: str, choices: list) -> str:
    choice: str = input(prompt).strip()
    while True:
        if choice in choices:
            break
        choice = input(f'Choice is not included in choices please select from the above list using the corresponding numbers included in: {", ".join(choices)}:\nChoice: ').strip()
    return choice

def get_choice_w_type(prompt: str, cast: Any) -> Any:
    while True:
        try:
            choice = cast(input(prompt).strip())
            return choice
        except ValueError:
            print(f'Please enter input of type {cast}\nChoice: ')

def diff_input() -> None:
    calc_type: str = get_choice(f"""
{"-"*11}Differential File Creation Options{"-"*11}
1. Electron Distribution Analysis using PDOS
2. Bader Charge Analysis
3. Charge Differential Analysis
Choice: """, ['1', '2', '3'])
    calc_type_dict: dict[str, str] = {
        '1': 'pdos',
        '2': 'bader',
        '3': 'chg'
}
    split_type: str = get_choice("""Which files should be produced?
1. Adsorbate only
2. Surface only
3. Both surface and adsorbates                                
4. Surface, Adsorbates and Surface + Adsorbates
Choice: """, ['1', '2', '3', '4'])
    split_type_dict: dict[str, str] = {
        '1': 'ads',
        '2': 'surf',
        '3': 'both',
        '4': 'all'
    }
    #this might be best in subroutine
    num_ads: int = get_choice_w_type("Enter number of separate adsorbates in POSCAR list: ", int)
    diff_input_maker.run_diff_input_maker(calc_type_dict[calc_type], split_type_dict[split_type], num_ads)

def icore_input() -> None:

    while True: 
        filename: str = input("Enter name of the CONTCAR to use: ")
        if 'CONTCAR' in filename and (Path.cwd() / filename).exists():
            break
        print('Filename must contain CONTCAR and be in the current working directory.')

    calc_type: str = get_choice(f"""
Does the contcar contain adsorbates?
1. Yes
2. No
Choice: """, ['1', '2'])
    calc_type_dict: dict[str, str] = {
        '1': 'ads',
        '2': 'surf'
    }
    if calc_type == '1':
        num_ads: int = get_choice_w_type("Enter number of separate adsorbates in POSCAR list: ", int)
        all_atoms: str = get_choice("""Would you like directories for all relaxed atoms within the surface or only specific atoms?
1. All relaxed atoms
2. Select atoms
Choice: """, ['1', '2'])
        all_atoms_dict: dict[str, bool] = {
        '1': True,
        '2': False
        }
        if not all_atoms_dict[all_atoms]: 
            spec_atom_num: int = get_choice_w_type("Enter the atom that directories for core level binding shifts are needed: ", int)
            num_surr_atoms: int = get_choice_w_type("Enter number of atoms around that atom that need to have directories: ", int)
            icore_input_maker.run_icore_input_maker(filename, calc_type_dict[calc_type], num_ads, partial=True, num_surr_atoms=num_surr_atoms, aoi=spec_atom_num)
        else:
           icore_input_maker.run_icore_input_maker(filename, calc_type_dict[calc_type], num_ads, partial=False) 
    else:
        num_ads: int = 0
        icore_input_maker.run_icore_input_maker(filename, calc_type_dict[calc_type], num_ads)

def visualize_doscar() -> None:
    pass

def color_atoms() -> None:
    calc_type = get_choice("""
Enter calculation type:
1. Bader
2. CLBES (in progress)
3. PDOS (in progress)
Choice: """,['1', '2', '3'])
    calc_type_dict = {
        '1': 'bader',
        '2': 'clbes',
        '3': 'pdos'
    }

    mode = get_choice("""
Enter mode of coloring:
1. All
2. By atom type
3. Differential 
Choice: """,['1', '2', '3'])
    width = get_choice_w_type('Increment width: ', float)
    atomic_data_visualizer.run_atomic_data_visualizer(calc_type_dict[calc_type], mode, width = width)

def freeze_atoms() -> None:
    while True: 
        filename: str = input("Enter name of the file to freeze atoms in: ")
        if (Path.cwd() / filename).exists():
            filepath = Path.cwd() / filename
            break
        print('File must exist and be in current working directory.')
    
    freeze_type: str = get_choice('What type of freeze do you want to do?\n1. By layer\n2. By z-position\nChoice: ', ['1', '2'])
    if freeze_type == '1':
        num_layers = get_choice_w_type("How many layers are contained in the unit cell: ", int)
        relaxed_layers = get_choice_w_type("How many layers should be relaxed: ", int)
        num_ads_present = get_choice_w_type("How many adsorbates species are present: ", int)
        tolerance = get_choice_w_type("What tolerance should be used to split layers: ", float)
        atom_freezer.run_freeze_atoms(filepath, 'layer', num_layer=num_layers, relaxed_layers=relaxed_layers, num_ads=num_ads_present, tolerance=tolerance)
    else:
        zpos = get_choice_w_type("What z-position do you want to freeze atoms below? ", float)
        atom_freezer.run_freeze_atoms(filepath, 'zpos', zpos=zpos)
        
DISPATCH: dict[str, Any]= {
    "11": diff_input,
    "12": icore_input,
    "21": visualize_doscar,
    "22": color_atoms,
    "31": freeze_atoms,
    "0": exit
}
    
def handle_function(choice):
    while True:
        action = DISPATCH[choice]
        if action:
            action()
            break
        else:
            choice = input("Invalid Selection! Please select a choice from above.")

def main():
    parser = argparse.ArgumentParser(description=f"""{'-'*60}
A wrapper that handles the operation of several subroutines.
                                     
File Creation for Differential Analyses 
- Bader Charge Analysis
- Charge Differential Analysis
- PDOS Electron Distribution Analysis
    - Allows for the definition of different shell sizes around atom's core for which to calculate the number of electrons in the shell

File Creation for Core Level Binding Energy Shifts
- Final State Approximation
- Initial State Approximation (In Progress)

Freezing atoms based on layers                                 
                                    

Post Processing to Visualize Atoms in VESTA
- Bader Charge Analysis (In Progress)
    - Visualize Atoms as a function of Total Electrons
- Core Level Binding Energy (In Progress)
    - Visualize Atoms as a function of Core Level Binding Energies
- PDOS Electron Distribution Analysis (In Progress)
    - Visualize Atoms as a function of Total Electrons in user defined shell
- VESTA Viewport (In Progress)
    - Normalizes view for all CONTCAR Files
{'-'*60}
""", formatter_class=argparse.RawTextHelpFormatter)
    args = parser.parse_args()
    print(f"""{'-' * 42}
|{' '*14}VASP Helper{' '*15}|
|{' '*12}By Ariel Whitten{' '*12}|
{'-' * 42}
Enter number of the option of your choice from below to start program:\n""")    

    user_choice: str = get_choice(f"""{'-' * 15} File Creation Utilities {'-' * 16}
11. Differential Analysis
12. ICORE

{'-' * 15} Post Processing Utilities {'-' * 14}
21. DOSCAR Visualization (In Progress)
22. Color by ___ (In Progress)

{'-' * 15} CONTCAR Utilities {'-' * 14}
31. Define Relaxed and Static Layers for Single Files

0. Exit Program
Choice: """, ['11', '12', '21', '22', '31','0'])
    handle_function(user_choice)

if __name__ == "__main__":
    main()