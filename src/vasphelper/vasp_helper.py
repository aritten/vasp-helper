#!/usr/bin/python3

import argparse
from vasphelper import diff_input_maker
from typing import Any
from pathlib import Path

def get_choice(prompt: str, choices: list) -> str:
    choice: str = input(prompt).strip()
    while True:
        if choice in choices:
            break
        choice = input(f'Choice is not included in choices please select from {", ".join(choices)}:\n').strip()
    return choice

def get_choice_w_type(prompt: str, cast: Any) -> Any:
    while True:
        try:
            choice = cast(input(prompt).strip())
            return choice
        except ValueError:
            print(f'Please enter input of type {cast}')

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
    num_ads: int = get_choice_w_type("Enter number of separate adsorbates in POSCAR list: ", int)
    diff_input_maker.run_diff_input_maker(calc_type_dict[calc_type], split_type_dict[split_type], num_ads)

def icore_input() -> None:
    pass

def visualize_doscar() -> None:
    pass

def color_atoms() -> None:
    pass

DISPATCH: dict[str, Any]= {
    "1": diff_input,
    "2": icore_input,
    "3": visualize_doscar,
    "4": color_atoms,
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
- Bader Charge Analysis (In Progress)
- Charge Differential Analysis (In Progress)
- PDOS Electron Distribution Analysis
    - Allows for the definition of different shell sizes around atom's core for which to calculate the number of electrons in the shell

File Creation for Core Level Binding Energy Shifts
- Final State Approximation (In Progress)
- Initial State Approximation (In Progress)

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

    user_choice: str = input(f"""{'-' * 15} File Creation Utilities {'-' * 16}
1. Differential
2. ICORE (In Progress)

{'-' * 15} Post Processing Utilities {'-' * 14}
3. DOSCAR Visualization (In Progress)
4. Color by ___ (In Progress)

0. Exit Program
Choice: """)
    handle_function(user_choice)

if __name__ == "__main__":
    main()