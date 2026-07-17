"""
A module that visualizes data from Bader charge, CLBES and Electron Distribution Analysis.

"""

#!/usr/bin/python3

import argparse
from pathlib import Path
from typing import Any
from vasphelper import vasp_file_manager as vfm
from vasphelper import file_manager as fm
from bisect import bisect_right
from os import listdir
import numpy as np

# CONTSTANTS

CUR_DIR = Path.cwd()

COLOR_LIST = [
["0", "0", "0"],
["21", "14", "56"],
["40", "17", "89"],
["62", "15", "114"],
["85", "19", "125"],
["105", "28", "128"],
["127", "36", "129"],
["148", "43", "128"],
["171", "51", "124"],
["192", "58", "117"],
["214", "68", "108"],
["231", "82", "98"],
["244", "104", "91"],
["250", "128", "94"],
["253", "155", "106"],
["254", "179", "123"],
["253", "205", "144"],
["252", "229", "166"],
["251", "252", "191"],
["255","255", "255"]]

# FUNCTIONS

def read_bader_file(filepath: Path, case_name: str, atom_type_list: list[str]):
    bader_content = fm.read_text(filepath)
    bader_content = [i.replace('\n', '').split() for i in bader_content][2:-4]
    dtype = [('Number', 'i4'), ('Case', 'U12'), ('Element', 'U2'), ('Data', 'f4')]
    bader_values = np.array(list(zip([i for i in range(1, len(bader_content) + 1)], [case_name] * (len(bader_content) + 1), atom_type_list, [float(i[4]) for i in bader_content])), dtype= dtype)
    return bader_values

def read_clbes_file(filepath: Path, case_name: str, atom_type_list: list[str]):
    #fix this??
    clbes_content = fm.read_text(filepath)
    dtype = [('Number', 'i4'), ('Case', 'U12'), ('Element', 'U2'), ('Data', 'f4')]
    clbes_values = np.array(list(zip([int(i[0]) for i in clbes_content], [case_name] * (len(clbes_content) + 1), atom_type_list, [float(i[1]) for i in clbes_content])), dtype= dtype)
    return clbes_values

def read_pdos_file(filepath: Path, case_name: str, atom_type_list: list[str]):
    # fix this
    pdos_content = fm.read_text(filepath)
    dtype = [('Number', 'i4'), ('Case', 'U12'), ('Element', 'U2'), ('Data', 'f4')]
    pdos_values = np.array(list(zip([int(i[0]) for i in pdos_content], [case_name] * (len(pdos_content) + 1), atom_type_list, [float(i[1]) for i in pdos_content])), dtype= dtype)
    return pdos_values

def get_values_by_type(nums: list[int], types: list[str]):
    values_by_type = []
    start = 0
    end = 0
    for key, val in enumerate(nums):
        end += val
        try:
            values_by_type.extend([types[key]]*val)
        except KeyError:
            values_by_type.append([types[key]]*val)
        start += val
    return values_by_type

def get_max_value(values_list):
    #take in all atoms sorted by atom type
    max_value = {}
    for key in values_list['Element']:
        mask = values_list['Element'] == key
        element_values = list(values_list[mask]['Data'])
        max_value[str(key)] = (float(max(element_values)))
    return max_value

def make_increments(width: float, values):
    max_list = get_max_value(values)
    increments_list = {}
    for key, max_val in max_list.items():
        increments_list[key] = [max_val - int(j) * width for j in range(len(COLOR_LIST))][::-1]
    return increments_list

def color_atoms_by_value(value_list, increments: dict[str, list[float]], color_info: list[list[str]], types):
    i = 0
    temp = []
    for key in types:
        element_mask = value_list['Element'] == key
        for value in value_list[element_mask]['Data']:
            idx = bisect_right(increments[key], value) - 1
            idx = max(0, min(idx, len(COLOR_LIST) - 1))
            color_info[i][3:6] = COLOR_LIST[idx]
            i += 1
    return color_info
        
def get_case_list(prefix: str, suffix: str):
    case_list = []
    for file in listdir(CUR_DIR):
        if file.startswith(prefix):
            name = file.removeprefix(prefix).removesuffix(suffix)
            case_list.append(name)
    return case_list

def read_vesta_file(filename: str):
    vesta_content = fm.read_text(CUR_DIR / f'CONTCAR_{filename}.vesta')
    start_index = vesta_content.index("SITET\n") + 1
    end_index = vesta_content.index("  0 0 0 0 0 0\n")
    color_info = [i.replace('\n', '').split() for i in vesta_content[start_index:end_index]]
    return start_index, end_index, color_info

def write_vesta_file(filename: str, color_info: list[list[str]], start: int, end: int):
    vesta_content = fm.read_text(CUR_DIR / f'CONTCAR_{filename}.vesta')
    for i, val in enumerate(range(start,end)):
        vesta_content[val] = '    '.join(color_info[i]) + '\n'
    fm.write_text(''.join(vesta_content), CUR_DIR / f'out_CONTCAR_{filename}.vesta')

def write_bader_outfile():
    ...

def handle_pdos():
    # add method to calculate pdos
    # format of data file needs be redone
    pass

def handle_bader(width: int, case_name: str, types_list: list[str]):

    values = read_bader_file(CUR_DIR / f'ACF_{case_name}.dat', case_name, types_list)
    return values

def handle_clbes():
    # add way to calculate clbes 
    # add atom number to front of calculation for clbes
    ...

CALC_DISPATCH: dict[str, Any] = {
    'pdos': handle_pdos,
    'bader': handle_bader,
    'clbes': handle_clbes
}

CASE_LIST_DISPATCH: dict[str, list[str]] = {
    'pdos': ['', ''],
    'bader': ['ACF_', '.dat'],
    'clbes': ['', '']
}

def run_atomic_data_visualizer(calc_type: str, mode: str, *, width: float):
    # add option of all, atom, diff
    
    calc_handler = CALC_DISPATCH[calc_type]

    case_list = get_case_list(CASE_LIST_DISPATCH[calc_type][0], CASE_LIST_DISPATCH[calc_type][1])

    overall_types = {}
    accumulator = []
    for case_name in case_list:
        start, end, color_info = read_vesta_file(case_name)
        contcar = vfm.ContcarClass(CUR_DIR / f'CONTCAR_{case_name}', 0)
        contcar.parse_atomic_data()
        types_list = get_values_by_type(contcar.nums, contcar.types)
        values = calc_handler(width, case_name, types_list)
        accumulator.append(values)
        overall_types[case_name] = contcar.types

    overall_table = np.concatenate(accumulator)

    increments = make_increments(width, overall_table)
    
    for case_name in case_list:
        start, end, color_info = read_vesta_file(case_name)
        case_mask = overall_table['Case'] == case_name
        color_info = color_atoms_by_value(overall_table[case_mask], increments, color_info, overall_types[case_name])
        write_vesta_file(case_name, color_info, start, end)

def main():
    parser = argparse.ArgumentParser(description=f"""Split CONTCAR into surface and adsorbant CONTCARs and makes files for one of the follow analyses:
- Bader Charge
- Charge Density
- Electron Distribution Analysis.\n
In order to run this program, you need:
{'-'*60}
- CONTCARs for each case
- POTCARs for each element
- INCAR from geometry calculation
- KPOINTS with correct accuracy for calculation type
{'-'*60}
WARNING:
This program will not change KPOINTS make sure correct KPOINTS are specified in KPOINTS file.
""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("calc_type", choices=['bader', 'clbes', 'pdos'], help="Choose the type of visualization between Bader Charge Analysis 'bader', Core Level Binding Energy Shifts 'clbes' or PDOS Electron Differential Analysis 'pdos'.")
    parser.add_argument("mode", choices=['surf', 'atom', 'diff'], help="Choose the type of CONTCAR split from the surface 'surf', adsorbate 'ads',both 'both' or 'all which will produce a directory for the surface, adsorbates and overall system.")
    parser.add_argument("width", help="Specifies width of increments.", type=float, default=0.1)
    args = parser.parse_args()

    run_atomic_data_visualizer(calc_type=args.calc_type, mode=args.mode, width=args.width)

if __name__ == '__main__':
    main()
