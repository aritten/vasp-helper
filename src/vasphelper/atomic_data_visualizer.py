"""
A module that visualizes data from Bader charge, CLBES and Electron Distribution Analysis.

"""

#!/usr/bin/python3

import argparse
from pathlib import Path
from typing import Any, List, Tuple, Union
from numpy.typing import NDArray
from vasphelper import vasp_file_manager as vfm
from vasphelper import file_manager as fm
from vasphelper import math_functions as mf
from bisect import bisect_right
from os import listdir
import numpy as np
import time


######CONSTANTS######

CUR_DIR = Path.cwd()

COLOR_NEUTRAL = ["220", "220", "220"]

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

######FILE FUNCTIONS######

def read_bader_file(filepath: Path, 
                    case_name: str, 
                    atom_type_list: list[str]) -> NDArray[Any]:
    bader_content = fm.read_text(filepath)
    bader_content = [i.replace('\n', '').split() for i in bader_content][2:-4]
    dtype = [('Number', 'i4'), ('Case', 'U20'), ('Element', 'U2'), ('Data', 'f4')]
    bader_values = np.array(list(zip([i for i in range(1, len(bader_content) + 1)], [case_name] * len(bader_content), atom_type_list, [float(i[4]) for i in bader_content])), dtype= dtype)
    return bader_values

def read_outfile(filepath: Path, 
                 case_name: str, 
                 atom_type_list: list[str], 
                 data_col: tuple[int, int], 
                 dtype: list[tuple[str, str]]) -> NDArray[Any]:
    content = fm.read_text(filepath)
    data: list[tuple] = []
    for idx, line in enumerate((content), start= 1):
        values = line.strip().split()
        values = [float(v) for v in values[data_col[0]:data_col[1]]]
        data.append((line[0], idx, line[2], *values))
    return np.array(data, dtype= dtype)

def read_vesta_file(filename: str) -> tuple[int, int, list[list[str]]]:
    vesta_content = fm.read_text(CUR_DIR / f'CONTCAR_{filename}.vesta')
    start_index = vesta_content.index("SITET\n") + 1
    end_index = vesta_content.index("  0 0 0 0 0 0\n")
    color_info = [i.replace('\n', '').split() for i in vesta_content[start_index:end_index]]
    return start_index, end_index, color_info

def write_vesta_file(filename: str, 
                     orb: str, 
                     calc_type: str, 
                     output_dir: Path, 
                     color_info: list[list[str]], 
                     start: int, 
                     end: int) -> None:
    vesta_content = fm.read_text(CUR_DIR / f'CONTCAR_{filename}.vesta')
    for i, val in enumerate(range(start,end)):
        vesta_content[val] = '    '.join(color_info[i]) + '\n'
    fm.write_text(''.join(vesta_content), output_dir / f'out_{calc_type}_CONTCAR_{orb}_{filename}.vesta')

def write_outfile(filepath: Path, 
                  output_table: NDArray[Any], 
                  calc_type: str, 
                  case_list: list[str]) -> None:
    contents = [f'{calc_type.upper()} Data']
    labels = '  '.join(list(output_table.dtype.names or []))
    for case_name in case_list:
        contents.append(case_name)
        contents.append(labels)
        mask = output_table['Case'] == case_name
        atom_data = []
        for entry in output_table[mask]:
            atom_data.append(' '.join([str(val) for val in list(entry)[1:]]))
        contents.append('\n'.join(atom_data))
    
    fm.write_text('\n'.join(contents), filepath)

######HELPER FUNCTIONS######

def get_values_by_type(nums: list[int], types: list[str]) -> list[str]:
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

def color_atoms_by_value(value_list: NDArray[Any], 
                         col: str, 
                         increments: dict[str, list[float]], 
                         color_info: list[list[str]], 
                         types: list[str]) -> list[list[str]]:
    i = 0
    for key in types:
        element_mask = value_list['Element'] == key
        for value in value_list[element_mask][col]:
            if value != -5:
                idx = bisect_right(increments[key], value) - 1
                idx = max(0, min(idx, len(COLOR_LIST) - 1))
                color_info[i][3:6] = COLOR_LIST[idx]
            else:
                color_info[i][3:6] = COLOR_NEUTRAL
            i += 1
    return color_info
   
def get_case_list(prefix: str, suffix: str) -> list[str]:
    case_list = []
    for file in listdir(CUR_DIR):
        path = CUR_DIR / file
        if file.startswith(prefix) and path.is_file():
            name = file.removeprefix(prefix).removesuffix(suffix)
            case_list.append(name)
    return case_list

######HANDLE FUNCTIONS######

def handle_pdos(contcar: vfm.ContcarClass,
                width: float, 
                case_name: str, 
                types_list: list[str]) -> NDArray[Any]:
    # pass orbital info to this function
    num, electron_list = mf.calculate_pdos(CUR_DIR / f'DOSCAR_{case_name}')
    dtype = []
    orbtial_list: list[tuple]= [('s', 'f4'), ('p', 'f4'), ('d', 'f4'), ('f', 'f4')]
    dtype = [('Case', 'U20'), ('Number', 'i4'), ('Element', 'U2'), *[orbtial_list[orb] for orb in range(len(electron_list[0]) - 1)], ('t', 'f4')]
    data = []
    for idx, (elem_type, orb) in enumerate(zip(types_list, electron_list), start= 1):
        data.append((case_name, idx, elem_type, *orb))
    pdos_values = np.array(data, dtype= dtype)
    return pdos_values


def handle_bader(contcar: vfm.ContcarClass,
                 width: int, 
                 case_name: str, 
                 types_list: list[str]) -> NDArray[Any]:

    bader_values = read_bader_file(CUR_DIR / f'ACF_{case_name}.dat', case_name, types_list)
    return bader_values

def handle_clbes(contcar: vfm.ContcarClass, 
                 width: float, 
                 case_name: str, 
                 types_list: list[str], 
                 surf_energy: float,  
                 ref_energy: float,
                 exp_energy: float) -> NDArray[Any]:

    clbes_list = mf.calculate_clbes(CUR_DIR / f'clbes_{case_name}.dat', surf_energy, ref_energy, exp_energy)
    clbes_data = []
    for atom in clbes_list:
        start = 0
        for atom_type in contcar.types_nums:
            if atom_type == atom[0]:
                clbes_data.append((case_name, start + int(atom[1]), atom_type, atom[2]))
            else:  
                start += contcar.types_nums[atom_type][1]
    data = []
    total_atom = 1
    for atom_type in contcar.types_nums:
        for i in range(contcar.types_nums[atom_type][1]):
            data.append((case_name, total_atom, atom_type, '-5'))
            total_atom += 1
    
    for atom in clbes_data:
        data[atom[1] - 1] = atom

    dtype = [('Case', 'U20'), ('Number', 'i4'), ('Element', 'U2'), ('Data', 'f4')]
    clbes_values = np.array(data, dtype=dtype)
    return clbes_values    


######DISPATCH FUNCTION######

CALC_DISPATCH: dict[str, Any] = {
    'pdos': handle_pdos,
    'bader': handle_bader,
    'clbes': handle_clbes
}

CASE_LIST_DISPATCH: dict[str, list[str]] = {
    'pdos': ['DOSCAR_', ''],
    'bader': ['ACF_', '.dat'],
    'clbes': ['clbes_', '.dat']
}

def run_atomic_data_visualizer(calc_type: str, 
                               mode: str, 
                               *, 
                               width: float, 
                               col_list: list[str] = ['t'], 
                               **calc_kwargs) -> None:
    # add option of all, atom, diff
    # add option to choose how to color based on orbital

    calc_handler = CALC_DISPATCH[calc_type]
    case_list = get_case_list(CASE_LIST_DISPATCH[calc_type][0], CASE_LIST_DISPATCH[calc_type][1])

    if calc_type == 'clbes' and calc_kwargs['ref_energy'] == None:
            calc_kwargs['ref_energy'] = mf.calculate_ref_energy(CUR_DIR, case_list)
        
    overall_types = {}
    accumulator = []
    print(f'VESTA files colored by {calc_type.upper()} created for:')
    for case_name in case_list:
        print(case_name)
        contcar = vfm.ContcarClass(CUR_DIR / f'CONTCAR_{case_name}', 0)
        contcar.parse_atomic_data()
        types_list = get_values_by_type(contcar.nums, contcar.types)
        values = calc_handler(contcar, width, case_name, types_list, **calc_kwargs)
        accumulator.append(values)
        overall_types[case_name] = list(dict.fromkeys(contcar.types))

    overall_table = np.concatenate(accumulator)
    output_dir: Path = fm.check_dir(CUR_DIR / f'{calc_type}_out')
    
    for col in col_list:
        increments = mf.make_increments(width, overall_table, col, 'Element', len(COLOR_LIST))
        for case_name in case_list:
            start, end, color_info = read_vesta_file(case_name)
            case_mask = overall_table['Case'] == case_name
            color_info = color_atoms_by_value(overall_table[case_mask], col, increments, color_info, overall_types[case_name])
            write_vesta_file(case_name, col, calc_type, output_dir, color_info, start, end)

    write_outfile(output_dir / f'{calc_type}_outfile.dat', overall_table, calc_type, case_list)


######MAIN FUNCTION######

def main() -> None:
    parser = argparse.ArgumentParser(description=f"""Visualize data in VESTA files from:
- Bader Charge
- Core Level Binding Energy
- Electron Distribution Analysis.\n
In order to run this program, you need:
{'-'*60}
- CONTCARs for each case
- .vesta for each case
- Input file for data for each calculation type:
    - Bader Charge Analysis: ACF_*.dat
    - Core Level Binding Energy: clbes_*.dat (accessed using batch submission script "clbes_job.sh")
    - Electron Distribution Analysis: DOSCARs for each case
{'-'*60}
""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "mode",
        choices=["all", "atom", "diff"],
        help="Choose the type of CONTCAR split: 'atom', 'all' or 'diff'."
    )
    parser.add_argument("width", help="Specifies width of increments.", type=float, default=0.1, nargs="?")
    calc_subparsers = parser.add_subparsers(dest="calc_type", help="Choose the type of visualization between Bader Charge Analysis 'bader', Core Level Binding Energy Shifts 'clbes' or PDOS Electron Differential Analysis 'pdos'.")
    pdos_subparser = calc_subparsers.add_parser("pdos")
    pdos_subparser.add_argument('-o', '--orbitals', default='t', type = str.lower, nargs = '+', choices=['s', 'p', 'd', 'f', 't', 'T'])
    bader_subparser = calc_subparsers.add_parser("bader")
    clbes_subparser = calc_subparsers.add_parser("clbes")
    clbes_subparser.add_argument("surf_energy", type=float)
    clbes_subparser.add_argument("-e", "--exp_energy", type=float, default=0)
    clbes_subparser.add_argument("-r", "--ref_energy", default = None, type=float)
    
    args = parser.parse_args()

    if args.calc_type == 'pdos':
        run_atomic_data_visualizer(calc_type=args.calc_type, 
                                   mode=args.mode, 
                                   width=args.width, 
                                   col_list=args.orbitals)
    elif args.calc_type == 'clbes':
        run_atomic_data_visualizer(calc_type=args.calc_type, 
                                  mode=args.mode, 
                                  width=args.width,
                                  col_list= ['Data'],
                                  surf_energy=args.surf_energy,
                                  ref_energy=args.ref_energy,
                                  exp_energy=args.exp_energy)
    else:
        run_atomic_data_visualizer(calc_type=args.calc_type, 
                            mode=args.mode, 
                            width=args.width,
                            col_list= ['Data'])

if __name__ == '__main__':
    main()
