#!/usr/bin/python3
from shutil import copy2
import pandas as pd
from pathlib import Path
import argparse
from vasphelper import file_manager as fm
from vasphelper import vasp_file_manager as vfm
from typing import Any

#####CONSTANTS######

CUR_DIR = Path.cwd()

SPLIT_DISPATCH: dict[str, list[str]]= {
    'ads': ['ads'],
    'surf': ['surf'],
    'both': ['ads', 'surf'],
    'all': ['ads', 'surf', 'overall']
}

######FUNCTIONS######

def create_list_of_cases () -> list[str]:
    data = []
    for file in CUR_DIR.iterdir():
        file = file.name
        if not (file.endswith(".vesta") or file.endswith(".png")):
            if file == "CONTCAR_CO2" or file.startswith("CONTCAR_CO2_0"):
                pass
            elif file.startswith("CONTCAR") == True:
                data.append(file.removeprefix("CONTCAR_"))
    return data

def build_dirs(case_list: list[str], dir_path: Path) -> list[Path]:
    dir_list: list[Path] = []
   
    for case in case_list:
        path = dir_path / case 
        path.mkdir()
        dir_list.append(Path(dir_path / case))
    print("Built Case Directories...")
    return dir_list

def build_rwigs_dirs(case_list: list[str], dir_path: Path, rwigs_data: list[str]) -> list[Path]:
    dir_list: list[Path] = []
    for directory in rwigs_data:
        path = dir_path / directory
        path.mkdir()
        for case in case_list:
            subpath = dir_path / directory / case 
            subpath.mkdir()
            dir_list.append(Path(dir_path / directory / case))
    print("Built Case Directories...")
    return dir_list

def get_rwigs_list(unique: set) -> pd.DataFrame:
    try:
        return pd.DataFrame(pd.read_csv(CUR_DIR / 'RWIGS_inputs.csv', index_col= 'Type'))
    except FileNotFoundError:
        print('RWIGS_inputs.csv not found.')
        while True: 
            try: 
                ratios: list[float] = [float(ratio) for ratio in input('What ratios of default RWIGS values would you like to use?\nRatios: ').split()]
                break
            except ValueError:
                print("Please enter only float values.")
    
        rwigs_list = pd.DataFrame(index = list(unique))

        for ratio in ratios:    
            col_name = f'Ratio_{str(ratio).replace('.', '_')}'
            for entry in unique:
                with open(CUR_DIR / f'POTCAR_{entry}' , 'r') as f:
                    line: int = 17
                    [next(f) for _ in range(line)]
                    rwigs_def: float = float(next(f).split()[2].replace(';', ''))
                
                rwigs_list.loc[entry, col_name] = rwigs_def * ratio
        rwigs_list.index.name = 'Type' 
        rwigs_list.to_csv(CUR_DIR / 'RWIGS_inputs.csv')

    return rwigs_list


def handle_pdos(dir_path: Path, case_list: list[str], num_ads: int, incar_parameter_dict: dict[str,str], atom_list: dict[str, str], unique_atoms: set) -> list[Path]:
    print(f'Creating files for PDOS Electon Differential Analysis...')
    rwigs_list: pd.DataFrame = get_rwigs_list(unique_atoms)
    incar_parameter_dict['LORBIT'] = '0'
    dir_list: list[Path] = build_rwigs_dirs(case_list, dir_path, rwigs_data=list(rwigs_list.columns))
    for directory in dir_list:
        rwigs_dir = directory.parent
        incar_parameter_dict['RWIGS']= ' '.join(str(rwigs_list[rwigs_dir.name].get(str(atom),'')) for atom in atom_list[directory.name])
        incar_parameter_dict['ROPT'] = f'{len(atom_list[directory.name])}*0.0005'
        vfm.populate_vasp_dirs(CUR_DIR, dir_path, directory, list(atom_list[directory.name]), incar_parameter_dict)
    return dir_list

def handle_bader(dir_path: Path, case_list: list[str], num_ads: int, incar_parameter_dict: dict[str,str], atom_list: dict[str, str], unique_atoms: set) -> list[Path]:
    print(f'Creating files for Bader Charge Analysis...')
    incar_parameter_dict['LCHARG'] = '.TRUE.'
    incar_parameter_dict['LAECHG'] = '.TRUE.'
    dir_list: list[Path] = build_dirs(case_list, dir_path)
    for directory in dir_list:
        incar_parameter_dict['ROPT'] = f'{len(atom_list[directory.name])}*0.0005'
        vfm.populate_vasp_dirs(CUR_DIR, dir_path, directory, list(atom_list[directory.name]), incar_parameter_dict)
    
    return dir_list

def handle_chg(dir_path: Path, case_list: list[str], num_ads: int, incar_parameter_dict: dict[str,str], atom_list: dict[str, str], unique_atoms: set):
    print(f'Creating files for Charge Density Analysis...')
    incar_parameter_dict['LCHARG'] = '.TRUE.'
    dir_list: list[Path] = build_dirs(case_list, dir_path)
    for directory in dir_list:
        incar_parameter_dict['ROPT'] = f'{len(atom_list[directory.name])}*0.0005'
        vfm.populate_vasp_dirs(CUR_DIR, dir_path, directory, list(atom_list[directory.name]), incar_parameter_dict)
    
    return dir_list

CALC_DISPATCH: dict[str, Any] = {
    'pdos': handle_pdos,
    'bader': handle_bader,
    'chg': handle_chg
}
    
def run_diff_input_maker(calc_type: str, split_type: str, num_ads: int):

    print(f'Progress:\n{"-"*50}\nMaking Case List...')
    case_list = create_list_of_cases()

    incar_parameter_dict: dict[str, str] = {
        'IBRION': '-1',
        'NSW': '0'
    }

    required_files: list[str] = [
        'INCAR',
        'KPOINTS',
    ]

    split_targets = SPLIT_DISPATCH[split_type]
    calc_handler = CALC_DISPATCH[calc_type]

    contcar_list = ["CONTCAR_" + case for case in case_list]

    for target in split_targets:
        dir_path = fm.check_dir(CUR_DIR / target)
        atom_dict: dict[str, list[str]] = vfm.create_contcars(num_ads, case_list, dir_path, CUR_DIR)

        unique_atoms = vfm.check_unique_atom_atom_types(atom_dict)
        required_files += [f'POTCAR_{atom}' for atom in unique_atoms]
        fm.check_files(CUR_DIR, required_files)

        dir_list = calc_handler(dir_path, case_list, num_ads, incar_parameter_dict, atom_dict, unique_atoms)
        print('Populated all input files into directories...')

        fm.remove_files(dir_path, contcar_list)

        vfm.write_calcfile(dir_path, dir_list)


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
    parser.add_argument("calc_type", choices=['bader', 'chg', 'pdos'], help="Choose the type of calculation to create files for between Bader Charge Analysis 'bader', Charge Differential Analysis 'chg' or PDOS Electron Differential Analysis 'pdos'.")
    parser.add_argument("split_type", choices=['ads', 'surf', 'both', 'all'], help="Choose the type of CONTCAR split from the surface 'surf', adsorbate 'ads',both 'both' or 'all which will produce a directory for the surface, adsorbates and overall system.")
    parser.add_argument("num_ads", help="Specify the number of adsorbates species in CONTCAR. In original CONTCARs, make sure elements that occur in the surface and in the adsorbates are separated in header.", type=int)
    args = parser.parse_args()

    run_diff_input_maker(args.calc_type, args.split_type, args.num_ads)

if __name__ == '__main__':
    main()