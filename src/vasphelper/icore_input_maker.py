#!/usr/bin/python3


######################################HOW TO###########################################
# This program will preform ICORE calculation file creation automatically for both    #
# surface and adsorbate calculations.                                                 #
# Necessary Files:                                                                    #
# - POTCAR_(extention of atom type) for each type                                     #
# - POSCAR file with Selective Dynamics                                               #
# - INCAR file without ICORELEVEL DATA but IBRION = -1 and NSW = 0                    #
# - KPOINTS                                                                           #
# - corelevel.dat containing core level data in order of POSCAR                       #
# Surface calculations are specified by                                               #
# python3 ICOREMAKER.py (filename) -surf                                              #
# This make files for each relaxed atoms and enter directories into a CalcFile.dat    #
# list allowing for easy use in batch submission.                                     #
# Adsorbate calculations are specified by:                                            #
# python3 ICOREMAKER.py (filename) -ads (number of atom of interest)                  #
# (number of adsorbates)                                                              #
# Optional tags:                                                                      #
# -num (number of atoms create files for)                                             #
# -tol (specifies tolerance for atoms near edge of periodic boundary to test across   #
# boundaries)                                                                         #
# How to run example file:                                                            #
# Surface:                                                                            #
# python3 ICOREMAKER.py CONTCAR_CH4_2 -surf                                           #
# Adsorbate:                                                                          #
# python3 ICOREMAKER.py CONTCAR_CH4_2 -ads 241 5 -num 7 -tol 0.3                      #
#######################################################################################

import sys
import argparse
from vasphelper import file_manager as fm
from vasphelper import vasp_file_manager as vfm
from pathlib import Path
from typing import Any
from itertools import accumulate
import pandas as pd

#####CONSTANTS######

CUR_DIR: Path = Path.cwd()
ICORE_DIR: Path =  fm.check_dir(CUR_DIR / 'icore')

def find_atoms_closest(atom, relax_atom, coords, num_closest, tolerance): 
        
    #intialize lists
    dist_list = [100]*num_closest
    close = [0]*num_closest

    #set tolerance
    low_tol = tolerance
    high_tol = 1 - tolerance

    for i in range(len(relax_atom)):
        atom_select = [float(x) for x in coords[atom]]
        atom_test = [float(x) for x in coords[relax_atom[i]]]
        dist = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
        #find atoms that are within same periodic cell distances
        if dist <= dist_list[num_closest - 1]:
            dist_list[num_closest - 1] = dist
            close[num_closest - 1] = relax_atom[i]
            dist_atom = zip(dist_list, close)
            dist_atom = list(dist_atom)
            temp_res = sorted(dist_atom, key = lambda x: x[0])
            dist_list, close = zip(*temp_res)
            dist_list = list(dist_list)
            close = list(close)
        #test atoms near periodic edges for their distance across boundary conditions
        elif (atom_select[0] > high_tol or atom_select[0] < low_tol) or (atom_select[1] > high_tol or atom_select[1] < low_tol):
            if (atom_test[0] > high_tol or atom_test[0] < low_tol) and (atom_test[1] > high_tol or atom_test[1] < low_tol):
                atom_test_x1 = atom_test[0] + 1
                dist1 = (( atom_select[0] - atom_test_x1 ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                atom_test_x2 = atom_test[0] - 1
                dist2 = (( atom_select[0] - atom_test_x2 ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                atom_test_y1 = atom_test[1] + 1
                dist3 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test_y1 ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                atom_test_y2 = atom_test[1] - 1
                dist4 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test_y2 ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                dist_order = [dist1, dist2, dist3, dist4]
                sorted_dist_order = sorted(dist_order)
                dist = sorted_dist_order[0]
                if dist <= dist_list[num_closest - 1]:
                    dist_list[num_closest - 1] = dist
                    close[num_closest - 1] = relax_atom[i]
                    dist_atom = zip(dist_list, close)
                    dist_atom = list(dist_atom)
                    temp_res = sorted(dist_atom, key = lambda x: x[0])
                    dist_list, close = zip(*temp_res)
                    dist_list = list(dist_list)
                    close = list(close)
            elif (atom_test[0] > high_tol or atom_test[0] < low_tol):
                atom_test[0] = atom_test[0] + 1
                dist1 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                atom_test[0] = atom_test[0] - 1
                dist2 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                dist_order = [dist1, dist2]
                sorted_dist_order = sorted(dist_order)
                dist = sorted_dist_order[0]
                if dist <= dist_list[num_closest - 1]:
                    dist_list[num_closest - 1] = dist
                    close[num_closest - 1] = relax_atom[i]
                    dist_atom = zip(dist_list, close)
                    dist_atom = list(dist_atom)
                    temp_res = sorted(dist_atom, key = lambda x: x[0])
                    dist_list, close = zip(*temp_res)
                    dist_list = list(dist_list)
                    close = list(close)
            elif (atom_test[1] > high_tol or atom_test[1] < low_tol):
                atom_test[1] = atom_test[1] + 1
                dist1 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                atom_test[1] = atom_test[1] - 1
                dist2 = (( atom_select[0] - atom_test[0] ) ** 2 + ( atom_select[1] - atom_test[1] ) ** 2 + ( atom_select[2] - atom_test[2] ) ** 2) ** (1/2)
                dist_order = [dist1, dist2]
                sorted_dist_order = sorted(dist_order)
                dist = sorted_dist_order[0]
                if dist <= dist_list[num_closest - 1]:
                    dist_list[num_closest - 1] = dist
                    close[num_closest - 1] = relax_atom[i]
                    dist_atom = zip(dist_list, close)
                    dist_atom = list(dist_atom)
                    temp_res = sorted(dist_atom, key = lambda x: x[0])
                    dist_list, close = zip(*temp_res)
                    dist_list = list(dist_list)
                    close = list(close)

    return (close)

def read_icore_data() -> pd.DataFrame:
    #get data from corelevel.dat
    icore_data: pd.DataFrame = pd.DataFrame(pd.read_csv(CUR_DIR / 'corelevel.csv', index_col= 'Type'))
    return icore_data

def make_icore_files(contcar: vfm.ContcarClass, atom_dict: dict[int, tuple[str, str, int]], parameter_dict: dict[str, str]) -> None:
    icore_data = read_icore_data()
    dir_list: list[Path] = []
    for atom in atom_dict.values():
        element, tag, pos = atom
        dir_name: Path = make_icore_directory(atom)
        #wont work for adsorbates yet need to add adsorbate tag to types_nums
        fm.write_text(contcar.create_icore_contcar(atom), ICORE_DIR / f'CONTCAR_{element}{tag}{pos}')
        #CLNT not correct
        parameter_dict.update({
                'ROPT': f'{len(contcar.icore_types)}*0.0005',
                'ICORELEVEL' : '2',
                'CLNT': str(contcar.icore_index),
                'CLN': icore_data['CLN'][element],
                'CLL': icore_data['CLL'][element],
                'CLZ': '1'
            })
        vfm.populate_vasp_dirs(CUR_DIR, ICORE_DIR / f'CONTCAR_{element}{tag}{pos}', dir_name, contcar.icore_types, parameter_dict)
        dir_list.append(dir_name)
    fm.remove_files(ICORE_DIR, [f'CONTCAR_{name.name}' for name in dir_list])
    vfm.write_calcfile(ICORE_DIR, dir_list)


def make_icore_directory(atom: tuple[str, str, int]) -> Path:
    element = atom[0] + atom[1]
    dir_name = f'{element}{atom[2]}'
    (ICORE_DIR / dir_name).mkdir(parents=True, exist_ok=True)
    return ICORE_DIR / dir_name


def handle_ads_icore(contcar: vfm.ContcarClass, parameter_dict: dict[str, str], *, all_atoms: bool = True, num_surr_atoms: int | None = 7, aoi: int | None = None, tolerance: float | None = 0.3) -> None:
    
    total_atoms: int = sum(contcar.nums)

    ads_atoms: int = sum(contcar.nums[-contcar.num_ads:])
    ads_dict: dict[int, tuple[str, str, int]] = {}
    surf_dict: dict[int, tuple[str, str, int]] = {}
    surf_total = total_atoms - ads_atoms - 1
    for key, value in contcar.relax_atoms.items():
        if key > surf_total:
            contcar.relax_atoms[key] = (value[0], '_ads', value[2])
    
    if all_atoms:
        make_icore_files(contcar, contcar.relax_atoms, parameter_dict)
    else:
        #separate the atoms by type
        contcar.clean_xyz_data()
        close_atoms: dict[int, tuple[str, str, int]] = {}
        for atom_type in contcar.types_nums.keys():
            print(atom_type)
            atoms_of_atom_type = [key for key, val in contcar.relax_atoms.items() if atom_type == val[0]+val[1]]
            atoms_surf_with_ads = find_atoms_closest(aoi, atoms_of_atom_type, contcar.coordinates, num_surr_atoms, tolerance)
            close_atoms.update({key: value for key, value in contcar.relax_atoms.items() if key in atoms_surf_with_ads})
        make_icore_files(contcar, close_atoms, parameter_dict)
    
    print("Directories made and names entered into the CalcFile.dat")
    print("Use batch submission script if you need WAVECAR add directly to submission script\nto reduce storage usage and allow for easier removal after running job.")


def handle_surf_icore(contcar: vfm.ContcarClass, parameter_dict: dict[str, str], **kwargs) -> None:
    print("Directories made and names entered into the CalcFile.dat")
    make_icore_files(contcar, contcar.relax_atoms, parameter_dict)
    print("Use batch submission script if you need WAVECAR add directly to submission script\nto reduce storage usage and allow for easier removal after running job.")


ICORE_DISPATCH: dict[str, Any] ={
    'ads': handle_ads_icore,
    'surf': handle_surf_icore
}

def run_icore_input_maker(filename: str, icore_type: str, num_ads: int, *, all_atoms: bool | None = None, num_surr_atoms: int | None = None, aoi: int | None = None, tolerance: float | None = None):
    
    file_path = CUR_DIR / filename

    contcar: vfm.ContcarClass = vfm.ContcarClass(file_path, num_ads)
    contcar.parse_atomic_data()
    relax_atoms = contcar.find_relax_atoms()

    icore_handler = ICORE_DISPATCH[icore_type]

    required_files: list[str] = [
        'INCAR',
        'KPOINTS',
        'corelevel.csv'
    ]

    parameter_dict: dict[str, str] = {
        'IBRION': '-1',
        'NSW': '0'
    }

    unique_atoms = vfm.check_unique_atom_atom_types(contcar.types)
    required_files += ['POTCAR_'+ unique for unique in unique_atoms]
    fm.check_files(CUR_DIR, required_files)

    icore_handler(contcar, parameter_dict, all_atoms=all_atoms, num_surr_atoms=num_surr_atoms, aoi=aoi, tolerance=tolerance)

def main():    
    
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=f"""######################################HOW TO###########################################
# This program will preform ICORE calculation file creation automatically for both    #
# surface and adsorbate calculations.                                                 #
# Necessary Files:                                                                    #
# - POTCAR_(extention of atom type) for each type                                     #
# - POSCAR file with Selective Dynamics                                               #
# - INCAR file without ICORELEVEL DATA but IBRION = -1 and NSW = 0                    #
# - KPOINTS                                                                           #
# - corelevel.dat containing core level data in order of POSCAR                       #
# Surface calculations are specified by                                               #
# python3 ICOREMAKER.py (filename) -surf                                              #
# This make files for each relaxed atoms and enter directories into a CalcFile.dat    #
# list allowing for easy use in batch submission.                                     #
# Adsorbate calculations are specified by:                                            #
# python3 ICOREMAKER.py (filename) -ads (number of atom of interest)                  #
# (number of adsorbates)                                                              #
# Optional tags:                                                                      #
# -num (number of atoms create files for)                                             #
# -tol (specifies tolerance for atoms near edge of periodic boundary to test across   #
# boundaries)                                                                         #
# How to run example file:                                                            #
# Surface:                                                                            #
# python3 ICOREMAKER.py CONTCAR_CH4_2 -surf                                           #
# Adsorbate:                                                                          #
# python3 ICOREMAKER.py CONTCAR_CH4_2 -ads 241 5 -num 7 -tol 0.3                      #
#######################################################################################
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("filename")
    subparsers = parser.add_subparsers(dest="icore_type", required=True,help="Choose the type of structure is contained within unit cell between Bulk 'bulk', Surface 'surf' or Surface with Adsorbates 'ads")
    surf_parser = subparsers.add_parser("surf", help="Runs file maker for surfaces")
    bulk_parser = subparsers.add_parser("bulk", help="Runs file maker for bulk")
    ads_parser = subparsers.add_parser("ads", help="Runs file maker for surfaces with adsorbates. Run icoreinputmaker ads -h to see all sub arguments needed for user specifying this ICORE type.")
    ads_subparser = ads_parser.add_subparsers(dest="mode", required=True)
    ads_all = ads_subparser.add_parser("all", help="Use all atoms")
    ads_all.add_argument("-n", "--num-ads", help="Specify the number of adsorbates species in CONTCAR. In original CONTCARs, make sure elements that occur in the surface and in the adsorbates are separated in header.", type=int, required = True)
    ads_partial = ads_subparser.add_parser("partial",  help="Use only selected atoms")
    ads_partial.add_argument("-n", "--num-ads", help="Specify the number of adsorbates species in CONTCAR. In original CONTCARs, make sure elements that occur in the surface and in the adsorbates are separated in header.", type=int, required = True)
    ads_partial.add_argument("-i", '--aoi', type = int, required = True)
    ads_partial.add_argument("-s", "--num-surr-atoms", type = int, default=7)
    ads_partial.add_argument('-t', '--tolerance', type = float, default = 0.3)

    args = parser.parse_args()

    if args.icore_type == 'ads':
        if args.mode == "all":
            run_icore_input_maker(
                filename=args.filename,
                icore_type = 'ads',
                num_ads = args.num_ads,
                all_atoms = True,
            )
        else:
                run_icore_input_maker(
                filename=args.filename,
                icore_type = 'ads',
                num_ads = args.num_ads,
                all_atoms = False,
                num_surr_atoms = args.num_surr_atoms,
                aoi = args.aoi,
                tolerance= args.tolerance
            )
    elif args.icore_type == 'surf':
        run_icore_input_maker(filename=args.filename,icore_type = 'surf', num_ads = 0)


if __name__ == '__main__':
    main()
