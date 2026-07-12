"""
A module that conducts all functions related to creating files for ICORE calculations.

"""

#!/usr/bin/python3

import argparse
from vasphelper import file_manager as fm
from vasphelper import vasp_file_manager as vfm
from pathlib import Path
from typing import Any
import pandas as pd

######CONSTANTS######

CUR_DIR: Path = Path.cwd()
ICORE_DIR: Path =  fm.check_dir(CUR_DIR / 'icore')

######HELPER FUNCTIONS######

def sp_pbc_dist(point, point_list, L):

    dist_list: list[float] = []

    for pt in point_list:
        x_dist: float = point[0] - pt[0]
        y_dist: float = point[1] - pt[1]
        z_dist: float = point[2] - pt[2]

        x_dist -= L[0] * round(x_dist/L[0])
        y_dist -= L[1] * round(y_dist/L[1])
        z_dist -= L[2] * round(z_dist/L[2])

        dist_list.append(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist)

    return dist_list

def find_atoms_closest(contcar: vfm.ContcarClass, aoi: int, num_surr_atoms: int):
    contcar.clean_xyz_data()
    contcar.find_all_atom()
    requested_atoms: dict[int, tuple[str, str, int]] = {}
    atom_coords = contcar.coordinates[aoi]
    contcar.find_box_dim()
    for atom_type in contcar.types_nums.keys():
        atom_list = [key for key, val in contcar.all_atoms.items() if atom_type == val[0]+val[1]]
        if len(atom_list) > num_surr_atoms:
            atoms_coords = [contcar.coordinates[i] for i in atom_list]
            dist_list: list[float] = sp_pbc_dist(atom_coords, atoms_coords, contcar.box_dim)
            atom_dist = zip(atom_list, dist_list)
            atom_dist = sorted(atom_dist, key= lambda item: item[1])
            atom_dist = atom_dist[:num_surr_atoms]
            close_atoms, _ = zip(*atom_dist)
        else:
            close_atoms = atom_list
        requested_atoms.update({key: value for key, value in contcar.all_atoms.items() if key in close_atoms})

    return requested_atoms

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
        fm.write_text(contcar.create_icore_contcar(atom), ICORE_DIR / f'CONTCAR_{element}{tag}{pos}')
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

#causing name issues
def make_icore_directory(atom: tuple[str, str, int]) -> Path:
    element = atom[0] + atom[1]
    dir_name = f'{element}{atom[2]}'
    (ICORE_DIR / dir_name).mkdir(parents=True, exist_ok=True)
    return ICORE_DIR / dir_name


def handle_ads_icore(contcar: vfm.ContcarClass, parameter_dict: dict[str, str], *, partial: bool = False, num_surr_atoms: int | None = 7, aoi: int | None = None) -> None:
    
    if not partial:
        contcar.find_relax_atoms()
        make_icore_files(contcar, contcar.relax_atoms, parameter_dict)
    else:
        close_atoms = find_atoms_closest(contcar, aoi, num_surr_atoms) # type: ignore
        make_icore_files(contcar, close_atoms, parameter_dict)


def handle_surf_icore(contcar: vfm.ContcarClass, parameter_dict: dict[str, str], **kwargs) -> None:
    contcar.find_relax_atoms
    make_icore_files(contcar, contcar.relax_atoms, parameter_dict)


ICORE_DISPATCH: dict[str, Any] ={
    'ads': handle_ads_icore,
    'surf': handle_surf_icore
}

def run_icore_input_maker(filename: str, icore_type: str, num_ads: int, *, partial: bool | None = None, num_surr_atoms: int | None = None, aoi: int | None = None):
    
    file_path = CUR_DIR / filename

    contcar: vfm.ContcarClass = vfm.ContcarClass(file_path, num_ads)
    contcar.parse_atomic_data()

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

    icore_handler(contcar, parameter_dict, partial=partial, num_surr_atoms=num_surr_atoms, aoi=aoi)

def main():    
        
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=f"""Creates files for core level binding energy calculations.

In order to run this program, you need:
{'-'*60}
- CONTCARs for calculation in fractional coordinates
- POTCARs for each element
- INCAR from geometry calculation
- KPOINTS with correct accuracy for calculation type
- corelevel.dat containing core level data in order of POSCAR
{'-'*60}
WARNING:
This program will not change KPOINTS make sure correct KPOINTS are specified in KPOINTS file.
""", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("filename", help='Filename of the CO')
    subparsers = parser.add_subparsers(dest="icore_type", required=True,help="Choose the type of structure is contained within unit cell between Bulk 'bulk', Surface 'surf' or Surface with Adsorbates 'ads")
    surf_parser = subparsers.add_parser("surf", help="Runs file maker for surfaces")
    bulk_parser = subparsers.add_parser("bulk", help="Runs file maker for bulk")
    # add help menu to Run icoreinputmaker ads -h to see all sub arguments needed for user specifying this ICORE type.
    ads_parser = subparsers.add_parser("ads", help="Runs file maker for surfaces with adsorbates.")
    ads_subparser = ads_parser.add_subparsers(dest="mode", required=True)
    ads_all = ads_subparser.add_parser("all", help="Use all atoms")
    ads_all.add_argument("-n", "--num-ads", help="Specify the number of adsorbates species in CONTCAR. In original CONTCARs, make sure elements that occur in the surface and in the adsorbates are separated in header.", type=int, required = True)
    ads_partial = ads_subparser.add_parser("partial",  help="Use only selected atoms")
    ads_partial.add_argument("-n", "--num-ads", help="Specify the number of adsorbates species in CONTCAR. In original CONTCARs, make sure elements that occur in the surface and in the adsorbates are separated in header.", type=int, required = True)
    ads_partial.add_argument("-i", '--aoi', type = int, required = True)
    ads_partial.add_argument("-s", "--num-surr-atoms", type = int, default=7)

    args = parser.parse_args()

    if args.icore_type == 'ads':
        if args.mode == "all":
            run_icore_input_maker(
                filename=args.filename,
                icore_type = 'ads',
                num_ads = args.num_ads,
                partial = False,
            )
        else:
                run_icore_input_maker(
                filename=args.filename,
                icore_type = 'ads',
                num_ads = args.num_ads,
                partial = True,
                num_surr_atoms = args.num_surr_atoms,
                aoi = args.aoi,

            )
    elif args.icore_type == 'surf':
        run_icore_input_maker(filename=args.filename,icore_type = 'surf', num_ads = 0)


if __name__ == '__main__':
    main()
