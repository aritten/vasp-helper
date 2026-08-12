
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Any, List, Tuple, Union
from os import listdir
import argparse

from vasphelper import math_functions as mf
from vasphelper import file_manager as fm

CUR_DIR = Path.cwd()

######HELPER FUNCTIONS######

def get_case_list(prefix: str, suffix: str) -> list[str]:
    case_list = []
    for file in listdir(CUR_DIR):
        path = CUR_DIR / file
        if file.startswith(prefix) and path.is_file():
            name = file.removeprefix(prefix).removesuffix(suffix)
            case_list.append(name)
    return case_list

def clbes_convolution(energies, sigma: float, filepath: Path):

    E_pts = np.array([energy[2] for energy in energies], dtype=float)
    E, I = mf.gaussian_convolution(E_pts, sigma)

    plt.figure(figsize = (7,4))
    plt.plot(E, I)
    plt.ylabel('Intensity')
    plt.yticks(ticks=[])
    plt.xlabel('Binding Energy (eV)')
    plt.gca().invert_xaxis()
    plt.savefig(filepath)
    #plt.show()



def plot_dos(filepath: Path, og_values: np.ndarray[Any], x: np.ndarray[Any], y: np.ndarray[Any]) -> None:


    plt.figure(figsize=(10,6))
    plt.plot(og_values[:,0], og_values[:,1], 'o', label='Original Values')
    plt.plot(x, y, '-.', label='Interpolated Values')
    plt.legend()
    plt.xlabel('E - $E_F$ (eV)')
    plt.ylabel('Density of States')
    plt.xticks()

    plt.savefig(filepath)
    #plt.show()

######HANDLE FUNCTIONS######

def handle_xps(case_name: str, sigma: float, surf_energy: float, ref_energy: float, exp_energy: float):

    inputpath: Path = CUR_DIR / f'clbes_{case_name}.dat'
    outputpath: Path = CUR_DIR / f'theoretical_xps_{case_name}.png'
    clbe = mf.calculate_clbes(inputpath, surf_energy, ref_energy, exp_energy)

    clbes_convolution(clbe, sigma, outputpath)
    

def handle_dos(case_name):

    inputpath = CUR_DIR / f'DOSCAR_{case_name}'
    outputpath: Path = CUR_DIR / f'dos_{case_name}.png'
    pts = 2000
    doscar_content = fm.read_text(inputpath)
    num_atoms = int(doscar_content[0].split()[1])
    dos_len, fermi_energy = [float(i) for i in doscar_content[5].split()[2:4]]
    dos_len = int(dos_len)
    overall_dos = np.array([i.split() for i in doscar_content[6:dos_len + 6]], dtype = float)
    overall_dos[:, 0] -= fermi_energy

    x_dense, y_dense = mf.cubic_spline(overall_dos[:, 0], overall_dos[:, 1], pts)

    plot_dos(outputpath, overall_dos, x_dense, y_dense)

def handle_pdos():
    ...

######DISPATCH FUNCTION######

MODE_DISPATCH: dict[str, Any] = {
    'xps': handle_xps,
    'dos': handle_dos,
    'pdos': handle_pdos
}

CASE_LIST_DISPATCH: dict[str, list[str]] = {
    'pdos': ['DOSCAR_', ''],
    'dos': ['DOSCAR_', ''],
    'xps': ['clbes_', '.dat']
}

def run_data_plotter(mode, **calc_kwargs):

    handler = MODE_DISPATCH[mode]

    case_list = get_case_list(CASE_LIST_DISPATCH[mode][0], CASE_LIST_DISPATCH[mode][1])
    if mode == 'xps' and calc_kwargs['ref_energy'] == None:
            calc_kwargs['ref_energy'] = mf.calculate_ref_energy(CUR_DIR, case_list)

    print("Graphs made for:")
    accumulator= []
    for case_name in case_list:
        print(case_name)
        handler(case_name, **calc_kwargs)


######MAIN FUNCTION######

def main():
    parser = argparse.ArgumentParser(description=f"""Graph:
    - Theoretical XPS spectra
    - Density of States
    - Partial Density of States\n
    In order to run this program, you need:
    {'-'*60}
    - Input file for data for each calculation type:
        - Core Level Binding Energy: clbes_*.dat (accessed using batch submission script "clbes_job.sh")
        - Electron Distribution Analysis: DOSCARs for each case
    {'-'*60}
    """, formatter_class=argparse.RawTextHelpFormatter)
    mode_subparser = parser.add_subparsers(dest = "mode", help="Choose the type of graph to make: 'atom', 'all' or 'diff'.")
    xps_parser = mode_subparser.add_parser('xps')
    xps_parser.add_argument("surf_energy", type=float)
    xps_parser.add_argument("width", help="Specifies width of increments.", type=float, default=0.1, nargs="?")
    xps_parser.add_argument("-e", "--exp_energy", type=float, default=0)
    xps_parser.add_argument("-r", "--ref_energy", default = None, type=float)
    dos_parser = mode_subparser.add_parser('dos')
    pdos_parser = mode_subparser.add_parser('pdos')
    pdos_parser.add_argument('-o', '--orbitals', default='t', type = str.lower, nargs = '+', choices=['s', 'p', 'd', 'f', 't', 'T', 'all'])
    
    args = parser.parse_args()

    if args.mode == 'xps':
        run_data_plotter(mode=args.mode, sigma=args.width, surf_energy=args.surf_energy, ref_energy=args.ref_energy, exp_energy= args.exp_energy)
    elif args.mode == 'dos':
        run_data_plotter(mode=args.mode)


if __name__ == '__main__':
    main()