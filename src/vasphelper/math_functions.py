#!/usr/bin/python3

from typing import Any, List, Tuple, Union
from numpy.typing import NDArray
import numpy as np
from scipy.integrate import trapezoid
from scipy import interpolate as inter
from pathlib import Path
from vasphelper import file_manager as fm

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

def get_max_value(values_list: NDArray[Any], col: str, sortby: str) -> dict[str, float]:
    # take in all atoms sorted by atom type
    max_value = {}
    for key in values_list[sortby]:
        mask = values_list[sortby] == key
        element_values = list(values_list[mask][col])
        max_value[str(key)] = (round(float(max(element_values)), 5))
    return max_value

def make_increments(width: float, values: NDArray[Any], col: str, sortby: str, num_increments: int) -> dict[str, list[float]]:
    max_list = get_max_value(values, col, sortby)
    increments_list = {}
    for key, max_val in max_list.items():
        increments_list[key] = [max_val - int(j) * width for j in range(num_increments)][::-1]
    return increments_list

def cubic_spline(x: NDArray[Any], y: NDArray, num_pts: int):
    x_dense = np.linspace(x.min(), 0, num_pts)
    #utilized cubic because industry standard and it runs faster since it is highly optimized to run the fortran code
    cubic = inter.make_interp_spline(x, y, k=3)
    y_dense = np.clip(cubic(x_dense), 0.0, None)

    return x_dense, y_dense

def inter_to_zero(x: NDArray[Any], y: NDArray[Any], num_pts: int) -> float:
    y_dense, x_dense = cubic_spline(x, y, num_pts)
    zero_value = trapezoid(y_dense, x_dense)

    return zero_value

def gaussian_convolution(E_pts: np.ndarray, sigma_eV: float):

    I_pts = np.ones_like(E_pts)

    E_grid = np.linspace(E_pts.min(), E_pts.max(), 600)
    diff = E_grid[:, None] - E_pts[None, :]
    kernel = np.exp(-(diff**2) / (2 * sigma_eV**2))
    I_smooth = kernel @ I_pts

    return E_grid, I_smooth

def calculate_orbitals(doscar_content: list[str], 
                       num_atoms: int, 
                       dos_len: int, 
                       pts: int) -> tuple[list[list[float]], float]:
    overall_from_pdos = 0
    total_electron_list = []
    for atom in range(1, num_atoms + 1):
        pdos_content = np.array([i.split() for i in doscar_content[1:int(dos_len) + 1]], dtype = float)
        pdos_content[:, 0] -= float(doscar_content[0].split()[3])
        electron_by_orbital = []
        for orb in range(1, len(pdos_content[atom])):
            electron_by_orbital.append(float(inter_to_zero(pdos_content[:,0], pdos_content[:,orb], pts)))
        total_electrons: float = sum(electron_by_orbital)
        overall_from_pdos += total_electrons
        electron_by_orbital.append(total_electrons)
        total_electron_list.append(electron_by_orbital)
        del doscar_content[0:int(dos_len) + 1]
    return total_electron_list, overall_from_pdos

def calculate_overall(doscar_content: list[str], 
                      dos_len: int, 
                      fermi_energy: float, 
                      pts: int) -> float:

    overall_dos = np.array([i.split() for i in doscar_content[6:dos_len + 6]], dtype = float)
    overall_dos[:, 0] -= fermi_energy
    overall_electrons = inter_to_zero(overall_dos[:, 0], overall_dos[:, 1], pts)

    return overall_electrons

def calculate_pdos(filepath: Path) -> tuple[int, list[list[float]]]:
    pts = 2000
    doscar_content = fm.read_text(filepath)
    num_atoms = int(doscar_content[0].split()[1])
    dos_len, fermi_energy = [float(i) for i in doscar_content[5].split()[2:4]]
    dos_len = int(dos_len)
    #make this route somewhere
    total_electrons = round(calculate_overall(doscar_content, dos_len, fermi_energy, pts), 0)
    total_electron_list, all_electron_count = calculate_orbitals(doscar_content[int(dos_len) + 6 :], num_atoms, dos_len, pts)

    return num_atoms, total_electron_list

def calculate_clbes(filepath: Path, 
                    surf_energy: float,  
                    e_ref: float,
                    exp_energy: float) -> list[list[Union[str, int, float]]]:
    clbes_content = fm.read_text(filepath)
    clbes_content = [i.replace('\n', '').split() for i in clbes_content]

    clbes_values = [[atom[0], round((float(atom[1]) - surf_energy) - (e_ref - surf_energy) + exp_energy,2)] for atom in clbes_content]

    for i, atom in enumerate(clbes_values):
        clbes_values[i] = ["".join([c for c in atom[0] if not c.isdigit()]), int("".join([c for c in atom[0] if c.isdigit()])), atom[1]]
    
    return clbes_values

def calculate_ref_energy(path: Path, case_list: list[str]) -> float:
    content = []
    ref_energy = 0
    for case_name in case_list:
       content = [i.replace('\n', '').split() for i in fm.read_text(path / f'clbes_{case_name}.dat')]
       ref_energy = min([float(i[1]) for i in content])
       #ref_energy = fm.read_text(CUR_DIR / f'clbes_{case_name}.dat').split()
    return ref_energy