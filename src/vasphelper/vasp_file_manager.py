
from shutil import copyfileobj
from pathlib import Path
from vasphelper import file_manager as fm
from itertools import groupby
from typing import Any
from itertools import accumulate
import sys

class ContcarClass:
     
    def __init__(self, path: Path, num_ads: int) -> None:
        self.path: Path = path
        self.num_ads: int = num_ads
        self.raw_content: list[str] = [line.strip('\n') for line in fm.read_text(path)]
        self.content: list[str] = self.remove_velocity_data()
        self.atomic_data, self.xyz = self.split_atomic_data_and_ionic_pos()

    def remove_velocity_data(self) -> list[str]:
        """
        Removes velocity data from CONTCAR retaining only data related to lattice parameters, atomic information and ionic positions.
        The delimiter is set at '' as is current convention in VASP for output CONTCARs

        Args:
            None
        Returns:
            contcar_contents (list[str]): contains all data relating to lattice parameters, atomic information and ionic positions
        """
        delimiter: list[str] = [' ', '']
        return [list(group) for k, group in groupby(self.raw_content, lambda x: x in delimiter) if not k][0]
    
    def split_atomic_data_and_ionic_pos(self) -> tuple[list[str], list[str]]:
        """
        Splits atomic data and ionic positions of a CONTCAR. Delimiter line number is dependent on the use of selective dynamics. "9" if selective dynamics and "8" if not.
        
        Args:
            None
        Returns:
            header (list[str]): contains list containing header from CONTCAR
            ionic positions (list[str]): contains all ionic positions of atoms within the surface
        """
        line_delimiter: int = 0
        if self.content[7].lower().startswith('s'):
            line_delimiter = 9
        else:
            line_delimiter = 8
        
        return self.content[:line_delimiter], self.content[line_delimiter:]

    def find_box_dim(self) -> None:
        if self.atomic_data[-1].lower().startswith('d'):
            self.box_dim = [1, 1, 1]
        else:
            print("Please use fractional coordinates.")
            sys.exit()
    
    def write_out_select_tags(self, zpos):
        for index, val in enumerate(self.coordinates):
            coords = '  ' + '  '.join([f'{entry:.16f}' for entry in self.coordinates[index]])
            if val[2] > zpos:
                self.xyz[index] = coords + '   T   T   T'
            else:
                self.xyz[index] = coords + '   F   F   F'
        
        if not self.atomic_data[7].lower().startswith('s'):
            self.atomic_data.insert(-1, 'Selective Dynamics')
        up_pos_file = self.atomic_data + self.xyz

        return '\n'.join(up_pos_file)


    def freezer_by_layer(self, num_layers, relaxed_layers, tolerance):
        surf_total: int = sum(self.nums[:self.num_ads + 1])
        surf_atom_coords = self.coordinates[:surf_total]
        num_surf_atom_coords = zip(range(1,surf_total+1), surf_atom_coords)

        sort_by_z_pos = sorted(num_surf_atom_coords, key= lambda item: item[1][2])
        sublayers = [0.0]
        for i in range(1,len(sort_by_z_pos)):
            diff = sort_by_z_pos[i][1][2] - sort_by_z_pos[i-1][1][2]
            if diff > tolerance:
                sublayers.append(diff/2 + sort_by_z_pos[i-1][1][2])
        
        sub_per_layer = len(sublayers) / num_layers
        layers = sublayers[::int(sub_per_layer)][::-1]
        split = layers[relaxed_layers-1]
        print(f'Atoms relaxed above {split}')
        return split

    def clean_xyz_data(self) -> None:
        """
        Removes True and False tags for freezing atoms within the unit cell.
        Args:
            None
        Side Effects:
            Adds coordinates attribute to ContcarClass.
        Returns:
            None
        """
        self.coordinates: list[list[float]] = []
        for i, line in enumerate(self.xyz):
            self.coordinates.append(list(map(float, line.split()[0:3])))

    def parse_atomic_data(self) -> None:
        """
        Parses atomic types and numbers of each atom type as well as a comined atom dictionary that contains atoms that are repeated in both the surface and 
        Args:
            None
        Side Effects:
            Adds types, nums and types_nums to the class attributes of ContcarClass.
            types (list[str]): contains a list of all atom types of atoms in the surface
            nums (list[str]): contains number of all 
        Returns:
            None
        """
        self.types: list[str] = self.atomic_data[5].split()
        self.nums: list[int] = list(map(int, self.atomic_data[6].split()))
        combined: list[tuple[str, int]] = list(zip(self.types, self.nums))
        if self.num_ads != 0:
            types_ads: list[str] = []
            num_surf_spec = len(self.types) - self.num_ads
            for i, entry in enumerate(self.types):
                if i >= num_surf_spec:
                    types_ads.append(f'{entry}_ads')
                else:
                    types_ads.append(entry)
            self.types_nums: dict[str, tuple[str, int]] = dict(zip(types_ads, combined))
        else: 
            self.types_nums: dict[str, tuple[str, int]] = dict(zip(self.types, combined))

    def diff_types_nums(self, *, species_start: int = 0, species_end: int | None = None) -> tuple[str, str]:
        num_line = '    ' + '   '.join(str(i) for i in self.nums[species_start:species_end])
        types_line = '   '.join(self.types[species_start:species_end])
        return num_line, types_line
    
    def icore_types_nums(self, atom: tuple[str, str, int], types: list[str], nums: list[int], index: int, last: int) -> tuple[str, str]:
        element, tag, pos = atom
        count = nums[index]
        line_pos: int = 0
        if pos == 1:
            types[index : index] = [element]
            nums[index : index + 1] = [1, count - 1]
            line_pos = 1 + index
        elif pos == last:
            types[index : index] = [element]
            nums[index : index + 1] = [count - 1, 1]
            line_pos = 2 + index
        else:
            types[index : index] = [element, element]
            greater = count - pos
            less = count - greater - 1
            nums[index: index + 1] = [less, 1, greater]
            line_pos = 2 + index
        
        self.icore_types = types
        self.icore_index = line_pos
        types_line = '   '.join(types)
        nums_line = '    ' + '   '.join(str(i) for i in nums)
        return types_line, nums_line
    
    def create_split_contcar(self, name: str) -> str:
        total_species: int  = len(self.types)
        surf_total: int = sum(self.nums[:self.num_ads + 1])
        num_surf: int = total_species - self.num_ads
        if name == 'surf':
            xyz = self.xyz[:surf_total]
            nums_line, types_line = self.diff_types_nums(species_end=(num_surf))
        else:
            xyz = self.xyz[surf_total:]
            nums_line, types_line = self.diff_types_nums(species_start=(num_surf))
        return '\n'.join(['\n'.join(self.atomic_data[:5]), types_line, nums_line, self.atomic_data[7], self.atomic_data[8], '\n'.join(xyz)])
    
    def create_icore_contcar(self, atom: tuple[str, str, int]) -> str:
        names: list[str] = list(self.types_nums.keys())
        types: list[str] = list(self.types)
        nums: list[int] = list(self.nums)
        name = atom[0] + atom[1]
        index = names.index(name)
        last = self.types_nums[name][1]
        types_line, nums_line = self.icore_types_nums(atom, types, nums, index, last)
        return '\n'.join(['\n'.join(self.atomic_data[:5]), types_line, nums_line, self.atomic_data[7], self.atomic_data[8], '\n'.join(self.xyz)])
    
    def find_relax_atoms(self) -> None:
        
        total_atoms: int = sum(self.nums)
        ads_atoms: int = sum(self.nums[-self.num_ads:])
        surf_total = total_atoms - ads_atoms - 1

        self.relax_atoms: dict[int, tuple[str, str, int]] = {}
        data = tuple(zip(self.types, list(accumulate(self.nums))))
        atom_index = 0
        prev_end = 0
    
        for atom, coordinate in enumerate(self.xyz):
            # this could be changed to provide context to xyz data
            if coordinate.endswith('T'):
                while atom >= data[atom_index][1]:
                    prev_end = data[atom_index][1]
                    atom_index += 1
            if atom > surf_total:
                self.relax_atoms[atom] = (data[atom_index][0], '_ads' ,atom - prev_end + 1)
            else:
                self.relax_atoms[atom] = (data[atom_index][0], '' ,atom - prev_end + 1)

    def find_all_atom(self) -> None:

        total_atoms: int = sum(self.nums)
        ads_atoms: int = sum(self.nums[-self.num_ads:])
        surf_total = total_atoms - ads_atoms - 1

        self.all_atoms: dict[int, tuple[str, str, int]] = {}
        data = tuple(zip(self.types, list(accumulate(self.nums))))
        atom_index = 0
        prev_end = 0

        for atom, coordinate in enumerate(self.xyz):
            while atom >= data[atom_index][1]:
                prev_end = data[atom_index][1]
                atom_index += 1
            if atom > surf_total:
                self.all_atoms[atom] = (data[atom_index][0], '_ads' ,atom - prev_end + 1)
            else:
                self.all_atoms[atom] = (data[atom_index][0], '' ,atom - prev_end + 1)

def check_unique_atom_atom_types(atom_list: list[str]) -> set:
    unique = set()
    for v in atom_list:
        unique.add(v)
    return unique

def potcar_concatentate(atom_list: list[str], src: Path, dest: Path) -> None:
    with open(dest / 'POTCAR', 'wb') as destination:
        for atom in atom_list:
            with open(src / f'POTCAR_{atom}', 'rb') as source:
                copyfileobj(source, destination)

def change_incar_parameters(source: Path, dest: Path, parameter_dict: dict[str, str]) -> None:
# change this but works
    content: list[str] = fm.read_text(source / 'INCAR')
    if not content[-1].endswith('\n'):
        content[-1] += '\n'
    
    exists = set()

    for i, line in enumerate(content):
            tag = line.split('=', 1)[0].strip()
            if tag in parameter_dict:
                content[i] = f'{tag:<9}= {parameter_dict[tag]}\n'
                exists.add(tag)

    for key in parameter_dict:
        if key not in exists:
            content.append(f'{key:<9}= {parameter_dict[key]}\n')
    
    fm.write_text("".join(content), dest / 'INCAR')

# should this be else where?
def copy_contcars_diff(num_ads: int, data: list[str], dir_name: Path, cur_dir: Path) -> dict[str, list[str]]:
    atom_dict: dict[str, list[str]] = {}

    print(f'Updating CONTCARs for {dir_name.name}...')

    for entry in data:
        contcar: ContcarClass = ContcarClass(cur_dir / f'CONTCAR_{entry}', num_ads)
        contcar.parse_atomic_data()
        if dir_name.name == 'overall':
            fm.copy_file(cur_dir / f'CONTCAR_{entry}', dir_name / f'CONTCAR_{entry}')
        else:
            split_contcar = contcar.create_split_contcar(dir_name.name)
            fm.write_text(split_contcar, dir_name / f'CONTCAR_{entry}')

        atom_dict[entry] = contcar.types

    return atom_dict

def populate_vasp_dirs(cur_dir: Path, contcar_path: Path, directory: Path, atom_types: list[str], parameter_dict: dict[str, str]) -> None:
    change_incar_parameters(cur_dir, directory, parameter_dict)
    #its possible that this should be done only once for each variation
    potcar_concatentate(atom_types, cur_dir, directory)
    fm.copy_file(cur_dir / 'KPOINTS', directory / 'KPOINTS')
    #build in CONTCAR copying and changing into this
    fm.copy_file(contcar_path, directory / 'POSCAR')

def write_calcfile(dir_path: Path, dir_list: list[Path]):
    with open(dir_path / 'CalcFile.dat', 'w') as f:
        total = 0
        for dir in dir_list:
            dir = str(dir)[len(str(dir_path)) + 1:]
            f.write(f"{dir}\n")
            total = total + 1
    print(f'Total Directories made for {dir_path.name}:  {total}\nDirectory names entered into the CalcFile.dat')