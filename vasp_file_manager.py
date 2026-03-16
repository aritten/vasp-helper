
from shutil import copyfileobj
from pathlib import Path
import file_manager as fm


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
            content.append(f'{tag:<9}= {parameter_dict[key]}\n')
    
    fm.write_text("".join(content), dest / 'INCAR')

def create_contcars(num_ads: int, data: list[str], dir_name: Path, cur_dir: Path) -> dict[str, str]:
    
    atom_dic = {}

    for entry in data:

        contcar_cont = [line.strip('\n') for line in fm.read_text(cur_dir / f"CONTCAR_{entry}")]

        atom_types = contcar_cont[5].split()

        atom_nums = contcar_cont[6].split()
        
        atom_nums = list(map(int, atom_nums))
        surf_total: int = sum(atom_nums[0:num_ads + 1])
        total_species: int  = len(atom_types)
        atom_total: int = sum(atom_nums)
        num_line: str = ''
        print(dir_name)
        xyz: list[str] = []
        if dir_name.name == 'surf':
            xyz = contcar_cont[9:(9 + surf_total)]
            num_line = '    ' + '   '.join(str(i) for i in atom_nums[0:(total_species - num_ads)])
            atom_types = atom_types[0:(total_species - num_ads)]
        elif dir_name.name == 'ads':
            xyz = contcar_cont[(9 + surf_total):(9 + atom_total)]
            num_line = '    ' + '   '.join(str(i) for i in atom_nums[(total_species - num_ads):total_species])
            atom_types = atom_types[(total_species - num_ads):total_species]
        elif dir_name.name == 'overall':
            xyz = contcar_cont[9:(9 + atom_total)]
            num_line = '    ' + '   '.join(str(i) for i in atom_nums)

        atom_dic[entry] = atom_types

        #remake CONTCAR
        cont_file_format: str = '\n'.join(['\n'.join(contcar_cont[:5]), "    ".join(atom_types), num_line, contcar_cont[7], contcar_cont[8], '\n'.join(xyz)])
        fm.write_text(cont_file_format, cur_dir / dir_name / f'CONTCAR_{entry}')
    
    print(f'Updating CONTCARs for {dir_name}...')
    return atom_dic

def populate_vasp_dirs(cur_dir: Path, dir_path: Path, directory: Path, atom_types: list[str], parameter_dict: dict[str, str]) -> None:
    change_incar_parameters(cur_dir, directory, parameter_dict)
    potcar_concatentate(atom_types, cur_dir, directory)
    fm.copy_file(cur_dir / 'KPOINTS', directory / 'KPOINTS')
    fm.copy_file(dir_path / f'CONTCAR_{directory.name}', directory / 'POSCAR')


def write_calcfile(dir_path: Path, dir_list: list[Path]) -> int:
    with open(dir_path / "CalcFile.dat", "w") as f:
        total = 0
        for dir in dir_list:
            dir = str(dir)[len(str(dir_path)) + 1:]
            f.write(f"{dir}\n")
            total = total + 1
    return total