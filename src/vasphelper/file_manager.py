from shutil import copyfile
from shutil import copyfileobj
from shutil import rmtree
from pathlib import Path

def read_text(file_path: Path, mode:str = 'r') -> list[str]:
    with open(file_path, mode) as f:
        data = f.readlines()
    return data

def write_text(contents: str, file_path: Path, mode: str = 'w') -> None:
    with open(file_path, mode) as f:
        f.write(contents)

def copy_file(src: Path, dst: Path) -> None:
    copyfile(src, dst)

def check_dir(path: Path) -> Path:
    # only needs to be done for outer most directory
    if path.exists():
        rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return(path)

def remove_files(path: Path, file_list: list[str]):
    for file in file_list:
        file_for_removal = Path(path / file)
        file_for_removal.unlink()
    print("Removing temporary files..")


