from shutil import copyfile
from shutil import copyfileobj
from shutil import rmtree
from pathlib import Path
import os
import sys

def read_text(file_path: Path, mode:str = 'r') -> list[str]:
    with open(file_path, mode) as f:
        data = f.readlines()
    return data

def write_text(contents: str, file_path: Path, mode: str = 'w') -> None:
    """
    Writes text from """
    with open(file_path, mode) as f:
        f.write(contents)

def copy_file(src: Path, dst: Path) -> None:
    copyfile(src, dst)

def check_dir(path: Path) -> Path:
    if path.exists():
        if path.is_file():
            path.unlink()
        else:
            try:
                rmtree(path)
            except PermissionError as exc:
                raise PermissionError(
                    f"Cannot remove existing directory {path}. It may be locked or in use."
                ) from exc
    path.mkdir(parents=True, exist_ok=True)
    return path

def remove_files(path: Path, file_list: list[str]):
    for file in file_list:
        file_for_removal = Path(path / file)
        file_for_removal.unlink()

def check_files(path: Path, file_required: list[str]):
    dnt_exist = [f for f in file_required if not os.path.exists(path / f)]
    if dnt_exist:
        for i in dnt_exist:
            print(f"Missing File: {i}")
        raise SystemExit(f"Please put required files into {path.name}")
    print("All required files present...")