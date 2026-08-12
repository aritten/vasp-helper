import pandas as pd
import pytest
from pathlib import Path

from vasphelper import icore_input_maker as iim


class DummyContcar:
    def __init__(self):
        self.coordinates = {
            0: (0.0, 0.0, 0.0),
            1: (0.9, 0.0, 0.0),
            2: (1.1, 0.0, 0.0),
        }
        self.box_dim = (1.0, 1.0, 1.0)
        self.all_atoms = {
            0: ("O", "", 1),
            1: ("H", "", 2),
            2: ("H", "", 3),
        }
        self.types_nums = {"O": 1, "H": 2}

    def clean_xyz_data(self):
        pass

    def find_all_atom(self):
        pass

    def find_box_dim(self):
        pass


def test_sp_pbc_dist_wraps_coordinates():
    point = (0.0, 0.0, 0.0)
    neighbors = [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0), (0.0, 0.9, 0.0)]
    distances = iim.sp_pbc_dist(point, neighbors, (1.0, 1.0, 1.0))
    assert distances == pytest.approx([0.25, 0.25, 0.01])


def test_find_atoms_closest_selects_one_per_type():
    contcar = DummyContcar()
    result = iim.find_atoms_closest(contcar, aoi=0, num_surr_atoms=1) #type: ignore

    assert len(result) == 2
    assert 0 in result
    assert 1 in result
    assert result[0] == ("O", "", 1)
    assert result[1] == ("H", "", 2)


def test_read_icore_data_reads_csv(tmp_path, monkeypatch):
    csv_file = tmp_path / "corelevel.csv"
    csv_file.write_text("Type,CLN,CLL\nC,1,2\n")
    monkeypatch.setattr(iim, "CUR_DIR", tmp_path)
    df = iim.read_icore_data()
    assert df.loc["C", "CLN"] == 1
    assert df.loc["C", "CLL"] == 2


def test_make_icore_directory_creates_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(iim, "ICORE_DIR", tmp_path)
    result = iim.make_icore_directory(("C", "", 1))
    assert result == tmp_path / "C1"
    assert result.exists()
    assert result.is_dir()


def test_handle_ads_icore_calls_partial_handler(monkeypatch):
    recorded = {}

    def fake_find_atoms_closest(contcar, aoi, num_surr_atoms):
        recorded["called"] = True
        recorded["args"] = (aoi, num_surr_atoms)
        return {0: ("C", "", 1)}

    monkeypatch.setattr(iim, "find_atoms_closest", fake_find_atoms_closest)

    class Dummy:
        pass

    result = iim.handle_ads_icore(Dummy(), {}, partial=True, num_surr_atoms=5, aoi=7) #type: ignore
    assert recorded["called"]
    assert recorded["args"] == (7, 5)
    assert result == {0: ("C", "", 1)}


def test_handle_surf_icore_returns_relax_atoms():
    class DummyContcar2:
        def __init__(self):
            self.relax_atoms = {0: ("O", "", 1)}

        def find_relax_atoms(self):
            self.relax_atoms = {0: ("O", "", 1)}

    contcar = DummyContcar2()
    result = iim.handle_surf_icore(contcar, {}) #type: ignore
    assert result == {0: ("O", "", 1)}


def test_handlge_bulk_icore_returns_all_atoms():
    class DummyContcar3:
        def __init__(self):
            self.all_atoms = {0: ("Fe", "", 1)}

        def find_all_atom(self):
            self.all_atoms = {0: ("Fe", "", 1)}

    contcar = DummyContcar3()
    result = iim.handle_bulk_icore(contcar, {}) #type: ignore
    assert result == {0: ("Fe", "", 1)}


def test_make_icore_files_invokes_vasp_helpers(tmp_path, monkeypatch):
    monkeypatch.setattr(iim, "ICORE_DIR", tmp_path)
    monkeypatch.setattr(
        iim,
        "read_icore_data",
        lambda: pd.DataFrame({"CLN": {"C": "1"}, "CLL": {"C": "2"}}),
    )

    written = []
    populated = []
    removed = []
    calc_written = []

    def fake_write_text(text, dest):
        written.append((text, Path(dest)))

    def fake_populate_vasp_dirs(root, contcar_path, directory, types, params):
        populated.append((Path(root), Path(contcar_path), Path(directory), types.copy(), params.copy()))

    def fake_remove_files(directory, files):
        removed.append((Path(directory), files.copy()))

    def fake_write_calcfile(directory, dirs):
        calc_written.append((Path(directory), [Path(x) for x in dirs]))

    monkeypatch.setattr(iim.fm, "write_text", fake_write_text)
    monkeypatch.setattr(iim.vfm, "populate_vasp_dirs", fake_populate_vasp_dirs)
    monkeypatch.setattr(iim.fm, "remove_files", fake_remove_files)
    monkeypatch.setattr(iim.vfm, "write_calcfile", fake_write_calcfile)

    class DummyContcar4:
        def __init__(self):
            self.icore_types = ["C"]
            self.icore_index = 1

        def create_icore_contcar(self, atom):
            return "ICORE_CONTCAR"

    contcar = DummyContcar4()
    atom_dict = {0: ("C", "", 1)}
    params = {"IBRION": "-1", "NSW": "0"}

    iim.make_icore_files(contcar, atom_dict, params) #type: ignore

    assert written
    assert populated
    assert removed
    assert calc_written
    assert written[0][1].name == "CONTCAR_C1"
    assert calc_written[0][0] == tmp_path