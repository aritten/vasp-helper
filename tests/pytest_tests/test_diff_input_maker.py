from pathlib import Path

import pandas as pd
import pytest

from vasphelper import diff_input_maker as dim

def test_create_list_of_cases_filters_contcar_files(tmp_path, monkeypatch):
    monkeypatch.setattr(dim, "CUR_DIR", tmp_path)

    (tmp_path / "CONTCAR_a").write_text("")
    (tmp_path / "CONTCAR_b").write_text("")
    (tmp_path / "ref_CONTCAR_CO2").write_text("")
    (tmp_path / "ref_CONTCAR_CO2_01").write_text("")
    (tmp_path / "CONTCAR_extra.vasp").write_text("")
    (tmp_path / "README.txt").write_text("")

    result = dim.create_list_of_cases()
    assert set(result) == {"a", "b"}


def test_build_dirs_creates_case_directories(tmp_path):
    case_list = ["a", "b", "c"]
    directories = dim.build_dirs(case_list, tmp_path)

    assert [p.name for p in directories] == case_list
    for directory in directories:
        assert directory.exists()
        assert directory.is_dir()


def test_build_rwigs_dirs_creates_nested_directories(tmp_path, monkeypatch):
    monkeypatch.setattr(dim, "CUR_DIR", tmp_path)
    case_list = ["a", "b"]
    rwigs_data = ["Ratio_1_0", "Ratio_2_0"]
    directories = dim.build_rwigs_dirs(case_list, tmp_path, rwigs_data)

    expected_paths = {str(tmp_path / ratio / case) for ratio in rwigs_data for case in case_list}
    assert {str(p) for p in directories} == expected_paths
    for directory in directories:
        assert directory.exists()
        assert directory.is_dir()


def test_get_rwigs_list_reads_existing_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(dim, "CUR_DIR", tmp_path)
    csv_text = "Type,Ratio_1_0,Ratio_2_0\nCu,0.10,0.20\nO,0.05,0.10\n"
    (tmp_path / "RWIGS_inputs.csv").write_text(csv_text)

    df = dim.get_rwigs_list({"Cu", "O"})
    assert list(df.index) == ["Cu", "O"]
    assert list(df.columns) == ["Ratio_1_0", "Ratio_2_0"]
    assert df.loc["Cu", "Ratio_1_0"] == 0.10
    assert df.loc["O", "Ratio_2_0"] == 0.10


def test_handle_pdos_creates_rwigs_directories_and_calls_populate(tmp_path, monkeypatch):
    monkeypatch.setattr(dim, "CUR_DIR", tmp_path)

    dummy_df = pd.DataFrame(
        {"Ratio_1_0": [0.10, 0.20], "Ratio_2_0": [0.15, 0.25]},
        index=["Cu", "O"],
    )
    dummy_df.index.name = "Type"
    monkeypatch.setattr(dim, "get_rwigs_list", lambda unique: dummy_df)

    calls = []

    def fake_populate_vasp_dirs(src_dir, contcar_path, directory, atoms, incar_dict):
        calls.append(
            {
                "src_dir": src_dir,
                "contcar_path": contcar_path,
                "directory": directory,
                "atoms": atoms,
                "incar": dict(incar_dict),
            }
        )

    monkeypatch.setattr(dim.vfm, "populate_vasp_dirs", fake_populate_vasp_dirs)

    case_list = ["case1"]
    dir_path = tmp_path / "pdos_input"
    dir_path.mkdir()

    (dir_path / "CONTCAR_case1").write_text("")

    incar_parameters = {"IBRION": "-1", "NSW": "0"}
    atom_list = {"case1": ["Cu", "O"]}
    unique_atoms = {"Cu", "O"}

    result_dirs = dim.handle_pdos(
        dir_path,
        case_list=case_list,
        incar_parameter_dict=incar_parameters,
        atom_list=atom_list,
        unique_atoms=unique_atoms,
    )

    expected_dirs = {
        dir_path / "Ratio_1_0" / "case1",
        dir_path / "Ratio_2_0" / "case1",
    }
    assert expected_dirs.issubset(set(result_dirs))
    assert all(p.exists() for p in expected_dirs)
    assert len(calls) == 2
    assert calls[0]["contcar_path"] == dir_path / "CONTCAR_case1"
    assert "RWIGS" in calls[0]["incar"]


def test_handle_bader_creates_directories_and_calls_populate(tmp_path):
    calls = []

    def fake_populate_vasp_dirs(src_dir, contcar_path, directory, atoms, incar_dict):
        calls.append((contcar_path, directory, atoms, dict(incar_dict)))

    dim.vfm.populate_vasp_dirs = fake_populate_vasp_dirs

    case_list = ["case1", "case2"]
    dir_path = tmp_path / "bader_input"
    dir_path.mkdir()
    incar_parameters = {"IBRION": "-1", "NSW": "0"}
    atom_list = {"case1": ["Cu"], "case2": ["O"]}
    unique_atoms = {"Cu", "O"}

    result_dirs = dim.handle_bader(
        tmp_path,
        case_list=case_list,
        incar_parameter_dict=incar_parameters,
        atom_list=atom_list,
        unique_atoms=unique_atoms,
    )

    assert len(result_dirs) == 2
    assert all(p.exists() for p in result_dirs)
    assert calls[0][0] == tmp_path / "CONTCAR_case1"
    assert calls[1][0] == tmp_path / "CONTCAR_case2"
    assert calls[0][2] == ["Cu"]
    assert calls[1][2] == ["O"]
    assert calls[0][1].name == "case1"
    assert "LCHARG" in incar_parameters


def test_handle_chg_creates_directories_and_calls_populate(tmp_path):
    calls = []

    def fake_populate_vasp_dirs(src_dir, contcar_path, directory, atoms, incar_dict):
        calls.append((contcar_path, directory, atoms, dict(incar_dict)))

    dim.vfm.populate_vasp_dirs = fake_populate_vasp_dirs

    case_list = ["case1"]
    dir_path = tmp_path / "chg_input"
    dir_path.mkdir()
    incar_parameters = {"IBRION": "-1", "NSW": "0"}
    atom_list = {"case1": ["Cu", "O"]}
    unique_atoms = {"Cu", "O"}

    result_dirs = dim.handle_chg(
        tmp_path,
        case_list=case_list,
        incar_parameter_dict=incar_parameters,
        atom_list=atom_list,
        unique_atoms=unique_atoms,
    )

    assert len(result_dirs) == 1
    assert result_dirs[0].exists()
    assert calls[0][2] == ["Cu", "O"]
    assert calls[0][0] == tmp_path / "CONTCAR_case1"
    assert "LCHARG" in incar_parameters