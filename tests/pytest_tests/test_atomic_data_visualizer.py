from pathlib import Path
import numpy as np
import pytest

from vasphelper import atomic_data_visualizer as adv
from vasphelper import math_functions as mf


def test_read_bader_file_parses_text(monkeypatch):
    sample_lines = [
        "header line 1\n",
        "header line 2\n",
        "1 C 0 0 0.123\n",
        "2 O 0 0 0.456\n",
        "3 H 0 0 0.789\n",
        "footer1\n",
        "footer2\n",
        "footer3\n",
        "footer4\n",
    ]
    monkeypatch.setattr(adv.fm, "read_text", lambda filepath: sample_lines)

    arr = adv.read_bader_file(Path("ACF_test.dat"), "test_case", ["C", "O", "H"])

    assert arr.dtype.names == ("Number", "Case", "Element", "Data")
    assert arr["Case"].tolist() == ["test_case", "test_case", "test_case"]
    assert arr["Element"].tolist() == ["C", "O", "H"]
    assert np.allclose(arr["Data"], [0.123, 0.456, 0.789])


def test_read_outfile_uses_columns(monkeypatch):
    sample_lines = [
        "A 1 2.0 3.5\n",
        "B 2 4.0 6.5\n",
    ]
    monkeypatch.setattr(adv.fm, "read_text", lambda filepath: sample_lines)

    dtype = [
        ("Case", "U20"),
        ("Idx", "i4"),
        ("Label", "U1"),
        ("X", "f4"),
        ("Y", "f4"),
    ]
    arr = adv.read_outfile(Path("dummy.dat"), "case0", ["A", "B"], (2, 4), dtype)

    assert arr.shape == (2,)
    assert arr["Case"].tolist() == ["A", "B"]
    assert arr["Idx"].tolist() == [1, 2]
    assert arr["Label"].tolist() == ["1", "2"]
    assert np.allclose(arr["X"], [2.0, 4.0])
    assert np.allclose(arr["Y"], [3.5, 6.5])


def test_get_case_list(tmp_path, monkeypatch):
    monkeypatch.setattr(adv, "CUR_DIR", tmp_path)
    (tmp_path / "ACF_a.dat").write_text("x")
    (tmp_path / "ACF_b.dat").write_text("x")
    (tmp_path / "DOSCAR_c").write_text("x")

    case_list = adv.get_case_list("ACF_", ".dat")
    assert sorted(case_list) == ["a", "b"]


def test_handle_clbes_builds_structured_array(monkeypatch):
    monkeypatch.setattr(
        mf,
        "calculate_clbes",
        lambda filepath, surf_energy, ref_energy, exp_energy: [
            ("C", "1", -783.2),
            ("O", "2", -784.1),
        ],
    )

    class FakeContcar:
        def __init__(self):
            self.types_nums = {"C": (None, 1), "O": (None, 2)}

    contcar = FakeContcar()

    arr = adv.handle_clbes(
        contcar, #type: ignore 
        0.1,
        "mycase",
        ["C", "O", "O"],
        -646.92,
        -137.8,
        500,
    )

    assert arr.dtype.names == ("Case", "Number", "Element", "Data")
    assert arr["Case"].tolist() == ["mycase", "mycase", "mycase"]
    assert arr["Element"].tolist() == ["C", "O", "O"]
    assert arr["Number"].tolist() == [1, 2, 3]
    assert np.allclose(arr["Data"], [-783.2, -5.0, -784.1])


def test_get_values_by_type_repeats_types():
    assert adv.get_values_by_type([2, 1, 3], ["C", "O", "H"]) == [
        "C",
        "C",
        "O",
        "H",
        "H",
        "H",
    ]


def test_inter_to_zero_integrates_positive_area():
    x = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
    y = np.array([0.5, 0.75, 1.0, 0.5, 0.0], dtype=float)
    value = mf.inter_to_zero(x, y, num_pts=100)
    assert value == pytest.approx(0.75, rel=1e-6)


def test_calculate_orbitals_returns_per_atom_totals():
    sample_doscar = [
        "1 1 1 0\n",
        "0.2132113E+02  0.1100190E-08  0.1107400E-08  0.2800000E-08  0.5000000E-15\n",
        "1.0\n",
        "CAR\n",
        "EM-structure-relaxtation\n",
        "1.037    -0.482    12    0.0      1.00000000\n",
        "-0.482  0.1254E+03  0.8561E+03\n",
        "-0.344  0.1296E+03  0.8740E+03\n",
        "-0.206  0.1365E+03  0.8929E+03\n",
        "-0.068  0.1440E+03  0.9128E+03\n",
        "0.070  0.1475E+03  0.9331E+03\n",
        "0.208  0.1437E+03  0.9530E+03\n",
        "0.346  0.1330E+03  0.9714E+03\n",
        "0.484  0.1180E+03  0.9877E+03\n",
        "0.623  0.9881E+02  0.1001E+04\n",
        "0.761  0.7482E+02  0.1012E+04\n",
        "0.899  0.4862E+02  0.1018E+04\n",
        "1.037  0.2608E+02  0.1022E+04\n",
        "1.037    -0.482  12  0.00000 1.0\n",                    
        "-0.482  0.1252E-01  0.4555E-01  0.1160E+00  0.6922E-01\n",
        "-0.344  0.1140E-01  0.4858E-01  0.1219E+00  0.8199E-01\n",
        "-0.206  0.1110E-01  0.5459E-01  0.1354E+00  0.9669E-01\n",
        "-0.068  0.1182E-01  0.6286E-01  0.1518E+00  0.1107E+00\n",
        "0.070  0.1293E-01  0.7173E-01  0.1631E+00  0.1191E+00\n",
        "0.208  0.1339E-01  0.7973E-01  0.1614E+00  0.1169E+00\n",
        "0.346  0.1267E-01  0.8684E-01  0.1451E+00  0.1038E+00\n",
        "0.484  0.1117E-01  0.9158E-01  0.1201E+00  0.8542E-01\n",
        "0.623  0.9290E-02  0.8794E-01  0.9221E-01  0.6596E-01\n",
        "0.761  0.6898E-02  0.7081E-01  0.6364E-01  0.4612E-01\n",
        "0.899  0.4209E-02  0.4481E-01  0.3726E-01  0.2730E-01\n",
        "1.037  0.1997E-02  0.2147E-01  0.1762E-01  0.1299E-01\n",
    ]
    totals, overall = mf.calculate_orbitals(sample_doscar[int(12) + 6 :], num_atoms=1, dos_len=4, pts=100)
    assert len(totals) == 1
    assert totals[0][0] == pytest.approx(0.005581677321693207, rel=1e-6)
    assert totals[0][1] == pytest.approx(0.026059419011388415, rel=1e-6)
    assert totals[0][2] == pytest.approx(0.0643933379150558, rel=1e-6)
    assert overall == pytest.approx(0.14081727311442724, rel=1e-6)


def test_handle_pdos_builds_structured_array(monkeypatch):
    sample_doscar = [
        "1 1 1 0\n",
        "0.2132113E+02  0.1100190E-08  0.1107400E-08  0.2800000E-08  0.5000000E-15\n",
        "1.0\n",
        "CAR\n",
        "EM-structure-relaxtation\n",
        "1.037    -0.482    12    0.0      1.00000000\n",
        "-0.482  0.1254E+03  0.8561E+03\n",
        "-0.344  0.1296E+03  0.8740E+03\n",
        "-0.206  0.1365E+03  0.8929E+03\n",
        "-0.068  0.1440E+03  0.9128E+03\n",
        "0.070  0.1475E+03  0.9331E+03\n",
        "0.208  0.1437E+03  0.9530E+03\n",
        "0.346  0.1330E+03  0.9714E+03\n",
        "0.484  0.1180E+03  0.9877E+03\n",
        "0.623  0.9881E+02  0.1001E+04\n",
        "0.761  0.7482E+02  0.1012E+04\n",
        "0.899  0.4862E+02  0.1018E+04\n",
        "1.037  0.2608E+02  0.1022E+04\n",
        "1.037    -0.482  12  0.00000 1.0\n",                    
        "-0.482  0.1252E-01  0.4555E-01  0.1160E+00  0.6922E-01\n",
        "-0.344  0.1140E-01  0.4858E-01  0.1219E+00  0.8199E-01\n",
        "-0.206  0.1110E-01  0.5459E-01  0.1354E+00  0.9669E-01\n",
        "-0.068  0.1182E-01  0.6286E-01  0.1518E+00  0.1107E+00\n",
        "0.070  0.1293E-01  0.7173E-01  0.1631E+00  0.1191E+00\n",
        "0.208  0.1339E-01  0.7973E-01  0.1614E+00  0.1169E+00\n",
        "0.346  0.1267E-01  0.8684E-01  0.1451E+00  0.1038E+00\n",
        "0.484  0.1117E-01  0.9158E-01  0.1201E+00  0.8542E-01\n",
        "0.623  0.9290E-02  0.8794E-01  0.9221E-01  0.6596E-01\n",
        "0.761  0.6898E-02  0.7081E-01  0.6364E-01  0.4612E-01\n",
        "0.899  0.4209E-02  0.4481E-01  0.3726E-01  0.2730E-01\n",
        "1.037  0.1997E-02  0.2147E-01  0.1762E-01  0.1299E-01\n",
    ]
    monkeypatch.setattr(adv.fm, "read_text", lambda filepath: sample_doscar)

    dummy_contcar = object()
    arr = adv.handle_pdos(dummy_contcar, 0.1, "testcase", ["C"]) #type: ignore
    print(arr['s'])
    print(arr['p'])
    print(arr['t'])
    assert arr.dtype.names == ("Case", "Number", "Element", "s", "p", "d", "f", "t")
    assert arr["Case"].tolist() == ["testcase"]
    assert arr["Element"].tolist() == ["C"]
    assert np.allclose(arr["s"], [0.00557865])
    assert np.allclose(arr["p"], [0.02605577])
    assert np.allclose(arr["t"], [0.14078967])


def test_get_max_value_and_make_increments():
    dtype = [("Case", "U10"), ("Number", "i4"), ("Element", "U2"), ("Data", "f4")]
    values = np.array(
        [("a", 1, "C", 2.1), ("a", 2, "O", 1.4), ("a", 3, "C", 0.9)],
        dtype=dtype,
    )
    
    max_map = mf.get_max_value(values, "Data", 'Element')
    assert max_map == {"C": 2.1, "O": 1.4}

    increments = mf.make_increments(0.5, values, "Data", 'Element', len(adv.COLOR_LIST))
    assert set(increments) == {"C", "O"}
    assert len(increments["C"]) == len(adv.COLOR_LIST)
    assert increments["C"][-1] == pytest.approx(2.1, rel=1e-6)
    assert increments["O"][-1] == pytest.approx(1.4, rel=1e-6)


def test_color_atoms_by_value_assigns_neutral_for_missing_data():
    dtype = [("Case", "U10"), ("Number", "i4"), ("Element", "U2"), ("Data", "f4")]
    values = np.array(
        [("a", 1, "C", 1.0), ("a", 2, "O", -5.0)],
        dtype=dtype,
    )
    increments = {"C": [0.0, 5.0], "O": [0.0, 5.0]}
    color_info = [["", "", "", "220", "220", "220"], ["", "", "", "220", "220", "220"]]

    colored = adv.color_atoms_by_value(values, "Data", increments, color_info, ["C", "O"])
    assert colored[0][3:6] != adv.COLOR_NEUTRAL
    assert colored[1][3:6] == adv.COLOR_NEUTRAL