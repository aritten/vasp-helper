"""
A module that visualizes data from Bader charge, CLBES and Electron Distribution Analysis.

"""

#!/usr/bin/python3

import argparse
from pathlib import Path
from typing import Any

# CONTSTANTS

CUR_DIR = Path.cwd()

# FUNCTIONS

def handle_pdos():
    pass

def handle_bader():
    ...

def handle_clbes():
    ...

CALC_DISPATCH: dict[str, Any] = {
    'pdos': handle_pdos,
    'bader': handle_bader,
    'clbes': handle_clbes
}

def main():
    parser = argparse.ArgumentParser(description=f"""Split CONTCAR into surface and adsorbant CONTCARs and makes files for one of the follow analyses:
- Bader Charge
- Charge Density
- Electron Distribution Analysis.\n
In order to run this program, you need:
{'-'*60}
- CONTCARs for each case
- POTCARs for each element
- INCAR from geometry calculation
- KPOINTS with correct accuracy for calculation type
{'-'*60}
WARNING:
This program will not change KPOINTS make sure correct KPOINTS are specified in KPOINTS file.
""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("calc_type", choices=['bader', 'clbes', 'pdos'], help="Choose the type of visualization between Bader Charge Analysis 'bader', Core Level Binding Energy Shifts 'clbes' or PDOS Electron Differential Analysis 'pdos'.")
    parser.add_argument("split_type", choices=['surf', 'atom', 'diff'], help="Choose the type of CONTCAR split from the surface 'surf', adsorbate 'ads',both 'both' or 'all which will produce a directory for the surface, adsorbates and overall system.")
    parser.add_argument("width", help="Specifies width of increments.", type=float, default=0.1)
    args = parser.parse_args()

if __name__ == '__main__':
    main()
