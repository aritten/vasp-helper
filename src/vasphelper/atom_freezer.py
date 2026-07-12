#!/usr/bin/python3

import sys
import argparse
from vasphelper import vasp_file_manager as vfm
from vasphelper import file_manager as fm
from pathlib import Path

CUR_DIR: Path = Path.cwd()

def run_freeze_atoms(filename: Path, freeze_type: str, *, zpos: float | None = None, num_layer: int | None = None, relaxed_layers: int | None = None, num_ads: int = 0, tolerance: float = 0) -> None:
    
    file_path = CUR_DIR / filename

    contcar: vfm.ContcarClass = vfm.ContcarClass(file_path, num_ads)
    contcar.parse_atomic_data()
    contcar.clean_xyz_data()
    if freeze_type == 'layer':
        split = contcar.freezer_by_layer(num_layer, relaxed_layers, tolerance)
        pos_content = contcar.write_out_select_tags(split)
    else:
        pos_content = contcar.write_out_select_tags(zpos)


    fm.write_text(pos_content, CUR_DIR / 'POSCAR')
    
def main():
    #add a way to just set zpos at which to freeze atoms
    parser = argparse.ArgumentParser(description=f"""This program is used to freeze atoms by layers.
In order to run this program, you need:
{'-'*60}
- CONTCAR that you wish to freeze
{'-'*60}
""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('filename', help='Filename to of VASP position file to freeze atoms within.')
    subparsers = parser.add_subparsers(dest="freeze_type", required=True,help="Choose whether to freeze atoms based on layer or z-position")
    layer_parser = subparsers.add_parser('layer', help='Freeze atoms based on layer position')
    layer_parser.add_argument('num_layers', help='Number of layers in the position file.', type=int)
    layer_parser.add_argument('relaxed_layers', help='Number of relaxed layers needed', type=int)
    layer_parser.add_argument('-n', '--num_ads', help='Specifies if adspecies are present.', type=int)
    layer_parser.add_argument('-t', '--tolerance', help='Adds tolerance to shift the layer split up or down by.', type=float, default=0.01)
    zpos_parser = subparsers.add_parser('zpos', help='Freeze based on z-positon')
    zpos_parser.add_argument('zpos', help='The z-position that all atoms will be relaxed above.', type=float)
    args = parser.parse_args()
    
    if args.freeze_type == 'layer':
        run_freeze_atoms(
            filename=args.filename,
            freeze_type=args.freeze_type,
            num_layer=args.num_layers,
            relaxed_layers=args.relaxed_layers,
            num_ads=args.num_ads,
            tolerance=args.tolerance
        )
    else:
        run_freeze_atoms(
            filename=args.filename,
            freeze_type=args.freeze_type,
            zpos=args.zpos
            )


if __name__ == '__main__':
    main()