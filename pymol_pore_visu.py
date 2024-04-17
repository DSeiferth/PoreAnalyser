from pymol import cmd                       # PyMOL visualisation commands
from pymol import __script__                # location of this script
from os import path                         # ability to find pathname
import sys                                  # ability to append to import path
sys.path.append(path.dirname(__script__))   # find wobj.py in script directory
import argparse                             # command line argument parsing

parser = argparse.ArgumentParser()
parser.add_argument(
    "-structure",
    nargs = "?",
    const = "7tvi_aligned_z.pdb",
    default = "7tvi_aligned_z.pdb")
parser.add_argument(
    "-surface",
    nargs = "?",
    const = "7tvi_aligned_z.pdb_ellipsoid.pdb",
    default = "7tvi_aligned_z.pdb_ellipsoid.pdb")

args = parser.parse_args()

# load structure into object named structure:
cmd.load(args.structure, "structure")
cmd.hide("all")
cmd.show("cartoon", 'structure')

# load pore surface
cmd.load(args.surface, "pore")
cmd.show('surface', 'pore')

