import os
from chimerax.core.commands import run

# Change to the folder with data files
path = '/Users/davidseiferth/Documents/PhD/PDB_PoreFinding/chimera_script' #'CHANGE_THIS_TO_YOUR_PATH'
fname = '7tu9_aligned_z.pdb' #'YOUR_INPUT_FILE_aligned_z.pdb'
os.chdir(path)
file_names = [fname + '_ellipsoid.pdb', fname]

# loop through the files, opening, processing, and closing each in turn
i = 0
for fn in file_names:
    print(f"Processing {fn}")  # show what file we're working on
    run(session, f"open {fn}")
    if i == 0:
        # Modify van der Waals radius using 'setattr' command in ChimeraX
        run(session, 'surface')
        run(session, 'setattr a radius +0.1')  # Adjust atomic radius
        i += 1
    else:
        run(session, 'cartoon')
