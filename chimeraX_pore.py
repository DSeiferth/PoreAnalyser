import os
from chimerax.core.commands import run

# Assuming you have access to the session, use this method:
def run_command(session, command):
    run(session, command)

# Change to the folder with data files
path = 'CHANGE_THIS_TO_YOUR_PATH'
fname = 'YOUR_INPUT_FILE_aligned_z.pdb'
os.chdir(path)
file_names = [fname + '_ellipsoid.pdb', fname]

# Access the current ChimeraX session
session = session  # Assuming this script is being run in an environment where the session is available

# loop through the files, opening, processing, and closing each in turn
i = 0
for fn in file_names:
    print(f"Processing {fn}")  # show what file we're working on
    run_command(session, f"open {fn}")
    if i == 0:
        # Modify van der Waals radius using 'setattr' command in ChimeraX
        run_command(session, 'surface')
        run_command(session, 'setattr a radius +0.1')  # Adjust atomic radius
        i += 1
    else:
        run_command(session, 'cartoon')
