import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj # for emitting status messages

# change to folder with data files
path = 'CHANGE_THIS_TO_YOUR_PATH'
fname = 'YOUR_INPUT_FILE_aligned_z.pdb'
os.chdir(path)
file_names = [ fname + '_circle.pdb', fname]

# loop through the files, opening, processing, and closing each in turn
i = 0 
for fn in file_names:
	replyobj.status("Processing " + fn)  # show what file we're working on
	rc("open " + fn)
	if i==0:
		#rc('molmap  #0  4  gridSpacing  0.5') # Atom Specification #0&protein 
		rc('vdwdefine  +.1')
		rc("surf")  # surface pore
		i += 1

#rc("turn z 90 90")
rc("preset apply publication 1") # make everything look nice

# save image to a file that ends in .png rather than .pdb
png_name = file_names[1][:-3] + "png"
rc("copy file " + png_name + " supersample 3")

### https://www.cgl.ucsf.edu/chimera/current/docs/ProgrammersGuide/processData.py